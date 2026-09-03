#!/usr/bin/env python3
"""Unit tests for the chrom_copies (constitutional ploidy) model - issue #40.

Covers group B of ploidy_state_table.md plus the guards that keep everything else
diploid. The state table rows are named in each test so a failure points back at the
spec rather than at the code.

The three numbers being tested apart:

    chrom_copies  copies of the region the sample was born with -> allele slots
    locus_copies  copies still present at a locus, ie after a deletion (Zygosity.HEM,
                  Zygosity.NO_COPIES)
    token_copies  how many of those carry this token

Only the first is new. The other two were already in variant_pool in enum form.
"""

import unittest

import pandas as pd

from rbceq2.core_logic.alleles import Allele, BloodGroup
from rbceq2.core_logic.constants import AlleleState, HAPLOID_SECOND_SLOT, PAR
from rbceq2.core_logic.data_procesing import (
    check_token_copies_fit_chrom_copies,
    chrom_copies_for_bg,
    get_fully_homozygous_alleles,
    get_ref,
    get_genotypes,
)
from rbceq2.core_logic.utils import BeyondLogicError, Zygosity
from rbceq2.IO.vcf import PloidyScan
from rbceq2.IO.vcf import VCF, gt_of, is_haploid_gt, is_single_copy


COMMON = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]


def one_row_df(chrom: str, pos: str, GT: str) -> pd.DataFrame:
    """A one-row VCF frame with an array-style FORMAT of just GT."""
    return pd.DataFrame(
        {
            "CHROM": [chrom],
            "POS": [pos],
            "ID": ["."],
            "REF": ["G"],
            "ALT": ["A"],
            "QUAL": ["."],
            "FILTER": ["PASS"],
            "INFO": ["."],
            "FORMAT": ["GT"],
            "SAMPLE": [GT],
        }
    )


class TestParCoordinates(unittest.TestCase):
    """PAR is data entry, so it gets checked against the blood groups that use it."""

    def test_every_db_chrX_gene_lands_on_the_expected_side(self) -> None:
        """XG and CD99 inside PAR1, XK, GATA1 and ATP11C outside, in both builds."""
        cases = [
            # build, pos, single copy?, what it is
            ("GRCh38", 2714453, False, "CD99 start"),
            ("GRCh38", 2717651, False, "CD99 end"),
            ("GRCh38", 2748343, False, "XG start"),
            ("GRCh38", 2776608, False, "XG*01N.03, starts in PAR1 and runs out of it"),
            ("GRCh38", 37680879, True, "XK start"),
            ("GRCh38", 37728261, True, "XK end"),
            ("GRCh38", 48794162, True, "GATA1"),
            ("GRCh38", 139726344, True, "ATP11C"),
            ("GRCh37", 2632494, False, "CD99 start"),
            ("GRCh37", 2694649, False, "XG end"),
            ("GRCh37", 37540132, True, "XK start"),
            ("GRCh37", 48652569, True, "GATA1"),
            ("GRCh37", 138808503, True, "ATP11C"),
        ]
        for build, pos, expected, what in cases:
            with self.subTest(build=build, pos=pos, what=what):
                self.assertEqual(is_single_copy("X", pos, build), expected)

    def test_par_boundaries_are_inclusive(self) -> None:
        for build in ("GRCh37", "GRCh38"):
            for start, end in PAR[build]["X"]:
                with self.subTest(build=build, interval=(start, end)):
                    self.assertFalse(is_single_copy("X", start, build))
                    self.assertFalse(is_single_copy("X", end, build))
                    self.assertTrue(is_single_copy("X", start - 1, build))
                    self.assertTrue(is_single_copy("X", end + 1, build))

    def test_par1_is_not_merged_across_builds(self) -> None:
        """GRCh37 2,699,521-2,781,479 is non-PAR in 37 and PAR in 38.

        A LANE-style merge of the two builds would get this wrong in one direction and
        no db.tsv coordinate currently sits there to catch it.
        """
        self.assertTrue(is_single_copy("X", 2_750_000, "GRCh37"))
        self.assertFalse(is_single_copy("X", 2_750_000, "GRCh38"))

    def test_autosomes_are_never_single_copy(self) -> None:
        self.assertFalse(is_single_copy("1", 25272548, "GRCh38"))  # RHD
        self.assertFalse(is_single_copy("4", 143999443, "GRCh38"))  # GYPB

    def test_unknown_build_is_diploid_rather_than_an_error(self) -> None:
        self.assertFalse(is_single_copy("X", 37686068, "hg19"))


class TestGtHelpers(unittest.TestCase):
    def test_gt_is_taken_from_the_first_format_field(self) -> None:
        self.assertEqual(gt_of("0/1:41,47:88:99"), "0/1")
        self.assertEqual(gt_of("1"), "1")
        self.assertEqual(gt_of(float("nan")), "")

    def test_haploid_detection(self) -> None:
        for GT in ("1", "0", "."):
            self.assertTrue(is_haploid_gt(GT), GT)
        for GT in ("0/1", "1|0", "./.", "1/1", ""):
            self.assertFalse(is_haploid_gt(GT), GT)


class TestInferHaploidChroms(unittest.TestCase):
    """The per-sample half of the answer."""

    def _vcf(self, df, build="GRCh38"):
        return VCF([df], lane_variants={}, unique_variants=set(), sample="s",
                   reference_genome=build)

    def test_haploid_gt_outside_par_marks_the_chromosome(self) -> None:
        """B1."""
        vcf = self._vcf(one_row_df("chrX", "37686068", "1"))
        self.assertEqual(vcf.haploid_chroms, frozenset({"X"}))

    def test_haploid_gt_inside_par_does_not(self) -> None:
        """B5. A haploid call in PAR is a caller error, not evidence of ploidy."""
        vcf = self._vcf(one_row_df("chrX", "2748343", "1"))
        self.assertEqual(vcf.haploid_chroms, frozenset())

    def test_haploid_gt_on_an_autosome_does_not(self) -> None:
        """The Axiom RHD encoding. Locus copy number, not ploidy."""
        vcf = self._vcf(one_row_df("chr1", "25272548", "0"))
        self.assertEqual(vcf.haploid_chroms, frozenset())

    def test_diploid_non_par_sample_is_not_marked(self) -> None:
        """B3. A female, or a ploidy-unaware caller - indistinguishable here."""
        vcf = self._vcf(one_row_df("chrX", "37686068", "1/1"))
        self.assertEqual(vcf.haploid_chroms, frozenset())

    def test_no_reference_genome_disables_inference(self) -> None:
        """Every pre-v2.4.4 construction of VCF keeps its old behaviour."""
        vcf = VCF([one_row_df("chrX", "37686068", "1")], lane_variants={},
                  unique_variants=set(), sample="s")
        self.assertEqual(vcf.haploid_chroms, frozenset())

    def test_a_haploid_no_call_is_not_evidence(self) -> None:
        """E4, resolved against the reading this test used to assert.

        It used to say "'.' is one slot that was not measured, not two slots" - the
        arity of the genotype carrying the ploidy even where the content is missing.
        Defensible, and the state table left it UNRESOLVED. There is now a real example
        and it goes the other way.

        A female sample in a jointly called cohort had 28,331 called diploid genotypes
        across non-PAR X and six bare '.'. Under the old reading those six decided it:
        the sample was reported single copy across X, and XK, ATP11C and GATA1 each lost
        a chromosome she has. Six degenerate no-calls outvoting twenty-eight thousand
        calls is not positive evidence, which is what this inference is supposed to
        require.

        The distinction that survives is between the two haploid genotypes rather than
        between haploid and diploid: a caller does not write '0' or '1' by accident, so
        those stay decisive on their own, and a bare '.' now says nothing.
        """
        vcf = self._vcf(one_row_df("chrX", "37686068", "."))
        self.assertEqual(vcf.haploid_chroms, frozenset())

    def test_a_haploid_called_allele_is_still_evidence(self) -> None:
        """The other half of E4 - what a caller does not write by accident."""
        for GT in ("0", "1"):
            with self.subTest(GT=GT):
                vcf = self._vcf(one_row_df("chrX", "37686068", GT))
                self.assertEqual(vcf.haploid_chroms, frozenset({"X"}))


class TestRemoveHomRefHaploid(unittest.TestCase):
    def test_haploid_zero_is_dropped_where_the_region_is_single_copy(self) -> None:
        """B2. One chromosome and it is reference - the same as '0/0'."""
        vcf = VCF([one_row_df("chrX", "37686068", "0")], lane_variants={},
                  unique_variants=set(), sample="s", reference_genome="GRCh38")
        self.assertEqual(len(vcf.df), 0)

    def test_haploid_zero_is_dropped_on_an_autosome_too(self) -> None:
        """D3. Reversed when locus_copies arrived.

        It used to be kept so the copy number case would be rejected by name rather than
        buried. Now that a gene reported consistently at one copy is read as one copy,
        there is nothing left to reject: the caller has said the copy it can see is
        reference, so the token has zero copies exactly as '0/0' does.
        """
        vcf = VCF([one_row_df("chr1", "25272548", "0")], lane_variants={},
                  unique_variants=set(), sample="s", reference_genome="GRCh38")
        self.assertEqual(len(vcf.df), 0)

    def test_dropping_the_row_keeps_the_ploidy_evidence(self) -> None:
        """The row goes, the fact that it was haploid does not.

        This is the whole reason _infer_locus_ploidy is a separate pass that runs first.
        Without it, a gene called entirely wildtype at one copy would be indistinguishable
        from a gene nobody typed by the time anything could ask.
        """
        vcf = VCF([one_row_df("chr1", "25272548", "0")], lane_variants={},
                  unique_variants=set(), sample="s", reference_genome="GRCh38")
        self.assertEqual(vcf.haploid_loci.get("1"), frozenset({25272548}))
        self.assertEqual(vcf.diploid_loci.get("1"), None)

    def test_diploid_hom_ref_still_dropped(self) -> None:
        for GT in ("0/0", "0|0"):
            with self.subTest(GT=GT):
                vcf = VCF([one_row_df("chr1", "25272548", GT)], lane_variants={},
                          unique_variants=set(), sample="s", reference_genome="GRCh38")
                self.assertEqual(len(vcf.df), 0)


class TestGetRefHaploid(unittest.TestCase):
    def test_haploid_alt_on_a_single_copy_region_is_hemizygous(self) -> None:
        """B1. One chromosome carrying one copy of the token."""
        self.assertEqual(
            get_ref({"GT": "1"}, "X:37686068_G_A", chrom_copies=1), Zygosity.HEM
        )

    def test_haploid_no_call_on_a_single_copy_region_is_no_data(self) -> None:
        """E4."""
        self.assertEqual(
            get_ref({"GT": "."}, "X:37686068_G_A", chrom_copies=1), Zygosity.NO_DATA
        )

    def test_haploid_gt_on_a_two_copy_region_is_rejected_by_name(self) -> None:
        """B5 and E1. Named error rather than a guess - table decisions 2 and 4."""
        with self.assertRaises(BeyondLogicError) as ctx:
            get_ref({"GT": "1"}, "1:25272548_G_A", chrom_copies=2)
        self.assertIn("pseudoautosomal", str(ctx.exception))
        self.assertIn("1:25272548_G_A", str(ctx.exception))

    def test_haploid_multiallelic_is_still_rejected(self) -> None:
        """E3 must keep failing loudly."""
        with self.assertRaises(BeyondLogicError) as ctx:
            get_ref({"GT": "2"}, "X:37686068_G_A", chrom_copies=1)
        self.assertIn("Multi-allelic", str(ctx.exception))

    def test_diploid_behaviour_is_untouched(self) -> None:
        """A1, A2, E2, E3 - the regression guard for every existing dataset."""
        self.assertEqual(get_ref({"GT": "0/1"}), Zygosity.HET)
        self.assertEqual(get_ref({"GT": "1/1"}), Zygosity.HOM)
        self.assertEqual(get_ref({"GT": "1|0"}), Zygosity.HET)
        self.assertEqual(get_ref({"GT": "./."}), Zygosity.NO_DATA)
        with self.assertRaises(BeyondLogicError):
            get_ref({"GT": "1/2"})


class TestChromCopiesForBg(unittest.TestCase):
    """Both halves are needed: the database half and the sample half."""

    class _FakeVcf:
        def __init__(self, haploid_chroms):
            self.haploid_chroms = haploid_chroms
            self.reference_genome = "GRCh38"

    def _bg(self, bg_type):
        return BloodGroup(type=bg_type, alleles={AlleleState.FILT: []}, sample="s")

    def test_non_par_blood_group_in_a_haploid_sample_is_one_slot(self) -> None:
        self.assertEqual(
            chrom_copies_for_bg(self._bg("XK"), self._FakeVcf(frozenset({"X"})),
                                {"XK": "X"}),
            1,
        )

    def test_par_blood_group_in_a_haploid_sample_stays_two_slots(self) -> None:
        """A6. XG and CD99 are inside PAR1, so a male carries them twice."""
        self.assertEqual(
            chrom_copies_for_bg(self._bg("XG"), self._FakeVcf(frozenset({"X"})),
                                {"XK": "X"}),
            2,
        )

    def test_non_par_blood_group_in_a_diploid_sample_stays_two_slots(self) -> None:
        """The database half alone would make XK haploid for females too."""
        self.assertEqual(
            chrom_copies_for_bg(self._bg("XK"), self._FakeVcf(frozenset()), {"XK": "X"}),
            2,
        )

    def test_autosomal_blood_group_is_always_two_slots(self) -> None:
        self.assertEqual(
            chrom_copies_for_bg(self._bg("FY"), self._FakeVcf(frozenset({"X"})),
                                {"XK": "X"}),
            2,
        )

    def test_empty_pool_still_resolves(self) -> None:
        """B2. The only XK row was dropped as hom ref, and the answer is still one slot.

        This is why the database half exists at all - there is nothing left in the
        variant pool at this point to work it out from.
        """
        bg = self._bg("XK")
        self.assertEqual(bg.alleles[AlleleState.FILT], [])
        self.assertEqual(
            chrom_copies_for_bg(bg, self._FakeVcf(frozenset({"X"})), {"XK": "X"}), 1
        )


class TestFullyHomozygousUsesChromCopies(unittest.TestCase):
    """'On every chromosome present' was hard coded to 2."""

    def setUp(self) -> None:
        self.allele = Allele(
            genotype="XK*N.16",
            phenotype="",
            genotype_alt="",
            phenotype_alt="",
            defining_variants=frozenset({"X:37686068_G_A"}),
            null=True,
            weight_geno=8,
        )

    def test_hemizygous_counts_as_fully_homozygous_when_single_copy(self) -> None:
        homs = get_fully_homozygous_alleles(
            [[self.allele]], {"X:37686068_G_A": 1}, chrom_copies=1
        )
        self.assertEqual(homs, [[self.allele]])

    def test_hemizygous_does_not_when_diploid(self) -> None:
        """C1. A SNV inside a het deletion is one copy of two - not fully homozygous."""
        homs = get_fully_homozygous_alleles(
            [[self.allele]], {"X:37686068_G_A": 1}, chrom_copies=2
        )
        self.assertEqual(homs, [[]])

    def test_default_is_diploid(self) -> None:
        homs = get_fully_homozygous_alleles([[self.allele]], {"X:37686068_G_A": 2})
        self.assertEqual(homs, [[self.allele]])


class TestHaploidGenotypeString(unittest.TestCase):
    """Table decision 1: 'XK*N.03/-'."""

    def _bg(self, chrom_copies):
        allele = Allele(
            genotype="XK*N.16",
            phenotype="",
            genotype_alt="",
            phenotype_alt="",
            defining_variants=frozenset({"X:37686068_G_A"}),
            null=True,
            weight_geno=8,
        )
        from rbceq2.core_logic.alleles import Pair

        return BloodGroup(
            type="XK",
            alleles={AlleleState.CO: None, AlleleState.NORMAL: [Pair(allele, allele)]},
            sample="s",
            chrom_copies=chrom_copies,
        )

    def test_single_copy_renders_the_sentinel(self) -> None:
        bg = get_genotypes({"XK": self._bg(1)})["XK"]
        self.assertEqual(bg.genotypes, [f"XK*N.16/{HAPLOID_SECOND_SLOT}"])

    def test_allele_comes_first(self) -> None:
        """'-' sorts before 'X', so the pair must not be re-sorted as strings."""
        bg = get_genotypes({"XK": self._bg(1)})["XK"]
        self.assertTrue(bg.genotypes[0].startswith("XK*N.16/"))

    def test_two_copies_still_duplicates(self) -> None:
        """B3. A female homozygote, or a ploidy-unaware caller, keeps both slots."""
        bg = get_genotypes({"XK": self._bg(2)})["XK"]
        self.assertEqual(bg.genotypes, ["XK*N.16/XK*N.16"])

    def test_sentinel_is_a_single_character_and_not_an_allele_name(self) -> None:
        self.assertEqual(HAPLOID_SECOND_SLOT, "-")


class TestTokenCopiesInvariant(unittest.TestCase):
    """ploidy_state_table.md section 4, weak form, per token."""

    def _bg(self, pool, chrom_copies):
        bg = BloodGroup(type="XK", alleles={AlleleState.FILT: []}, sample="s",
                        chrom_copies=chrom_copies)
        bg.variant_pool = pool
        return bg

    def test_mixed_ploidy_coding_in_one_sample_is_rejected(self) -> None:
        """'1' at one non-PAR X locus and '1/1' at another - the file says both."""
        bg = self._bg(
            {"X:37686068_G_A": Zygosity.HEM, "X:37694437_C_T": Zygosity.HOM}, 1
        )
        with self.assertRaises(BeyondLogicError) as ctx:
            check_token_copies_fit_chrom_copies(bg)
        self.assertIn("mixes ploidy codings", str(ctx.exception))
        self.assertIn("X:37694437_C_T", str(ctx.exception))

    def test_consistent_haploid_pool_passes(self) -> None:
        bg = self._bg(
            {"X:37686068_G_A": Zygosity.HEM, "X:37694437_C_T": Zygosity.HEM}, 1
        )
        check_token_copies_fit_chrom_copies(bg)

    def test_diploid_pool_is_never_affected(self) -> None:
        """Vacuous while chrom_copies is 2 - Zygosity tops out at HOM, which is 2."""
        bg = self._bg(
            {"1:159205564_G_A": Zygosity.HOM, "1:159204893_T_C": Zygosity.HET}, 2
        )
        check_token_copies_fit_chrom_copies(bg)

    def test_no_data_is_skipped_not_scored(self) -> None:
        """NO_DATA has no len_dict entry, so variant_pool_numeric omits it."""
        bg = self._bg({"X:37686068_G_A": Zygosity.NO_DATA}, 1)
        check_token_copies_fit_chrom_copies(bg)


class TestPloidyScan(unittest.TestCase):
    """Constitutional haploidy read off the raw line, before the filters discard it.

    A sample's haploid genotypes are usually not at a database locus, so row_is_wanted
    drops every one of them. Measured on one per-sample file: 480 haploid called
    genotypes outside PAR, 338 surviving the interval test, none surviving row_is_wanted
    - and the sample, a male, was reported diploid on XK, ATP11C and GATA1.
    """

    SAMPLES = ["s1", "s2"]

    def _row(self, chrom, pos, *gts):
        return [chrom, str(pos), ".", "G", "A", ".", "PASS", ".", "GT", *gts]

    def _scan(self, rows, genome="GRCh38"):
        scan = PloidyScan(genome)
        for row in rows:
            scan.observe(row[0].removeprefix("chr"), int(row[1]), self.SAMPLES, row)
        return scan

    def test_a_haploid_call_outside_par_is_recorded(self) -> None:
        scan = self._scan([self._row("chrX", 37686068, "1", "0/1")])
        self.assertEqual(scan.for_sample("s1"), frozenset({"X"}))
        self.assertEqual(scan.for_sample("s2"), frozenset())

    def test_inside_par_is_not(self) -> None:
        """A haploid GT there is a caller error, not a ploidy statement."""
        scan = self._scan([self._row("chrX", 2000000, "1", "1")])
        self.assertEqual(scan.for_sample("s1"), frozenset())

    def test_an_autosome_is_not(self) -> None:
        """Gene copy number on an autosome is a different question - see get_ref."""
        scan = self._scan([self._row("chr1", 25272548, "1", "1")])
        self.assertEqual(scan.for_sample("s1"), frozenset())

    def test_a_no_call_is_not(self) -> None:
        """The same rule as _infer_haploid_chroms - see gt_names_an_allele."""
        scan = self._scan([self._row("chrX", 37686068, ".", "./.")])
        self.assertEqual(scan.for_sample("s1"), frozenset())

    def test_a_haploid_reference_call_still_counts(self) -> None:
        """'0' is a ploidy statement even though it carries no ALT."""
        scan = self._scan([self._row("chrX", 37686068, "0", "0/0")])
        self.assertEqual(scan.for_sample("s1"), frozenset({"X"}))

    def test_each_sample_is_judged_separately(self) -> None:
        """A cohort mixes both codings in one file, one column each."""
        scan = self._scan([self._row("chrX", 37686068, "1", "0/1")])
        self.assertEqual(scan.for_sample("s1"), frozenset({"X"}))
        self.assertEqual(scan.for_sample("s2"), frozenset())

    def test_a_short_row_does_not_raise(self) -> None:
        """Fewer sample columns than the header promised."""
        scan = PloidyScan("GRCh38")
        scan.observe("X", 37686068, self.SAMPLES, self._row("chrX", 37686068, "1"))
        self.assertEqual(scan.for_sample("s1"), frozenset({"X"}))

    def test_an_unscanned_sample_is_empty_not_missing(self) -> None:
        self.assertEqual(PloidyScan("GRCh38").for_sample("nobody"), frozenset())


class TestObservedHaploidChromsReachesTheVCF(unittest.TestCase):
    """The scan and the frame are unioned, so either one is enough."""

    def _vcf(self, df, observed=None):
        return VCF([df], lane_variants={}, unique_variants=set(), sample="s",
                   reference_genome="GRCh38", observed_haploid_chroms=observed)

    def test_the_scan_supplies_what_the_frame_lost(self) -> None:
        """The whole point - the frame has a diploid row and nothing else."""
        vcf = self._vcf(one_row_df("chr1", "25272548", "0/1"),
                        observed=frozenset({"X"}))
        self.assertEqual(vcf.haploid_chroms, frozenset({"X"}))

    def test_the_frame_still_works_without_a_scan(self) -> None:
        """Every caller written before the scan keeps its behaviour."""
        vcf = self._vcf(one_row_df("chrX", "37686068", "1"))
        self.assertEqual(vcf.haploid_chroms, frozenset({"X"}))

    def test_neither_source_finding_anything_is_diploid(self) -> None:
        vcf = self._vcf(one_row_df("chr1", "25272548", "0/1"))
        self.assertEqual(vcf.haploid_chroms, frozenset())


if __name__ == "__main__":
    unittest.main()
