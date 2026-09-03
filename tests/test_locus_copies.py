#!/usr/bin/env python3
"""Unit tests for the locus_copies (gene copy number) model - issue #40.

Covers group D of ploidy_state_table.md, the rows where a caller encodes how many copies
of a *gene* are present as genotype ploidy. The state table rows are named in each test so
a failure points back at the spec rather than at the code.

The distinction under test throughout is the one that makes D2 differ from B1 on an
identical genotype:

    chrom_copies  copies of the region the sample was born with -> allele slots
    locus_copies  copies of the gene still present

A gene at one copy on one chromosome is ordinary hemizygosity and reports one slot. A gene
at one copy on *two* chromosomes reports two, the second holding no gene at all. Reading
the second as the first loses a chromosome the sample has.
"""

import unittest

import pandas as pd
from loguru import logger

from rbceq2.core_logic.alleles import Allele, BloodGroup, Pair
from rbceq2.core_logic.constants import (
    AlleleState,
    HAPLOID_SECOND_SLOT,
    NOVEL_DELETION_SLOT,
)
from rbceq2.core_logic.data_procesing import (
    add_refs,
    get_genotypes,
    get_ref,
    locus_copies_for_bg,
)
from rbceq2.core_logic.utils import BeyondLogicError, Zygosity
from rbceq2.db.db import Db, prepare_db, subtype_of
from rbceq2.IO.vcf import VCF


COMMON = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]


def rows_df(rows: list[tuple[str, str, str]]) -> pd.DataFrame:
    """A VCF frame with an array-style FORMAT of just GT. rows are (chrom, pos, GT)."""
    return pd.DataFrame(
        {
            "CHROM": [chrom for chrom, _, _ in rows],
            "POS": [pos for _, pos, _ in rows],
            "ID": ["."] * len(rows),
            "REF": ["G"] * len(rows),
            "ALT": ["A"] * len(rows),
            "QUAL": ["."] * len(rows),
            "FILTER": ["PASS"] * len(rows),
            "INFO": ["."] * len(rows),
            "FORMAT": ["GT"] * len(rows),
            "SAMPLE": [GT for _, _, GT in rows],
        }
    )


def a_vcf(rows: list[tuple[str, str, str]]) -> VCF:
    return VCF(
        [rows_df(rows)],
        lane_variants={},
        unique_variants=set(),
        sample="s",
        reference_genome="GRCh38",
    )


def a_bg(bg_type: str, sample: str = "s") -> BloodGroup:
    return BloodGroup(type=bg_type, alleles={AlleleState.FILT: []}, sample=sample)


class TestSubtypeOf(unittest.TestCase):
    """A copy number carries no breakpoints, so the answer is given at subtype level."""

    def test_strips_the_allele_number(self) -> None:
        self.assertEqual(subtype_of("RHD*01N.01"), "RHD*01N")

    def test_strips_every_trailing_group(self) -> None:
        """The eighteen XK entries that share one token collapse to one string."""
        self.assertEqual(subtype_of("XK*N.01.001"), "XK*N")

    def test_leaves_a_bare_subtype_alone(self) -> None:
        self.assertEqual(subtype_of("GYPB*01N"), "GYPB*01N")

    def test_survives_a_name_with_no_star(self) -> None:
        self.assertEqual(subtype_of("nonsense"), "nonsense")


class TestGeneAbsentSubtypes(unittest.TestCase):
    """Read off the real database - this is the half that decides what gets named."""

    @classmethod
    def setUpClass(cls) -> None:
        cls.db = Db(ref="GRCh38", df=prepare_db())

    def test_rhd_is_named(self) -> None:
        """The case in hand. One deletion defined null, so one answer."""
        self.assertEqual(self.db.gene_absent_subtypes.get("RHD"), "RHD*01N")

    def test_xk_collapses_its_eighteen_entries(self) -> None:
        self.assertEqual(self.db.gene_absent_subtypes.get("XK"), "XK*N")

    def test_glycophorins_are_not_named(self) -> None:
        """A deletion here can fuse two genes rather than remove one.

        GYPB*01N is a deletion that creates a chimera, so the gene is not absent at all.
        Copy number alone cannot tell that from GYPB*05N, which is a genuine whole gene
        deletion, and guessing would report a hybrid as an absence.
        """
        for bg_type in ("GYP", "GYPA", "GYPB"):
            with self.subTest(bg_type=bg_type):
                self.assertNotIn(bg_type, self.db.gene_absent_subtypes)

    def test_gene_conversion_does_not_disqualify(self) -> None:
        """The distinction that keeps RH working.

        RHD has hybrid alleles too, but they are gene conversions - RHD*03N.01 and
        RHD*03N.02 - and a conversion rewrites sequence without changing how many copies
        are present. So it cannot be what a missing copy is. The database says so
        structurally: each of them carries an insertion token beside the deletion, so
        none is defined by a deletion alone.
        """
        self.assertIn("RHD", self.db.gene_absent_subtypes)

    def test_every_answer_is_a_subtype_not_an_allele(self) -> None:
        for bg_type, subtype in self.db.gene_absent_subtypes.items():
            with self.subTest(bg_type=bg_type):
                self.assertNotIn(".", subtype)


class TestInferLocusPloidy(unittest.TestCase):
    """Raw observation, before remove_home_ref can drop the wildtype rows."""

    def test_haploid_and_diploid_are_recorded_apart(self) -> None:
        vcf = a_vcf([("chr1", "25272548", "1"), ("chr1", "25272595", "0/1")])
        self.assertEqual(vcf.haploid_loci.get("1"), frozenset({25272548}))
        self.assertEqual(vcf.diploid_loci.get("1"), frozenset({25272595}))

    def test_a_dropped_hom_ref_row_still_leaves_its_ploidy(self) -> None:
        """The reason this is a separate pass rather than read off the pool later."""
        vcf = a_vcf([("chr1", "25272548", "0")])
        self.assertEqual(len(vcf.df), 0)
        self.assertEqual(vcf.haploid_loci.get("1"), frozenset({25272548}))


class TestLocusCopiesForBg(unittest.TestCase):
    """D2/D3 versus a caller having a bad day."""

    loci = {"RHD": {"1": frozenset({25272548, 25272595, 25284574})}}

    def test_agreement_across_the_gene_reads_as_one_copy(self) -> None:
        """D2. Every locus the VCF reported for this gene is haploid."""
        vcf = a_vcf(
            [
                ("chr1", "25272548", "1"),
                ("chr1", "25272595", "1"),
                ("chr1", "25284574", "1"),
            ]
        )
        self.assertEqual(locus_copies_for_bg(a_bg("RHD"), vcf, self.loci), 1)

    def test_all_wildtype_still_reads_as_one_copy(self) -> None:
        """D3. Every row gets dropped, and this still has to come out at one."""
        vcf = a_vcf(
            [
                ("chr1", "25272548", "0"),
                ("chr1", "25272595", "0"),
                ("chr1", "25284574", "0"),
            ]
        )
        self.assertEqual(len(vcf.df), 0)
        self.assertEqual(locus_copies_for_bg(a_bg("RHD"), vcf, self.loci), 1)

    def test_a_mixed_coding_is_refused(self) -> None:
        """One haploid genotype among diploid ones claims one copy and two at once."""
        vcf = a_vcf([("chr1", "25272548", "1"), ("chr1", "25272595", "0/1")])
        self.assertIsNone(locus_copies_for_bg(a_bg("RHD"), vcf, self.loci))

    def test_a_diploid_gene_says_nothing(self) -> None:
        """None rather than 2 - nothing here can tell two copies from untyped."""
        vcf = a_vcf([("chr1", "25272548", "0/1")])
        self.assertIsNone(locus_copies_for_bg(a_bg("RHD"), vcf, self.loci))

    def test_an_unlisted_blood_group_says_nothing(self) -> None:
        vcf = a_vcf([("chr1", "25272548", "1")])
        self.assertIsNone(locus_copies_for_bg(a_bg("FY"), vcf, self.loci))

    def test_loci_outside_the_gene_are_not_counted(self) -> None:
        """A haploid genotype elsewhere on the chromosome is not this gene's business."""
        vcf = a_vcf([("chr1", "159205564", "1")])
        self.assertIsNone(locus_copies_for_bg(a_bg("RHD"), vcf, self.loci))


class TestGetRefWithLocusCopies(unittest.TestCase):
    """Same zygosity as the constitutional case, different shape of result."""

    def test_haploid_alt_at_one_gene_copy_is_hemizygous(self) -> None:
        """D2."""
        self.assertEqual(
            get_ref({"GT": "1"}, "1:25272548_G_A", chrom_copies=2, locus_copies=1),
            Zygosity.HEM,
        )

    def test_haploid_no_call_at_one_gene_copy_is_no_data(self) -> None:
        self.assertEqual(
            get_ref({"GT": "."}, "1:25272548_G_A", chrom_copies=2, locus_copies=1),
            Zygosity.NO_DATA,
        )

    def test_haploid_with_neither_count_at_one_still_raises(self) -> None:
        """The mixed coding case. Nothing to prefer between one copy and two."""
        with self.assertRaises(BeyondLogicError):
            get_ref({"GT": "1"}, "1:25272548_G_A", chrom_copies=2, locus_copies=None)

    def test_the_message_names_both_counts(self) -> None:
        with self.assertRaises(BeyondLogicError) as caught:
            get_ref({"GT": "1"}, "1:25272548_G_A", chrom_copies=2)
        self.assertIn("locus_copies", str(caught.exception))

    def test_diploid_genotypes_are_untouched_by_it(self) -> None:
        for GT, expected in (("0/1", Zygosity.HET), ("1/1", Zygosity.HOM)):
            with self.subTest(GT=GT):
                self.assertEqual(
                    get_ref({"GT": GT}, "1:25272548_G_A", 2, locus_copies=1), expected
                )


class TestMissingCopyGenotypeString(unittest.TestCase):
    """Table decisions 2 and 4: 'RHD*09.01/RHD*01N' and 'GYPB*03/Novel_gene_deletion'."""

    def _bg(self, bg_type, genotype, locus_copies, chrom_copies=2):
        allele = Allele(
            genotype=genotype,
            phenotype="",
            genotype_alt="",
            phenotype_alt="",
            defining_variants=frozenset({"1:25272548_G_A"}),
            null=False,
            weight_geno=1000,
        )
        reference = Allele(
            genotype=f"{bg_type}*01",
            phenotype="",
            genotype_alt="",
            phenotype_alt="",
            defining_variants=frozenset({"1:25272548_ref"}),
            null=False,
            reference=True,
            weight_geno=1000,
        )
        bg = BloodGroup(
            type=bg_type,
            alleles={
                AlleleState.CO: None,
                AlleleState.NORMAL: [Pair(allele, reference)],
            },
            sample="s",
            chrom_copies=chrom_copies,
        )
        bg.locus_copies = locus_copies
        return bg, reference

    def test_the_reference_slot_is_replaced_by_the_subtype(self) -> None:
        """D2. Pairing with the reference asserts wildtype where the gene is absent."""
        bg, reference = self._bg("RHD", "RHD*09.01", locus_copies=1)
        out = get_genotypes(
            {"RHD": bg},
            reference_alleles={"RHD": reference},
            gene_absent_subtypes={"RHD": "RHD*01N"},
        )["RHD"]
        self.assertEqual(out.genotypes, ["RHD*09.01/RHD*01N"])

    def test_an_unnamable_partner_gets_the_glyph(self) -> None:
        """E1. Honest where the database cannot say which deletion it was."""
        bg, reference = self._bg("GYPB", "GYPB*03", locus_copies=1)
        out = get_genotypes(
            {"GYPB": bg},
            reference_alleles={"GYPB": reference},
            gene_absent_subtypes={},
        )["GYPB"]
        self.assertEqual(out.genotypes, [f"GYPB*03/{NOVEL_DELETION_SLOT}"])

    def test_two_copies_of_the_gene_are_untouched(self) -> None:
        bg, reference = self._bg("RHD", "RHD*09.01", locus_copies=None)
        out = get_genotypes(
            {"RHD": bg},
            reference_alleles={"RHD": reference},
            gene_absent_subtypes={"RHD": "RHD*01N"},
        )["RHD"]
        self.assertEqual(out.genotypes, ["RHD*01/RHD*09.01"])

    def test_one_chromosome_still_wins_over_one_gene_copy(self) -> None:
        """XK in a male is both. One chromosome is the stronger statement - one slot.

        A gene at one copy on one chromosome is ordinary hemizygosity and has no missing
        second slot to name, so it must keep rendering '-' and not 'XK*N.16/XK*N'.
        """
        allele = Allele(
            genotype="XK*N.16",
            phenotype="",
            genotype_alt="",
            phenotype_alt="",
            defining_variants=frozenset({"X:37686068_G_A"}),
            null=True,
            weight_geno=8,
        )
        bg = BloodGroup(
            type="XK",
            alleles={AlleleState.CO: None, AlleleState.NORMAL: [Pair(allele, allele)]},
            sample="s",
            chrom_copies=1,
        )
        bg.locus_copies = 1
        out = get_genotypes(
            {"XK": bg}, gene_absent_subtypes={"XK": "XK*N"}
        )["XK"]
        self.assertEqual(out.genotypes, [f"XK*N.16/{HAPLOID_SECOND_SLOT}"])

    def test_the_named_allele_comes_first(self) -> None:
        """The absent copy marker must not be sorted into the first slot."""
        bg, reference = self._bg("GYPB", "GYPB*03", locus_copies=1)
        out = get_genotypes(
            {"GYPB": bg}, reference_alleles={"GYPB": reference}, gene_absent_subtypes={}
        )["GYPB"]
        self.assertTrue(out.genotypes[0].startswith("GYPB*03/"))

    def test_the_two_glyphs_are_different(self) -> None:
        """'-' says there is no second chromosome, 'Novel_gene_deletion' says there is one
        carrying no copy of the gene."""
        self.assertNotEqual(HAPLOID_SECOND_SLOT, NOVEL_DELETION_SLOT)

    def test_without_the_database_maps_nothing_changes(self) -> None:
        """Every call written before the database maps were added keeps its old
        behaviour."""
        bg, _ = self._bg("RHD", "RHD*09.01", locus_copies=None)
        out = get_genotypes({"RHD": bg})["RHD"]
        self.assertEqual(out.genotypes, ["RHD*01/RHD*09.01"])


class TestAddRefsMissingCopy(unittest.TestCase):
    """D3 through the path that has no alleles left to reason from.

    The common shape rather than an edge case. A sample missing one copy of a gene is
    usually ordinary at the loci that remain, so every row is a haploid '0', every row is
    dropped, and the blood group arrives here with nothing in it.
    """

    @classmethod
    def setUpClass(cls) -> None:
        cls.db = Db(ref="GRCh38", df=prepare_db())

    def _vcf_all_wildtype(self, bg_type: str) -> VCF:
        positions = sorted(self.db.loci_by_type[bg_type]["1"])
        return VCF(
            [rows_df([("chr1", str(pos), "0") for pos in positions])],
            lane_variants=self.db.lane_variants,
            unique_variants=self.db.unique_variants,
            sample="s",
            reference_genome="GRCh38",
        )

    def test_one_gene_copy_all_wildtype_names_the_absence(self) -> None:
        res = add_refs(self.db, {}, [], self._vcf_all_wildtype("RHD"))
        self.assertEqual(res["RHD"].genotypes, ["RHD*01/RHD*01N"])
        self.assertEqual(res["RHD"].chrom_copies, 2)
        self.assertEqual(res["RHD"].locus_copies, 1)

    def test_an_untouched_blood_group_keeps_both_slots(self) -> None:
        res = add_refs(self.db, {}, [], self._vcf_all_wildtype("RHD"))
        self.assertEqual(res["FY"].genotypes, ["FY*01/FY*01"])

    def test_no_vcf_keeps_the_old_behaviour(self) -> None:
        """Every call written before add_refs took a vcf, and the three-argument test
        calls."""
        self.assertEqual(add_refs(self.db, {}, [], None)["RHD"].genotypes,
                         ["RHD*01/RHD*01"])


class TestLociByTypeExcludesStructuralTokens(unittest.TestCase):
    """A breakpoint is not a locus, read off the real database.

    A structural token's position is where an event starts or where its replacement
    sequence came from. Nobody genotypes it, and on a paralogue pair it names the other
    gene, so it must not be one of the positions that votes on copy number.
    """

    @classmethod
    def setUpClass(cls) -> None:
        cls.db = Db(ref="GRCh38", df=prepare_db())

    def test_the_paralogues_share_no_position(self) -> None:
        """Counting every token position made these two sets share eight."""
        rhd = self.db.loci_by_type["RHD"]["1"]
        rhce = self.db.loci_by_type["RHCE"]["1"]
        self.assertEqual(rhd & rhce, frozenset())

    def test_the_paralogues_do_not_straddle_each_other(self) -> None:
        """Nothing here knows where a gene starts, and it does not need to.

        Dropping structural tokens is enough on its own to separate the two genes into
        adjacent blocks - which is the check that the rule is picking up real biology
        rather than just happening to delete the dissenters.
        """
        rhd = self.db.loci_by_type["RHD"]["1"]
        rhce = self.db.loci_by_type["RHCE"]["1"]
        self.assertLess(max(rhd), min(rhce))

    def test_a_donor_coordinate_does_not_vote_for_the_borrower(self) -> None:
        """RHD*01N.43 is 25272547_DEL_18244 in RHD plus 25402595_INS_18269 in RHCE.

        25402595 is an RHCE locus - RHCE*01N.13 and RHCE*02N.03 both define an ordinary
        SNP there. It was RHD's single diploid dissenter.
        """
        self.assertNotIn(25402595, self.db.loci_by_type["RHD"]["1"])
        self.assertIn(25402595, self.db.loci_by_type["RHCE"]["1"])

    def test_the_mirror_case_holds_too(self) -> None:
        """RHCE*01.44 carries 25301518_INS_1940, and 25301518 is an RHD locus."""
        self.assertNotIn(25301518, self.db.loci_by_type["RHCE"]["1"])
        self.assertIn(25301518, self.db.loci_by_type["RHD"]["1"])

    def test_an_ordinary_snp_still_votes(self) -> None:
        """25317062_T_C, on RHD*01N.43 beside the two structural tokens."""
        self.assertIn(25317062, self.db.loci_by_type["RHD"]["1"])

    def test_a_ref_token_still_votes(self) -> None:
        """'no change from reference' is a locus like any other, and D3 needs it."""
        self.assertIn(25408711, self.db.loci_by_type["RHCE"]["1"])

    def test_a_gene_defined_only_by_a_deletion_has_no_voters(self) -> None:
        """ABCC1's whole database entry is one whole gene deletion.

        Absent rather than present-and-empty, and that is the honest answer: one
        breakpoint cannot show agreement across a gene, and with a single position
        'every locus agrees' is satisfied by whatever one row lands on it.
        """
        for bg_type in ("ABCC1", "ATP11C", "CD99"):
            with self.subTest(bg_type=bg_type):
                self.assertNotIn(bg_type, self.db.loci_by_type)


class TestRhIsNoLongerRefusedOutright(unittest.TestCase):
    """The regression. Both RH genes were refused on every sample of a cohort.

    The log was 'RHD has 204 haploid and 1 diploid genotypes' and 'RHCE has 3 haploid
    and 95 diploid'. The single RHD dissenter was 25402595, and one of the three RHCE
    ones was 25301518 - each borrowed from the other gene, each enough on its own to
    refuse a gene that wants agreement across all of it.

    Built here as the shape that produced it: one gene at one copy, its neighbour at
    two. The cohort itself is not needed to reproduce it - the database is.
    """

    @classmethod
    def setUpClass(cls) -> None:
        cls.db = Db(ref="GRCh38", df=prepare_db())

    def _vcf_rhd_haploid_rhce_diploid(self) -> VCF:
        rows = [("chr1", str(pos), "1") for pos in self.db.loci_by_type["RHD"]["1"]]
        rows += [("chr1", str(pos), "0/1") for pos in self.db.loci_by_type["RHCE"]["1"]]
        return a_vcf(sorted(rows, key=lambda row: int(row[1])))

    def test_the_single_copy_gene_is_read_as_one_copy(self) -> None:
        vcf = self._vcf_rhd_haploid_rhce_diploid()
        self.assertEqual(
            locus_copies_for_bg(a_bg("RHD"), vcf, self.db.loci_by_type), 1
        )

    def test_its_two_copy_neighbour_is_still_not(self) -> None:
        """RHCE's shape is the opposite one, and refusing it is the right answer."""
        vcf = self._vcf_rhd_haploid_rhce_diploid()
        self.assertIsNone(
            locus_copies_for_bg(a_bg("RHCE"), vcf, self.db.loci_by_type)
        )


class TestTheContradictionWarning(unittest.TestCase):
    """1b. Which contradictions are worth saying out loud.

    Unanimity still decides the result. This is only about the warning, which fired on
    every sample of a cohort and so trained people to ignore it.
    """

    loci = {"RHD": {"1": frozenset({25272548, 25272595, 25284574, 25284578, 25290660})}}

    def _run(self, rows: list[tuple[str, str, str]]) -> tuple[int | None, list[str]]:
        messages: list[str] = []
        sink = logger.add(
            lambda m: messages.append(m.record["message"]), level="WARNING"
        )
        try:
            result = locus_copies_for_bg(a_bg("RHD"), a_vcf(rows), self.loci)
        finally:
            logger.remove(sink)
        return result, messages

    def test_a_gene_that_is_mostly_haploid_is_worth_saying(self) -> None:
        """Unanimity is all that stands in the way, so it may be costing a call."""
        result, messages = self._run(
            [
                ("chr1", "25272548", "1"),
                ("chr1", "25272595", "1"),
                ("chr1", "25284574", "1"),
                ("chr1", "25284578", "0/1"),
            ]
        )
        self.assertIsNone(result)
        self.assertEqual(len(messages), 1)

    def test_and_it_names_the_dissenters(self) -> None:
        """They are the whole diagnosis - the point is to go and look at them."""
        _, messages = self._run(
            [
                ("chr1", "25272548", "1"),
                ("chr1", "25272595", "1"),
                ("chr1", "25284574", "1"),
                ("chr1", "25284578", "0/1"),
            ]
        )
        self.assertIn("1:25284578", messages[0])
        self.assertNotIn("1:25272548", messages[0])

    def test_a_gene_that_is_mostly_diploid_is_not(self) -> None:
        """A diploid gene and a caller having a bad day. The ordinary case."""
        result, messages = self._run(
            [
                ("chr1", "25272548", "0/1"),
                ("chr1", "25272595", "0/1"),
                ("chr1", "25284574", "0/1"),
                ("chr1", "25284578", "1"),
            ]
        )
        self.assertIsNone(result)
        self.assertEqual(messages, [])

    def test_a_tie_is_not_evidence_of_a_copy_number(self) -> None:
        """Deliberate, and the judgement call in this rule.

        A tie says the two readings are equally supported, which is a mess rather than a
        copy number. On a sparsely called input a tie is usually one locus against one,
        and warning there would put the noise straight back.
        """
        result, messages = self._run(
            [("chr1", "25272548", "1"), ("chr1", "25272595", "0/1")]
        )
        self.assertIsNone(result)
        self.assertEqual(messages, [])

    def test_the_result_does_not_depend_on_any_of_this(self) -> None:
        """Whether it speaks or not, unanimity decides, exactly as before."""
        loud, _ = self._run(
            [
                ("chr1", "25272548", "1"),
                ("chr1", "25272595", "1"),
                ("chr1", "25284574", "0/1"),
            ]
        )
        quiet, _ = self._run(
            [
                ("chr1", "25272548", "1"),
                ("chr1", "25272595", "0/1"),
                ("chr1", "25284574", "0/1"),
            ]
        )
        self.assertIsNone(loud)
        self.assertIsNone(quiet)

    def test_agreement_still_says_nothing_at_all(self) -> None:
        result, messages = self._run(
            [("chr1", "25272548", "1"), ("chr1", "25272595", "1")]
        )
        self.assertEqual(result, 1)
        self.assertEqual(messages, [])


if __name__ == "__main__":
    unittest.main()
