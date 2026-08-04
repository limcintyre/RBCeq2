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

from rbceq2.core_logic.alleles import Allele, BloodGroup, Pair
from rbceq2.core_logic.constants import (
    AlleleState,
    HAPLOID_SECOND_SLOT,
    UNNAMED_SECOND_SLOT,
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
    """Table decisions 2 and 4: 'RHD*09.01/RHD*01N' and 'GYPB*03/?'."""

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
        self.assertEqual(out.genotypes, [f"GYPB*03/{UNNAMED_SECOND_SLOT}"])

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
        """'?' must not be sorted into the first slot."""
        bg, reference = self._bg("GYPB", "GYPB*03", locus_copies=1)
        out = get_genotypes(
            {"GYPB": bg}, reference_alleles={"GYPB": reference}, gene_absent_subtypes={}
        )["GYPB"]
        self.assertTrue(out.genotypes[0].startswith("GYPB*03/"))

    def test_the_two_glyphs_are_different(self) -> None:
        """'-' says there is no second chromosome, '?' says its allele is unnamed."""
        self.assertNotEqual(HAPLOID_SECOND_SLOT, UNNAMED_SECOND_SLOT)

    def test_without_the_database_maps_nothing_changes(self) -> None:
        """Every pre-v2.4.5 call of get_genotypes keeps its old behaviour."""
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
        """Every pre-v2.4.5 call of add_refs, and the three-argument test calls."""
        self.assertEqual(add_refs(self.db, {}, [], None)["RHD"].genotypes,
                         ["RHD*01/RHD*01"])


if __name__ == "__main__":
    unittest.main()
