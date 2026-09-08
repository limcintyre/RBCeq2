"""Generic geometry controls for late deletion-based allele placement.

TEST alleles below are synthetic token layouts, not biological definitions.
The companion RHAG pipeline regressions use the curated database.
"""
import unittest

from rbceq2.core_logic.alleles import Allele, BloodGroup, Pair
from rbceq2.core_logic.constants import AlleleState
from rbceq2.core_logic.utils import Zygosity
from rbceq2.filters.geno import (
    cant_pair_with_ref_cuz_a_deletion_names_the_missing_copy,
    cant_split_HEM_SNPs_across_alleles,
)


class TestDeletionPlacement(unittest.TestCase):
    SPLIT = "cant_split_HEM_SNPs_across_alleles"
    REFERENCE = "cant_pair_with_ref_cuz_a_deletion_names_the_missing_copy"

    @staticmethod
    def allele(name, variants=(), reference=False):
        return Allele(
            genotype=f"TEST*{name}", phenotype=".", genotype_alt=".",
            phenotype_alt=".", defining_variants=frozenset(variants),
            null=False, weight_geno=1000, reference=reference, sub_type="TEST*01",
        )

    def fixture(self, a_token="1:150_A_G", b_token="1:160_C_T",
                deletion="1:100_del_100", deletion_state=Zygosity.HET,
                b_state=Zygosity.HEM, counterpart=True, chrom_copies=2, co=None,
                deletion_reference=False):
        ref = self.allele("R", reference=True)
        deleted = self.allele("DEL", [deletion], reference=deletion_reference)
        a, b = self.allele("A", [a_token]), self.allele("B", [b_token])
        named = Pair(deleted, b)
        split = Pair(a, b)
        fallback = Pair(ref, b)
        pairs = ([named] if counterpart else []) + [split, fallback]
        bg = BloodGroup(
            type="TEST", sample="geometry", chrom_copies=chrom_copies,
            locus_copies=None,
            alleles={AlleleState.NORMAL: pairs, AlleleState.CO: co},
            variant_pool={deletion: deletion_state, a_token: Zygosity.HEM, b_token: b_state},
        )
        return bg, named, split, fallback

    @staticmethod
    def apply(bg):
        cant_split_HEM_SNPs_across_alleles({"TEST": bg})
        cant_pair_with_ref_cuz_a_deletion_names_the_missing_copy({"TEST": bg})

    def test_same_deletion_places_both_snps_on_the_surviving_copy(self):
        for state in (Zygosity.HET, Zygosity.HEM):
            with self.subTest(deletion_state=state):
                bg, named, split, fallback = self.fixture(deletion_state=state)
                pool = dict(bg.variant_pool)
                self.apply(bg)
                self.assertEqual(bg.alleles[AlleleState.NORMAL], [named])
                self.assertEqual(bg.filtered_out[self.SPLIT], [split])
                self.assertEqual(bg.filtered_out[self.REFERENCE], [fallback])
                self.assertIsNone(bg.locus_copies)
                self.assertEqual(bg.chrom_copies, 2)
                self.assertEqual(bg.variant_pool, pool)

    def test_without_a_surviving_counterpart_the_existing_paths_remain(self):
        bg, _, split, fallback = self.fixture(counterpart=False)
        self.apply(bg)
        self.assertEqual(bg.alleles[AlleleState.NORMAL], [split, fallback])

    def test_uncalled_or_homozygous_deletion_does_not_supply_placement(self):
        for state in (Zygosity.NO_DATA, Zygosity.NO_COPIES, Zygosity.HOM):
            with self.subTest(state=state):
                bg, named, split, fallback = self.fixture(deletion_state=state)
                self.apply(bg)
                self.assertEqual(bg.alleles[AlleleState.NORMAL], [named, split, fallback])

    def test_heterozygous_snp_does_not_supply_hemizygous_placement(self):
        bg, named, split, fallback = self.fixture(b_state=Zygosity.HET)
        self.apply(bg)
        self.assertEqual(bg.alleles[AlleleState.NORMAL], [named, split, fallback])

    def test_different_chromosome_or_boundary_snp_is_not_split_by_this_rule(self):
        for token in ("2:150_A_G", "1:100_A_G", "1:200_A_G", "1:250_A_G"):
            with self.subTest(token=token):
                bg, named, split, _ = self.fixture(a_token=token)
                cant_split_HEM_SNPs_across_alleles({"TEST": bg})
                self.assertIn(split, bg.alleles[AlleleState.NORMAL])
                self.assertIn(named, bg.alleles[AlleleState.NORMAL])

    def test_allele_outside_the_deletion_keeps_its_reference_alternative(self):
        bg, named, split, fallback = self.fixture(b_token="1:250_C_T")
        self.apply(bg)
        self.assertEqual(bg.alleles[AlleleState.NORMAL], [named, split, fallback])

    def test_synthesised_reference_tokens_are_not_explicit_snp_evidence(self):
        bg, named, split, fallback = self.fixture(b_token="1:160_ref")
        self.apply(bg)
        self.assertEqual(bg.alleles[AlleleState.NORMAL], [named, split, fallback])

    def test_other_structural_types_and_reference_definitions_are_not_replacements(self):
        for options in ({"deletion": "1:100_DUP_100"}, {"deletion_reference": True}):
            with self.subTest(options=options):
                bg, named, split, fallback = self.fixture(**options)
                self.apply(bg)
                self.assertEqual(bg.alleles[AlleleState.NORMAL], [named, split, fallback])

    def test_single_chromosome_and_coexisting_paths_are_untouched(self):
        for options in ({"chrom_copies": 1}, {"co": []}):
            with self.subTest(options=options):
                bg, named, split, fallback = self.fixture(**options)
                self.apply(bg)
                self.assertEqual(bg.alleles[AlleleState.NORMAL], [named, split, fallback])

    def test_missing_deletion_token_cannot_supply_evidence(self):
        bg, named, split, fallback = self.fixture()
        del bg.variant_pool["1:100_del_100"]
        self.apply(bg)
        self.assertEqual(bg.alleles[AlleleState.NORMAL], [named, split, fallback])


if __name__ == "__main__":
    unittest.main()
