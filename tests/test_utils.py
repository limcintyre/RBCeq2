import unittest

from rbceq2.core_logic.alleles import Allele, BloodGroup
from rbceq2.core_logic.utils import (
    apply_to_dict_values,
    check_available_variants,
    chunk_geno_list_by_rank,
    get_non_refs,
    sub_alleles_relationships,
)


# class TestChunklistByRank(unittest.TestCase):
#     def test_empty_list(self):
#         """Test chunking an empty list."""

#     def test_same_weight(self):
#         """Test all alleles with the same weight."""
#         alleles = [
#             Allele("A1", "M", ".", ".", frozenset(), 1, 1),
#             Allele("A2", "M", ".", ".", frozenset(), 1, 1),
#             Allele("A3", "M", ".", ".", frozenset(), 1, 1),
#         ]
#         expected = [alleles]
#         self.assertEqual(chunk_list_by_rank(alleles), expected)

#     def test_different_weights(self):
#         """Test alleles with different weights."""
#         alleles = [
#             Allele("A1", "M", ".", ".", frozenset(), 1, 1),
#             Allele("B1", "N", ".", ".", frozenset(), 1, 2),
#             Allele("C1", "O", ".", ".", frozenset(), 1, 3),
#         ]
#         expected = [[alleles[0]], [alleles[1]], [alleles[2]]]
#         self.assertEqual(chunk_list_by_rank(alleles), expected)

#     def test_mixed_weights(self):
#         """Test alleles with mixed and repeated weights."""
#         alleles = [
#             Allele("A1", "M", ".", ".", frozenset(), 1, 1),
#             Allele("B1", "N", ".", ".", frozenset(), 1, 2),
#             Allele("A2", "M", ".", ".", frozenset(), 1, 1),
#             Allele("C1", "O", ".", ".", frozenset(), 1, 3),
#             Allele("B2", "N", ".", ".", frozenset(), 1, 2),
#         ]
#         expected = [[alleles[0], alleles[2]], [alleles[1], alleles[4]], [alleles[3]]]
#         self.assertEqual(chunk_list_by_rank(alleles), expected)


class TestCheckAvailableVariants(unittest.TestCase):
    def test_all_variants_available(self):
        """Test all variants available"""
        allele = Allele("A1", "M", ".", ".", frozenset({"A", "B", "C"}), False, 1, 1)
        variant_pool = {"A": 2, "B": 2, "C": 2}
        result = check_available_variants(1, variant_pool, lambda x, y: x >= y, allele)
        self.assertTrue(all(result))

    def test_some_variants_unavailable(self):
        """Test some variants unavailable"""
        allele = Allele("A2", "N", ".", ".", frozenset({"A", "B", "C"}),False, 1, 1)
        variant_pool = {"A": 0, "B": 2, "C": 1}
        result = check_available_variants(1, variant_pool, lambda x, y: x >= y, allele)
        self.assertFalse(all(result))

    def test_no_variants_available(self):
        allele = Allele("A3", "O", ".", ".", frozenset({"A", "B", "C"}),False, 1, 1)
        variant_pool = {"A": 0, "B": 0, "C": 0}
        result = check_available_variants(1, variant_pool, lambda x, y: x >= y, allele)
        self.assertFalse(all(result))

    def test_variants_exceeding_count(self):
        allele = Allele("A4", "P", ".", ".", frozenset({"A", "B"}),False, 1, 1)
        variant_pool = {"A": 3, "B": 5}
        result = check_available_variants(2, variant_pool, lambda x, y: x > y, allele)
        self.assertTrue(all(result))


class TestApplyToDictValues(unittest.TestCase):
    def test_modify_blood_groups(self):
        @apply_to_dict_values
        def update_genotypes(
            blood_groups: list[BloodGroup], new_genotype: str
        ) -> list[BloodGroup]:
            for bg in blood_groups:
                bg.genotypes.append(new_genotype)
            return blood_groups

        blood_group1 = BloodGroup("A", alleles={}, sample="Sample1")
        blood_group2 = BloodGroup("B", alleles={}, sample="Sample2")
        input_dict = {"group1": [blood_group1], "group2": [blood_group2]}

        expected_group1_genotypes = ["newType"]
        expected_group2_genotypes = ["newType"]

        updated_dict = update_genotypes(input_dict, "newType")
        self.assertEqual(updated_dict["group1"][0].genotypes, expected_group1_genotypes)
        self.assertEqual(updated_dict["group2"][0].genotypes, expected_group2_genotypes)

    def test_empty_dictionary(self):
        @apply_to_dict_values
        def update_phenotypes(
            blood_groups: list[BloodGroup], new_phenotype: str
        ) -> list[BloodGroup]:
            for bg in blood_groups:
                bg.phenotypes.append(new_phenotype)
            return blood_groups

        self.assertEqual(update_phenotypes({}), {})


class TestGetNonRefs(unittest.TestCase):
    def setUp(self):
        self.alleles = [
            Allele(
                "genotype1",
                "phenotype1",
                ".",
                ".",
                frozenset(["variant1"]),
                False,
                1,
                False,
                "subtype1",
                None,
            ),
            Allele(
                "genotype2",
                "phenotype2",
                ".",
                ".",
                frozenset(["variant2"]),
                False,
                1,
                True,
                "subtype2",
                None,
            ),
            Allele(
                "genotype3",
                "phenotype3",
                ".",
                ".",
                frozenset(["variant3"]),
                False,
                1,
                False,
                "subtype3",
                None,
            ),
        ]

    def test_no_reference_alleles(self):
        options = [self.alleles[0], self.alleles[2]]
        result = get_non_refs(options)
        self.assertEqual(
            result, options, "Should return all alleles when no allele is a reference"
        )

    def test_some_reference_alleles(self):
        result = get_non_refs(self.alleles)
        expected = [self.alleles[0], self.alleles[2]]
        self.assertEqual(result, expected, "Should filter out reference alleles")

    def test_all_reference_alleles(self):
        options = [
            Allele(
                genotype="genotype4",
                phenotype="phenotype4",
                genotype_alt=".",
                phenotype_alt=".",
                defining_variants=frozenset(["variant4"]),
                null=False,
                weight_geno=1,
                reference=True,
                sub_type="subtype4",
                phases=None,
            )
        ]
        result = get_non_refs(options)
        self.assertEqual(
            result, [], "Should return an empty list when all alleles are references"
        )

    def test_empty_list(self):
        result = get_non_refs([])
        self.assertEqual(result, [], "Should handle empty input list gracefully")


class TestChunkGenoListByRank(unittest.TestCase):
    def setUp(self):
        self.alleles = [
            Allele(
                genotype="genotype1",
                phenotype="phenotype1",
                genotype_alt=".",
                phenotype_alt=".",
                defining_variants=frozenset(["variant1"]),
                null=False,
                weight_geno=1,
                reference=False,
                sub_type="A",
            ),
            Allele(
                genotype="genotype2",
                phenotype="phenotype2",
                genotype_alt=".",
                phenotype_alt=".",
                defining_variants=frozenset(["variant2"]),
                null=False,
                weight_geno=2,
                reference=False,
                sub_type="A",
            ),
            Allele(
                genotype="genotype3",
                phenotype="phenotype3",
                genotype_alt=".",
                phenotype_alt=".",
                defining_variants=frozenset(["variant3"]),
                null=False,
                weight_geno=1,
                reference=False,
                sub_type="B",
            ),
            Allele(
                genotype="genotype4",
                phenotype="phenotype4",
                genotype_alt=".",
                phenotype_alt=".",
                defining_variants=frozenset(["variant4"]),
                null=False,
                weight_geno=1,
                reference=False,
                sub_type="B",
            ),
        ]

    def test_empty_input(self):
        self.assertEqual(chunk_geno_list_by_rank([]), [])

    def test_same_weight_and_subtype(self):
        alleles = self.alleles[:2]  # Same subtype "A" with different weights
        result = chunk_geno_list_by_rank(alleles)
        self.assertEqual(len(result), 2)
        self.assertTrue(alleles[0] == result[0][0])
        self.assertTrue(alleles[1] == result[1][0])

    def test_different_weights_same_subtype(self):
        result = chunk_geno_list_by_rank(self.alleles[:2])
        self.assertEqual(len(result), 2)
        self.assertIn(self.alleles[0], result[0])
        self.assertIn(self.alleles[1], result[1])

    def test_different_subtypes_and_weights(self):
        result = chunk_geno_list_by_rank(self.alleles)
        self.assertEqual(len(result), 2)
        self.assertEqual(len(result[0]), 3)
        self.assertEqual(len(result[1]), 1)
        self.assertTrue(
            all(a in result[0] for a in [self.alleles[0]] + self.alleles[2:])
        )
        self.assertTrue(result[1][0] == self.alleles[1])

    def test_mixed_subtypes_identical_weights(self):
        alleles = [
            Allele(
                "genotype5",
                "phenotype5",
                ".",
                ".",
                frozenset(["variant5"]),
                False,
                1,
                False,
                "C",
                None,
            ),
            Allele(
                "genotype6",
                "phenotype6",
                ".",
                ".",
                frozenset(["variant6"]),
                False,
                1,
                False,
                "C",
                None,
            ),
            *self.alleles[:2],
        ]
        result = chunk_geno_list_by_rank(alleles)
        self.assertEqual(len(result), 2)
        self.assertIn(self.alleles[0], result[0])
        self.assertIn(self.alleles[1], result[1])
        self.assertIn(alleles[0], result[0])
        self.assertIn(alleles[1], result[0])


class TestSubAllelesRelationships(unittest.TestCase):
    def setUp(self):
        self.alleles = [
            Allele(
                "genotype1",
                "phenotype1",
                ".",
                ".",
                frozenset(["variant1"]),
                False,
                1,
                False,
                "subtype1",
                None,
            ),
            Allele(
                "genotype2",
                "phenotype2",
                ".",
                ".",
                frozenset(["variant1", "variant2"]),
                False,
                2,
                False,
                "subtype2",
                None,
            ),
        ]

    def test_multiple_alleles(self):
        all_alleles = {"key1": self.alleles}
        relationships, key = sub_alleles_relationships(all_alleles, "key1")
        self.assertTrue(
            relationships[f"{self.alleles[0].genotype}_isin_{self.alleles[1].genotype}"]
        )
        self.assertFalse(
            relationships[f"{self.alleles[1].genotype}_isin_{self.alleles[0].genotype}"]
        )
        self.assertEqual(key, "key1")

    def test_single_allele(self):
        all_alleles = {"key2": [self.alleles[0]]}
        relationships, key = sub_alleles_relationships(all_alleles, "key2")
        self.assertFalse(
            relationships[f"{self.alleles[0].genotype}_isin_{self.alleles[0].genotype}"]
        )
        self.assertEqual(key, "key2")

    def test_no_alleles(self):
        all_alleles = {"key3": []}
        relationships, key = sub_alleles_relationships(all_alleles, "key3")
        self.assertEqual(relationships, {})
        self.assertEqual(key, "key3")

    def test_nonexistent_key(self):
        all_alleles = {"key4": self.alleles}
        with self.assertRaises(KeyError, msg="Expected KeyError for non-existent key"):
            sub_alleles_relationships(all_alleles, "nonexistent_key")


if __name__ == "__main__":
    unittest.main()
