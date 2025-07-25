import unittest

from rbceq2.core_logic.alleles import Allele, BloodGroup, Pair


class Zygosity:
    HOM = "Homozygous"
    HET = "Heterozygous"


class TestAllele(unittest.TestCase):
    def setUp(self) -> None:
        self.allele1 = Allele(
            genotype="KN*01.1",
            phenotype="positive",
            genotype_alt="KN*02",
            phenotype_alt="negative",
            defining_variants=frozenset({"variant1", "variant2"}),
            null=False,
            weight_geno=1,
            sub_type="KN*01",
        )
        self.allele2 = Allele(
            genotype="KN*01.2",
            phenotype="negative",
            genotype_alt="KN*03",
            phenotype_alt="positive",
            defining_variants=frozenset({"variant1"}),
            null=False,
            weight_geno=1,
            sub_type="KN*01",
        )
        self.allele3 = Allele(
            genotype="KN*01.3",
            phenotype="positive",
            genotype_alt="KN*04",
            phenotype_alt="negative",
            defining_variants=frozenset({"variant1", "variant2", "variant3"}),
            null=False,
            weight_geno=2,
            sub_type="KN*01",
        )

    def test_post_init(self) -> None:
        self.assertEqual(self.allele1.number_of_defining_variants, 2)
        self.assertEqual(self.allele2.number_of_defining_variants, 1)
        self.assertEqual(self.allele3.number_of_defining_variants, 3)

    def test_contains(self) -> None:
        self.assertIn(self.allele2, self.allele1)
        self.assertNotIn(self.allele1, self.allele2)
        self.assertIn(self.allele1, self.allele3)

    def test_comparison(self) -> None:
        self.assertTrue(self.allele1 > self.allele3)
        self.assertFalse(self.allele1 < self.allele3)
        self.assertTrue(self.allele2 < self.allele1)
        self.assertFalse(self.allele1 < self.allele2)

        unsorted_alleles = [self.allele3, self.allele1, self.allele2]
        sorted_alleles = sorted(unsorted_alleles)
        correct_sorting = [self.allele3, self.allele2, self.allele1]
        self.assertEqual(sorted_alleles, correct_sorting)

    def test_blood_group_property(self) -> None:
        self.assertEqual(self.allele1.blood_group, "KN")
        self.assertEqual(self.allele2.blood_group, "KN")
        self.assertEqual(self.allele3.blood_group, "KN")


class TestBloodGroup(unittest.TestCase):
    def setUp(self) -> None:
        self.allele1 = Allele(
            genotype="KN*01",
            phenotype="positive",
            genotype_alt="KN*02",
            phenotype_alt="negative",
            defining_variants=frozenset({"variant1", "variant2"}),
            null=False,
            weight_geno=2,
        )
        self.allele2 = Allele(
            genotype="KN*02",
            phenotype="negative",
            genotype_alt="KN*03",
            phenotype_alt="positive",
            defining_variants=frozenset({"variant1"}),
            null=False,
            weight_geno=1,
        )

        self.blood_group = BloodGroup(
            type="A",
            alleles={"raw": [self.allele1, self.allele2]},
            sample="sample1",
            variant_pool={"variant1": Zygosity.HOM, "variant2": Zygosity.HET},
            genotypes=["KN*01", "KN*02"],
            phenotypes=["positive", "negative"],
        )

    def test_variant_pool_numeric(self) -> None:
        expected = {"variant1": 2, "variant2": 1}
        self.assertEqual(self.blood_group.variant_pool_numeric, expected)

    def test_number_of_putative_alleles(self) -> None:
        self.assertEqual(self.blood_group.number_of_putative_alleles, 2)

    def test_remove_pairs(self) -> None:
        p1 = Pair(self.allele1, self.allele2)
        self.blood_group.alleles["pairs"] = [p1]
        self.blood_group.remove_pairs([p1], "test_filter", "pairs")
        self.assertEqual(len(self.blood_group.alleles["pairs"]), 0)
        self.assertEqual(len(self.blood_group.filtered_out["test_filter"]), 1)

    def test_remove_alleles(self) -> None:
        self.blood_group.remove_alleles([self.allele1], "test_filter")
        self.assertEqual(len(self.blood_group.alleles["raw"]), 1)
        self.assertEqual(len(self.blood_group.filtered_out["test_filter"]), 1)


class TestPair(unittest.TestCase):
    def setUp(self) -> None:
        self.allele1 = Allele(
            genotype="KN*01",
            phenotype="positive",
            genotype_alt="KN*02",
            phenotype_alt="negative",
            defining_variants=frozenset({"variant1", "variant2"}),
            null=False,
            weight_geno=2,
            reference=True,
        )
        self.allele2 = Allele(
            genotype="KN*02",
            phenotype="negative",
            genotype_alt="KN*03",
            phenotype_alt="positive",
            defining_variants=frozenset({"variant1"}),
            null=False,
            weight_geno=1,
            reference=False,
        )
        self.pair = Pair(allele1=self.allele1, allele2=self.allele2)

    def test_post_init(self) -> None:
        self.assertEqual(self.pair.alleles, frozenset([self.allele1, self.allele2]))

    def test_eq(self) -> None:
        pair2 = Pair(allele1=self.allele1, allele2=self.allele2)
        self.assertEqual(self.pair, pair2)

    def test_contains(self) -> None:
        self.assertIn(self.allele1, self.pair)
        self.assertIn(self.allele2, self.pair)

    def test_iter(self) -> None:
        alleles = list(iter(self.pair))
        self.assertEqual(alleles, [self.allele1, self.allele2])

    def test_repr(self) -> None:
        expected_repr = "Pair(Genotype: KN*01/KN*02 Phenotype: numeric positive/negative alphanumeric: negative/positive) #"
        self.assertEqual(repr(self.pair), expected_repr)

    def test_str(self) -> None:
        expected_str = "Pair(Genotype: KN*01/KN*02 Phenotype: numeric positive/negative alphanumeric: negative/positive) #"
        self.assertEqual(str(self.pair), expected_str)

    def test_ordered(self) -> None:
        self.assertEqual(self.pair._ordered(), [self.allele1, self.allele2])

    def test_genotypes(self) -> None:
        self.assertEqual(self.pair.genotypes, ["KN*01", "KN*02"])

    def test_phenotypes(self) -> None:
        self.assertEqual(self.pair.phenotypes, ["positive", "negative"])

    def test_contains_reference(self) -> None:
        self.assertTrue(self.pair.contains_reference)

    def test_all_reference(self) -> None:
        pair_all_ref = Pair(
            allele1=Allele(
                genotype="KN*01",
                phenotype="positive",
                genotype_alt="KN*02",
                phenotype_alt="negative",
                defining_variants=frozenset({"variant1", "variant2"}),
                null=False,
                weight_geno=2,
                reference=True,
            ),
            allele2=Allele(
                genotype="KN*02",
                phenotype="negative",
                genotype_alt="KN*03",
                phenotype_alt="positive",
                defining_variants=frozenset({"variant1"}),
                null=False,
                weight_geno=1,
                reference=True,
            ),
        )
        self.assertTrue(pair_all_ref.all_reference)

    def test_comparable(self) -> None:
        expected_comparable = [frozenset(["KN*01"]), frozenset(["KN*02"])]
        self.assertEqual(self.pair.comparable, expected_comparable)


if __name__ == "__main__":
    unittest.main()
