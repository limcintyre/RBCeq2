import unittest
from collections import defaultdict

from rbceq2.core_logic.alleles import Allele, BloodGroup, Pair, Zygosity
from rbceq2.core_logic.constants import AlleleState
from rbceq2.filters.geno import (
    ABO_cant_pair_with_ref_cuz_261delG_HET,
    cant_pair_with_ref_cuz_SNPs_must_be_on_other_side,
    cant_pair_with_ref_cuz_trumped,
    ensure_co_existing_HET_SNP_used,
    filter_co_existing_in_other_allele,
    filter_co_existing_pairs,
    filter_co_existing_subsets,
    filter_co_existing_with_normal,
    filter_coexisting_pairs_on_antithetical_zygosity,
    filter_HET_pairs_by_weight,
    filter_pairs_by_context,
    filter_pairs_by_phase,
    filter_pairs_on_antithetical_modyfying_SNP,
    filter_pairs_on_antithetical_zygosity,
    flatten_alleles,
    parse_bio_info2,
    split_pair_by_ref,
)


class TestFlattenAlleles(unittest.TestCase):
    def test_flatten_alleles(self):
        allele1 = Allele(
            genotype="A1",
            phenotype="M",
            genotype_alt=".",
            phenotype_alt=".",
            defining_variants=frozenset({"var1"}),
            weight_geno=1,
            reference=False,
            sub_type="subtype1",
        )
        allele2 = Allele(
            genotype="A2",
            phenotype="N",
            genotype_alt=".",
            phenotype_alt=".",
            defining_variants=frozenset({"var2"}),
            weight_geno=2,
            reference=False,
            sub_type="subtype2",
        )
        pair1 = Pair(allele1=allele1, allele2=allele2)
        pair2 = Pair(allele1=allele2, allele2=allele1)
        # Should not add duplicates due to set behavior

        expected = {allele1, allele2}
        result = flatten_alleles([pair1, pair2])
        self.assertEqual(
            result, expected, "Should return a unique set of alleles from pairs"
        )

    def test_empty_list(self):
        expected = set()
        result = flatten_alleles([])
        self.assertEqual(
            result, expected, "Should return an empty set when input is an empty list"
        )

    def test_all_identical_pairs(self):
        allele = Allele(
            genotype="A1",
            phenotype="M",
            genotype_alt=".",
            phenotype_alt=".",
            defining_variants=frozenset({"var1"}),
            weight_geno=1,
            reference=False,
            sub_type="subtype1",
        )
        pair = Pair(allele1=allele, allele2=allele)
        expected = {allele}
        result = flatten_alleles([pair, pair])
        self.assertEqual(
            result,
            expected,
            "Should return a set with a single allele when all pairs are identical",
        )


class TestSplitPairByRef(unittest.TestCase):
    def test_normal_case(self):
        allele1 = Allele(
            genotype="A1",
            phenotype="M",
            genotype_alt=".",
            phenotype_alt=".",
            defining_variants=frozenset({"var1"}),
            weight_geno=1,
            reference=True,
            sub_type="subtype1",
        )
        allele2 = Allele(
            genotype="A2",
            phenotype="N",
            genotype_alt=".",
            phenotype_alt=".",
            defining_variants=frozenset({"var2"}),
            weight_geno=2,
            reference=False,
            sub_type="subtype2",
        )
        pair = Pair(allele1=allele1, allele2=allele2)
        ref, non_ref = split_pair_by_ref(pair)
        self.assertEqual(ref, allele1)
        self.assertEqual(non_ref, allele2)

    def test_both_reference(self):
        allele1 = Allele(
            genotype="A1",
            phenotype="M",
            genotype_alt=".",
            phenotype_alt=".",
            defining_variants=frozenset({"var1"}),
            weight_geno=1,
            reference=True,
            sub_type="subtype1",
        )
        allele2 = Allele(
            genotype="A2",
            phenotype="N",
            genotype_alt=".",
            phenotype_alt=".",
            defining_variants=frozenset({"var2"}),
            weight_geno=2,
            reference=True,
            sub_type="subtype2",
        )
        pair = Pair(allele1=allele1, allele2=allele2)
        with self.assertRaises(ValueError):
            split_pair_by_ref(pair)

    def test_neither_reference(self):
        allele1 = Allele(
            genotype="A1",
            phenotype="M",
            genotype_alt=".",
            phenotype_alt=".",
            defining_variants=frozenset({"var1"}),
            weight_geno=1,
            reference=False,
            sub_type="subtype1",
        )
        allele2 = Allele(
            genotype="A2",
            phenotype="N",
            genotype_alt=".",
            phenotype_alt=".",
            defining_variants=frozenset({"var2"}),
            weight_geno=2,
            reference=False,
            sub_type="subtype2",
        )
        pair = Pair(allele1=allele1, allele2=allele2)
        with self.assertRaises(ValueError):
            split_pair_by_ref(pair)


class TestParseBioInfo2(unittest.TestCase):
    def test_multiple_pairs(self):
        allele1 = Allele(
            genotype="A1+A2",
            phenotype="M",
            genotype_alt=".",
            phenotype_alt=".",
            defining_variants=frozenset({"var1"}),
            weight_geno=1,
            reference=False,
            sub_type="subtype1",
        )
        allele2 = Allele(
            genotype="B1+B2",
            phenotype="N",
            genotype_alt=".",
            phenotype_alt=".",
            defining_variants=frozenset({"var2"}),
            weight_geno=2,
            reference=False,
            sub_type="subtype2",
        )
        pair1 = Pair(allele1=allele1, allele2=allele2)
        pair2 = Pair(allele1=allele2, allele2=allele1)
        result = parse_bio_info2([pair1, pair2])
        expected = [
            [frozenset({"A1", "A2"}), frozenset({"B1", "B2"})],
            [frozenset({"B1", "B2"}), frozenset({"A1", "A2"})],
        ]
        self.assertEqual(result, expected)

    def test_empty_list(self):
        result = parse_bio_info2([])
        expected = []
        self.assertEqual(result, expected)


class TestFilterCoExistingPairs(unittest.TestCase):
    def test_normal_case(self):
        allele1 = Allele(
            genotype="A1",
            phenotype="M",
            genotype_alt="mushed",
            phenotype_alt=".",
            defining_variants=frozenset({"var1"}),
            weight_geno=1,
            reference=False,
            sub_type="A/B",
        )
        allele2 = Allele(
            genotype="A2",
            phenotype="N",
            genotype_alt=".",
            phenotype_alt=".",
            defining_variants=frozenset({"var2"}),
            weight_geno=2,
            reference=False,
            sub_type="C/D",
        )
        pair = Pair(allele1=allele1, allele2=allele2)
        bg = BloodGroup(
            type="ExampleType",
            alleles={AlleleState.CO: [pair]},
            sample="sample1",
            variant_pool={},
        )

        bg_filtered = list(filter_co_existing_pairs({1: bg}).values())[0]
        self.assertFalse(bg_filtered.alleles[AlleleState.CO])

    def test_no_removal_required(self):
        allele1 = Allele(
            genotype="A1",
            phenotype="M",
            genotype_alt=".",
            phenotype_alt=".",
            defining_variants=frozenset({"var1"}),
            weight_geno=1,
            reference=False,
            sub_type="A/B",
        )
        allele2 = Allele(
            genotype="A2",
            phenotype="N",
            genotype_alt=".",
            phenotype_alt=".",
            defining_variants=frozenset({"var2"}),
            weight_geno=2,
            reference=False,
            sub_type="C/D",
        )
        pair = Pair(allele1=allele1, allele2=allele2)
        bg = BloodGroup(
            type="ExampleType",
            alleles={AlleleState.CO: [pair]},
            sample="sample1",
            variant_pool={},
        )

        bg_filtered = list(filter_co_existing_pairs({1: bg}).values())[0]
        self.assertTrue(bg_filtered.alleles[AlleleState.CO])

    def test_empty_co_existing_list(self):
        bg = BloodGroup(
            type="ExampleType",
            alleles={AlleleState.CO: []},
            sample="sample1",
            variant_pool={},
        )

        bg_filtered = list(filter_co_existing_pairs({1: bg}).values())[0]
        self.assertEqual(bg_filtered.alleles[AlleleState.CO], [])


class TestFilterCoExistingInOtherAllele(unittest.TestCase):
    def test_pairs_removed(self):
        allele1 = Allele(
            genotype="FUT2*01.03.01",
            phenotype="Type1",
            genotype_alt="mushed",
            phenotype_alt=".",
            defining_variants=frozenset({"var1"}),
            weight_geno=1,
            reference=False,
            sub_type="subtype1",
        )
        allele2 = Allele(
            genotype="FUT2*01.06",
            phenotype="Type2",
            genotype_alt=".",
            phenotype_alt=".",
            defining_variants=frozenset({"var1"}),
            weight_geno=2,
            reference=False,
            sub_type="subtype2",
        )
        pair = Pair(allele1=allele1, allele2=allele2)
        raw_allele = Allele(
            genotype="FUT2*01.03.03",
            phenotype="Type3",
            genotype_alt=".",
            phenotype_alt=".",
            defining_variants=frozenset({"var1"}),
            weight_geno=3,
            reference=False,
            sub_type="subtype3",
        )
        bg = BloodGroup(
            type="ExampleType",
            alleles={AlleleState.CO: [pair], AlleleState.FILT: [raw_allele]},
            sample="sample1",
            variant_pool={},
        )
        bg_filtered = list(filter_co_existing_in_other_allele({1: bg}).values())[0]
        self.assertFalse(bg_filtered.alleles[AlleleState.CO])

    def test_no_removal_required(self):
        allele1 = Allele(
            genotype="FUT2*01.03.01",
            phenotype="Type1",
            genotype_alt=".",
            phenotype_alt=".",
            defining_variants=frozenset({"var2"}),
            weight_geno=1,
            reference=False,
            sub_type="subtype1",
        )
        allele2 = Allele(
            genotype="FUT2*01.06",
            phenotype="Type2",
            genotype_alt=".",
            phenotype_alt=".",
            defining_variants=frozenset({"var3"}),
            weight_geno=2,
            reference=False,
            sub_type="subtype2",
        )
        pair = Pair(allele1=allele1, allele2=allele2)
        raw_allele = Allele(
            genotype="FUT2*01.03.03",
            phenotype="Type3",
            genotype_alt=".",
            phenotype_alt=".",
            defining_variants=frozenset({"var1"}),
            weight_geno=3,
            reference=False,
            sub_type="subtype3",
        )
        bg = BloodGroup(
            type="ExampleType",
            alleles={AlleleState.CO: [pair], "raw": [raw_allele]},
            sample="sample1",
            variant_pool={},
        )

        bg_filtered = list(filter_co_existing_in_other_allele({1: bg}).values())[0]
        self.assertTrue(bg_filtered.alleles[AlleleState.CO])

    def test_empty_co_existing_list(self):
        bg = BloodGroup(
            type="ExampleType",
            alleles={AlleleState.CO: [], "raw": []},
            sample="sample1",
            variant_pool={},
        )

        bg_filtered = list(filter_co_existing_in_other_allele({1: bg}).values())[0]
        self.assertEqual(bg_filtered.alleles[AlleleState.CO], [])


class TestFilterCoExistingWithNormal(unittest.TestCase):
    def setUp(self):
        self.allele1 = Allele(
            genotype="KN1",
            phenotype="Type1",
            genotype_alt=".",
            phenotype_alt=".",
            defining_variants=frozenset({"var1"}),
            weight_geno=1,
            reference=False,
            sub_type="subtype1",
        )

        self.allele2 = Allele(
            genotype="KN2",
            phenotype="Type2",
            genotype_alt=".",
            phenotype_alt=".",
            defining_variants=frozenset({"var2"}),
            weight_geno=2,
            reference=False,
            sub_type="subtype2",
        )
        # Normal pair
        self.pair1 = Pair(allele1=self.allele1, allele2=self.allele2)
        self.allele3 = Allele(
            genotype="KN1+KN2",
            phenotype="Type2",
            genotype_alt="mushed",
            phenotype_alt=".",
            defining_variants=frozenset({"var1, var2"}),
            weight_geno=2,
            reference=False,
            sub_type="subtype2",
        )
        # Pair with coexisting alleles
        self.pair2 = Pair(allele1=self.allele1, allele2=self.allele3)

    def test_normal_case(self):
        bg = BloodGroup(
            type="ExampleType",
            alleles={AlleleState.NORMAL: [self.pair1], AlleleState.CO: [self.pair1]},
            sample="sample230",
            variant_pool={},
            filtered_out=defaultdict(list),
        )
        self.assertTrue(self.pair1 in bg.alleles[AlleleState.NORMAL])
        self.assertTrue(self.pair1 in bg.alleles[AlleleState.CO])
        bg.remove_pairs([self.pair1], "normal_filter", AlleleState.NORMAL)
        self.assertTrue(self.pair1 not in bg.alleles[AlleleState.NORMAL])
        filtered_bg = list(filter_co_existing_with_normal({1: bg}).values())[0]
        self.assertTrue(
            self.pair1 in filtered_bg.filtered_out["filter_co_existing_with_normal"]
        )
        self.assertEqual(filtered_bg.alleles, {"co_existing": [], "pairs": []})

    def test_no_removal_needed(self):
        # Modify pair2 to meet exclusion criteria
        bg = BloodGroup(
            type="ExampleType",
            alleles={AlleleState.NORMAL: [self.pair1], AlleleState.CO: [self.pair2]},
            sample="sample230",
            variant_pool={},
            filtered_out=defaultdict(list),
        )
        filtered_bg = list(filter_co_existing_with_normal({1: bg}).values())[0]
        self.assertTrue(self.pair2 in filtered_bg.alleles[AlleleState.CO])
        self.assertEqual(len(filtered_bg.alleles[AlleleState.CO]), 1)
        self.assertEqual(len(filtered_bg.alleles[AlleleState.NORMAL]), 1)

    def test_empty_co_existing_list(self):
        bg = BloodGroup(
            type="ExampleType",
            alleles={AlleleState.NORMAL: [self.pair1], AlleleState.CO: []},
            sample="sample230",
            variant_pool={},
            filtered_out=defaultdict(list),
        )
        filtered_bg = list(filter_co_existing_with_normal({1: bg}).values())[0]
        self.assertEqual(filtered_bg.alleles[AlleleState.CO], [])


class TestFilterCoExistingSubsets(unittest.TestCase):
    def setUp(self):
        self.allele1 = Allele(
            genotype="KN*01.06",
            phenotype="Type2",
            genotype_alt=".",
            phenotype_alt=".",
            defining_variants=frozenset(
                {"207782856_A_G, 207782916_A_T, 207782889_A_G, 207782931_A_G"}
            ),
            weight_geno=1,
            reference=False,
            sub_type="subtype1",
        )
        self.allele2 = Allele(
            genotype="KN*01.07",
            phenotype="Type2",
            genotype_alt=".",
            phenotype_alt=".",
            defining_variants=frozenset({"207782916_A_T, 207782889_A_G"}),
            weight_geno=1,
            reference=False,
            sub_type="subtype1",
        )
        self.allele3 = Allele(
            genotype="KN*01.10",
            phenotype="Type2",
            genotype_alt=".",
            phenotype_alt=".",
            defining_variants=frozenset({" 207782916_A_T, 207782931_A_G"}),
            weight_geno=1,
            reference=False,
            sub_type="subtype1",
        )
        self.allele4 = Allele(
            genotype="KN*01.07+KN*01.10",
            phenotype="Type2",
            genotype_alt="mushed",
            phenotype_alt=".",
            defining_variants=frozenset(
                {"207782916_A_T, 207782889_A_G, 207782931_A_G"}
            ),
            weight_geno=1,
            reference=False,
            sub_type="subtype1",
        )
        self.pair1 = Pair(allele1=self.allele1, allele2=self.allele4)  # Not a subset
        self.pair2 = Pair(allele1=self.allele4, allele2=self.allele4)  # Subset

    def test_pairs_removed(self):
        bg = BloodGroup(
            type="ExampleType",
            alleles={AlleleState.CO: [self.pair1, self.pair2], AlleleState.RAW: []},
            sample="sample230",
            variant_pool={
                "1:207782856_A_G": "Heterozygous",
                "1:207782889_A_G": "Homozygous",
                "1:207782916_A_T": "Homozygous",
                "1:207782931_A_G": "Homozygous",
            },
            filtered_out=defaultdict(list),
        )
        self.assertTrue(self.pair2 in bg.alleles[AlleleState.CO])
        filtered_bg = list(filter_co_existing_subsets({1: bg}).values())[0]

        self.assertTrue(self.pair2 not in filtered_bg.alleles[AlleleState.CO])
        self.assertTrue(
            self.pair2 in filtered_bg.filtered_out["filter_co_existing_subsets"]
        )

    def test_no_pairs_removed(self):
        bg = BloodGroup(
            type="ExampleType",
            alleles={AlleleState.CO: [self.pair1, self.pair2], AlleleState.RAW: []},
            sample="sample230",
            variant_pool={},
            filtered_out=defaultdict(list),
        )
        # Reverse the condition so no subsets exist
        filtered_bg = list(filter_co_existing_subsets({1: bg}).values())[0]
        self.assertTrue(self.pair1 in filtered_bg.alleles[AlleleState.CO])
        self.assertTrue(
            self.pair1 not in filtered_bg.filtered_out["filter_co_existing_subsets"]
        )  ###TODO wtf - is 'in' not working for Pair objects?

    def test_empty_co_existing_list(self):
        bg = BloodGroup(
            type="ExampleType",
            alleles={AlleleState.CO: []},
            sample="sample230",
            variant_pool={},
            filtered_out=defaultdict(list),
        )
        filtered_bg = list(filter_co_existing_subsets({1: bg}).values())[0]
        self.assertEqual(filtered_bg.alleles[AlleleState.CO], [])

    def test_permutation_pairs(self):
        # Create Allele objects
        alleleA = Allele(
            genotype="AlleleA",
            phenotype="TypeA",
            genotype_alt=".",
            phenotype_alt=".",
            defining_variants=frozenset({"varA"}),
            weight_geno=1,
            reference=False,
            sub_type="subtype1",
        )
        alleleB = Allele(
            genotype="AlleleB",
            phenotype="TypeB",
            genotype_alt=".",
            phenotype_alt=".",
            defining_variants=frozenset({"varB"}),
            weight_geno=1,
            reference=False,
            sub_type="subtype1",
        )

        # Create pairs that are permutations of each other
        pair1 = Pair(allele1=alleleA, allele2=alleleB)
        pair2 = Pair(allele1=alleleB, allele2=alleleA)

        # Prepare the BloodGroup object
        bg = BloodGroup(
            type="ExampleType",
            alleles={
                AlleleState.CO: [pair1, pair2],
                AlleleState.RAW: [alleleA, alleleB],
            },
            sample="sample231",
            variant_pool={"varA": Zygosity.HOM, "varB": Zygosity.HOM},
            filtered_out=defaultdict(list),
        )

        # Now run the function
        filtered_bg = list(filter_co_existing_subsets({1: bg}).values())[0]

        # Both pairs should remain because they are permutations
        self.assertIn(pair1, filtered_bg.alleles[AlleleState.CO])
        self.assertIn(pair2, filtered_bg.alleles[AlleleState.CO])

    def test_alleles_with_one_het_variant_not_same_subtype(self):
        # Create Allele objects
        allele1 = Allele(
            genotype="Allele1",
            phenotype="Type1",
            genotype_alt=".",
            phenotype_alt=".",
            defining_variants=frozenset({"var1"}),
            weight_geno=1,
            reference=False,
            sub_type="subtype1",
        )
        allele2 = Allele(
            genotype="Allele2",
            phenotype="Type2",
            genotype_alt=".",
            phenotype_alt=".",
            defining_variants=frozenset({"var2", "var3"}),
            weight_geno=1,
            reference=False,
            sub_type="subtype2",  # Different subtype
        )

        pair = Pair(allele1=allele1, allele2=allele2)

        # Prepare the BloodGroup object
        bg = BloodGroup(
            type="ExampleType",
            alleles={AlleleState.CO: [pair], AlleleState.RAW: [allele1, allele2]},
            sample="sample232",
            variant_pool={
                "var1": Zygosity.HET,
                "var2": Zygosity.HOM,
                "var3": Zygosity.HOM,
            },
            filtered_out=defaultdict(list),
        )

        # Run the function
        filtered_bg = list(filter_co_existing_subsets({1: bg}).values())[0]

        # The pair should not be removed due to the 'continue' in the branch
        self.assertIn(pair, filtered_bg.alleles[AlleleState.CO])
        self.assertNotIn(pair, filtered_bg.filtered_out["filter_co_existing_subsets"])

    def test_contains_reference_pair_removed(self):
        # Create Reference Allele
        ref_allele = Allele(
            genotype="RefAllele",
            phenotype="RefType",
            genotype_alt=".",
            phenotype_alt=".",
            defining_variants=frozenset({"ref_var"}),
            weight_geno=1,
            reference=True,  # Mark as reference
            sub_type="subtype1",
        )

        # Create Non-Reference Alleles
        allele1 = Allele(
            genotype="Allele1",
            phenotype="Type1",
            genotype_alt=".",
            phenotype_alt=".",
            defining_variants=frozenset({"var1"}),
            weight_geno=1,
            reference=False,
            sub_type="subtype1",
        )
        allele2 = Allele(
            genotype="Allele2",
            phenotype="Type2",
            genotype_alt=".",
            phenotype_alt=".",
            defining_variants=frozenset({"var2"}),
            weight_geno=1,
            reference=False,
            sub_type="subtype1",
        )

        # Create Combined Allele
        combined_allele = Allele(
            genotype="RefAllele+Allele1",
            phenotype="CombinedType",
            genotype_alt=".",
            phenotype_alt=".",
            defining_variants=frozenset({"ref_var", "var1"}),
            weight_geno=1,
            reference=False,
            sub_type="subtype1",
        )

        # Create Pairs
        pair_with_ref = Pair(allele1=ref_allele, allele2=allele1)
        other_pair = Pair(allele1=combined_allele, allele2=allele2)

        # Prepare the BloodGroup object
        bg = BloodGroup(
            type="ExampleType",
            alleles={
                AlleleState.CO: [pair_with_ref, other_pair],
                AlleleState.RAW: [ref_allele, allele1, allele2, combined_allele],
            },
            sample="sample233",
            variant_pool={
                "ref_var": Zygosity.HOM,
                "var1": Zygosity.HOM,
                "var2": Zygosity.HOM,
            },
            filtered_out=defaultdict(list),
        )

        # Run the function
        filtered_bg = list(filter_co_existing_subsets({1: bg}).values())[0]

        # Assert that pair_with_ref is removed
        self.assertNotIn(pair_with_ref, filtered_bg.alleles[AlleleState.CO])
        self.assertIn(
            pair_with_ref, filtered_bg.filtered_out["filter_co_existing_subsets"]
        )

        # Assert that other_pair remains
        self.assertIn(other_pair, filtered_bg.alleles[AlleleState.CO])


class TestFilterPairsOnAntitheticalZygosity(unittest.TestCase):
    def setUp(self):
        allele2 = Allele(
            genotype="FY*02",
            phenotype="FY:2",
            genotype_alt="FY*B",
            phenotype_alt="Fy(b+)",
            defining_variants=frozenset({"1:159175354_G_A"}),
            weight_geno=1000,
            reference=False,
            sub_type="FY*02",
            phases=None,
        )
        allele3 = Allele(
            genotype="FY*01",
            phenotype="FY:1",
            genotype_alt="FY*A",
            phenotype_alt="Fy(a+)",
            defining_variants=frozenset({"1:159175354_ref"}),
            weight_geno=1000,
            reference=True,
            sub_type="FY*01",
            phases=None,
        )
        allele4 = Allele(
            genotype="FY*01N.01",
            phenotype=".",
            genotype_alt=".",
            phenotype_alt="Fy(a-b-)",
            defining_variants=frozenset({"1:159175354_ref", "1:159174683_T_C"}),
            weight_geno=7,
            reference=False,
            sub_type="FY*01",
            phases=None,
        )
        # have to have both subtypes in pair
        self.pair1 = Pair(allele1=allele2, allele2=allele4)  # ok
        self.pair2 = Pair(allele1=allele3, allele2=allele4)  # not ok
        self.bg = BloodGroup(
            type="FY",
            alleles={AlleleState.NORMAL: [self.pair1, self.pair2]},
            sample="013Kenya",
            variant_pool={},
            filtered_out=defaultdict(list),
        )

        self.antitheticals = {
            "KN": ["207782916_A_T", "207782769_G_A,207782916_A_T,207782931_A_G"],
            "LU": ["45315445_G_A", "45315445_ref"],
            "LW": ["10397987_ref", "10397987_A_G"],
            "SC": ["43296522_ref", "43296522_G_A"],
            "YT": ["100490797_ref", "100490797_G_T"],
            "FY": ["159175354_ref", "159175354_G_A"],
        }
        filter_pairs_on_antithetical_zygosity({1: self.bg}, self.antitheticals)

    def test_pairs_removed(self):
        self.assertTrue(
            self.pair2 in self.bg.filtered_out["filter_pairs_on_antithetical_zygosity"]
        )
        self.assertTrue(self.pair2 not in self.bg.alleles[AlleleState.NORMAL])

    def test_no_pairs_removed(self):
        self.assertTrue(
            self.pair1
            not in self.bg.filtered_out["filter_pairs_on_antithetical_zygosity"]
        )
        self.assertTrue(self.pair1 in self.bg.alleles[AlleleState.NORMAL])

    def test_empty_normal_list(self):
        bg = BloodGroup(
            type="FY",
            alleles={AlleleState.NORMAL: []},
            sample="013Kenya",
            variant_pool={},
            filtered_out=defaultdict(list),
        )
        filtered_bg = list(
            filter_pairs_on_antithetical_zygosity({1: bg}, self.antitheticals).values()
        )[0]
        self.assertEqual(filtered_bg.alleles[AlleleState.NORMAL], [])


class TestFilterPairsOnAntitheticalModifyingSNP(unittest.TestCase):
    def setUp(self):
        allele1 = Allele(
            genotype="LU*02",
            phenotype="LU:2",
            genotype_alt="LU*B",
            phenotype_alt="Lu(a-b+)",
            defining_variants=frozenset({"19:45315445_ref"}),
            weight_geno=1,
            reference=True,
            sub_type="LU*02",
        )
        allele2 = Allele(
            genotype="LU*02.19",
            phenotype="LU:-18,19",
            genotype_alt=".",
            phenotype_alt="Au(a-b+)",
            defining_variants=frozenset({"19:45315445_ref", "19:45322744_A_G"}),
            weight_geno=1,
            reference=False,
            sub_type="LU*02",
        )
        allele3 = Allele(
            genotype="LU*01.19",
            phenotype="LU:..",
            genotype_alt=".",
            phenotype_alt="",
            defining_variants=frozenset({"19:45315445_G_A", "19:45322744_A_G"}),
            weight_geno=1,
            reference=False,
            sub_type="LU*01",
        )  #'LU*01.19/LU*02' not possible because modifying SNP (45322744_A_G) is hom
        self.pair1 = Pair(allele1=allele1, allele2=allele2)
        self.pair2 = Pair(allele1=allele1, allele2=allele3)

        self.antitheticals = {
            "KN": ["207782916_A_T", "207782769_G_A,207782916_A_T,207782931_A_G"],
            "LU": ["45315445_G_A", "45315445_ref"],
            "LW": ["10397987_ref", "10397987_A_G"],
            "SC": ["43296522_ref", "43296522_G_A"],
            "YT": ["100490797_ref", "100490797_G_T"],
        }

    def test_pairs_removed_due_to_homozygous_snp(self):
        bg = BloodGroup(
            type="LU",
            alleles={AlleleState.NORMAL: [self.pair1, self.pair2]},
            sample="128",
            variant_pool={
                "19:45315445_G_A": "Heterozygous",
                "19:45315445_ref": "Heterozygous",
                "19:45322744_A_G": "Homozygous",
            },
            filtered_out=defaultdict(list),
        )
        filtered_bg = list(
            filter_pairs_on_antithetical_modyfying_SNP(
                {1: bg}, self.antitheticals
            ).values()
        )[0]
        self.assertTrue(self.pair2 not in filtered_bg.alleles[AlleleState.NORMAL])
        self.assertTrue(
            self.pair2
            in filtered_bg.filtered_out["filter_pairs_on_antithetical_modyfying_SNP"]
        )

    def test_no_pairs_removed_due_to_heterozygous_snp(self):
        bg = BloodGroup(
            type="LU",
            alleles={AlleleState.NORMAL: [self.pair1]},
            sample="128",
            variant_pool={
                "19:45315445_G_A": "Heterozygous",
                "19:45315445_ref": "Heterozygous",
                "19:45322744_A_G": "Homozygous",
            },
            filtered_out=defaultdict(list),
        )
        filtered_bg = list(
            filter_pairs_on_antithetical_modyfying_SNP(
                {1: bg}, self.antitheticals
            ).values()
        )[0]
        self.assertTrue(self.pair1 in filtered_bg.alleles[AlleleState.NORMAL])
        self.assertTrue(
            self.pair1
            not in filtered_bg.filtered_out[
                "filter_pairs_on_antithetical_modyfying_SNP"
            ]
        )

    def test_empty_normal_list(self):
        bg = BloodGroup(
            type="LU",
            alleles={AlleleState.NORMAL: []},
            sample="128",
            variant_pool={},
            filtered_out=defaultdict(list),
        )
        filtered_bg = list(
            filter_pairs_on_antithetical_modyfying_SNP(
                {1: bg}, self.antitheticals
            ).values()
        )[0]
        self.assertEqual(filtered_bg.alleles[AlleleState.NORMAL], [])



class TestCantPairWithRefCuzSNPsMustBeOnOtherSide(unittest.TestCase):
    def setUp(self):
        self.allele1 = Allele(
            genotype="JK*01",
            phenotype="WIP",
            genotype_alt=".",
            phenotype_alt=".",
            defining_variants=frozenset({"18:43319519_ref"}),
            weight_geno=1000,
            reference=True,
            sub_type="JK*01",
        )
        self.allele2 = Allele(
            genotype="JK*01W.03",
            phenotype=".",
            genotype_alt=".",
            phenotype_alt=".",
            defining_variants=frozenset({"18:43310313_G_A"}),
            weight_geno=1000,
            reference=False,
            sub_type="JK*01",
        )
        self.allele3 = Allele(
            genotype="JK*01W.04",
            phenotype=".",
            genotype_alt=".",
            phenotype_alt=".",
            defining_variants=frozenset({"18:43311054_G_A"}),
            weight_geno=1000,
            reference=False,
            sub_type="JK*01",
        )
        self.allele4 = Allele(
            genotype="JK*01W.11",
            phenotype=".",
            genotype_alt=".",
            phenotype_alt=".",
            defining_variants=frozenset({"18:43310313_G_A", "18:43311054_G_A"}),
            weight_geno=1000,
            reference=False,
            sub_type="JK*01",
        )

        self.pair1 = Pair(
            allele1=self.allele3, allele2=self.allele2
        )  # JK*01W.03/4 can be on oposite strands - ok
        self.pair2 = Pair(
            allele1=self.allele1, allele2=self.allele2
        )  # JK*01W.03/4 can't be paired with ref as that means 18:43310313_G_A and
        # 18:43311054_G_A are together, which equals JK*01W.11 - not ok
        self.pair3 = Pair(
            allele1=self.allele1, allele2=self.allele4
        )  # JK*01W.11 can be with ref - ok
        self.bg = BloodGroup(
            type="JK",
            alleles={AlleleState.NORMAL: [self.pair1, self.pair2, self.pair3]},
            sample="003Kenya",
            variant_pool={
                "18:43310313_G_A": "Heterozygous",
                "18:43311054_G_A": "Heterozygous",
                "18:43311131_G_A": "Heterozygous",
                "18:43316538_A_G": "Heterozygous",
            },
            filtered_out=defaultdict(list),
        )
        cant_pair_with_ref_cuz_SNPs_must_be_on_other_side({1: self.bg})

    def test_pairs_removed(self):
        self.assertTrue(self.pair2 not in self.bg.alleles[AlleleState.NORMAL])
        self.assertTrue(
            self.pair2
            in self.bg.filtered_out["cant_pair_with_ref_cuz_SNPs_must_be_on_other_side"]
        )

    def test_pairs_not_removed(self):
        self.assertTrue(self.pair1 in self.bg.alleles[AlleleState.NORMAL])
        self.assertTrue(
            self.pair1
            not in self.bg.filtered_out[
                "cant_pair_with_ref_cuz_SNPs_must_be_on_other_side"
            ]
        )
        self.assertTrue(self.pair3 in self.bg.alleles[AlleleState.NORMAL])
        self.assertTrue(
            self.pair3
            not in self.bg.filtered_out[
                "cant_pair_with_ref_cuz_SNPs_must_be_on_other_side"
            ]
        )


class TestABOCantPairWithRefCuz261delGHET(unittest.TestCase):
    def setUp(self):
        self.allele1 = Allele(
            genotype="ABO*A1.01",
            phenotype=".",
            genotype_alt=".",
            phenotype_alt=".",
            defining_variants=frozenset({"9:136132908_T_TC"}),
            weight_geno=1000,
            reference=True,
            sub_type="ABO*A",
        )
        self.allele2 = Allele(
            genotype="ABO*AW.25",
            phenotype=".",
            genotype_alt=".",
            phenotype_alt=".",
            defining_variants=frozenset(
                {
                    "9:136131056_CG_C",
                    "9:136131289_C_T",
                    "9:136131651_G_A",
                    "9:136132908_T_TC",
                }
            ),
            weight_geno=1000,
            reference=False,
            sub_type="ABO*A",
        )
        self.allele3 = Allele(
            genotype="ABO*O.01.05",
            phenotype=".",
            genotype_alt=".",
            phenotype_alt=".",
            defining_variants=frozenset({"9:136132873_T_C", "9:136132908_ref"}),
            weight_geno=1000,
            reference=False,
            sub_type="ABO*O",
        )

        self.pair1 = Pair(
            allele1=self.allele1, allele2=self.allele2
        )  # Not possible as 136132908_T_TC is in defining vars so need an O - not ok
        self.pair2 = Pair(allele1=self.allele1, allele2=self.allele3)  # ok

        self.bg = BloodGroup(
            type="ABO",
            alleles={AlleleState.NORMAL: [self.pair1, self.pair2]},
            sample="192",
            variant_pool={
                "9:136132908_T_TC": "Heterozygous",
                "9:136132908_ref": "Heterozygous",
                "9:136131056_CG_C": "Heterozygous",
                "9:136131289_C_T": "Heterozygous",
                "9:136131651_G_A": "Heterozygous",
                "9:136132873_T_C": "Heterozygous",
            },
            filtered_out=defaultdict(list),
        )
        ABO_cant_pair_with_ref_cuz_261delG_HET({1: self.bg})

    def test_pairs_removed(self):
        self.assertTrue(self.pair2 in self.bg.alleles[AlleleState.NORMAL])
        self.assertTrue(
            self.pair2
            not in self.bg.filtered_out["ABO_cant_pair_with_ref_cuz_261delG_HET"]
        )

    def test_pairs_not_removed(self):
        self.assertTrue(self.pair1 not in self.bg.alleles[AlleleState.NORMAL])
        self.assertTrue(
            self.pair1 in self.bg.filtered_out["ABO_cant_pair_with_ref_cuz_261delG_HET"]
        )


class TestABOCantPairWithRefCuzTrumped(unittest.TestCase):
    def setUp(self):
        self.allele1 = Allele(
            genotype="FUT3*01",
            phenotype=".",
            genotype_alt=".",
            phenotype_alt=".",
            defining_variants=frozenset({}),
            weight_geno=1000,
            reference=True,
            sub_type="FUT3*01",
        )
        self.allele2 = Allele(
            genotype="FUT3*01.16",
            phenotype=".",
            genotype_alt=".",
            phenotype_alt=".",
            defining_variants=frozenset(
                {"19:5844043_C_T", "19:5844184_C_T", "19:5844367_C_T"}
            ),
            weight_geno=1000,
            reference=False,
            sub_type="FUT3*01",
        )
        self.allele3 = Allele(
            genotype="FUT3*01N.01.02",
            phenotype=".",
            genotype_alt=".",
            phenotype_alt=".",
            defining_variants=frozenset(
                {"19:5844184_C_T", "19:5844367_C_T", "19:5844838_C_T"}
            ),
            weight_geno=1,
            reference=False,
            sub_type="FUT3*01",
        )

        self.pair1 = Pair(
            allele1=self.allele1, allele2=self.allele2
        )  # Not possible - not ok
        self.pair2 = Pair(allele1=self.allele1, allele2=self.allele3)  # ok

        self.bg = BloodGroup(
            type="FUT3",
            alleles={AlleleState.NORMAL: [self.pair1, self.pair2]},
            sample="126",
            variant_pool={
                "19:5843883_C_G": "Heterozygous",
                "19:5844043_C_T": "Heterozygous",
                "19:5844184_C_T": "Heterozygous",
                "19:5844367_C_T": "Homozygous",
                "19:5844838_C_T": "Homozygous",
            },
            filtered_out=defaultdict(list),
        )
        cant_pair_with_ref_cuz_trumped({1: self.bg})

    def test_pairs_not_removed(self):
        self.assertTrue(self.pair2 in self.bg.alleles[AlleleState.NORMAL])
        self.assertTrue(
            self.pair2 not in self.bg.filtered_out["cant_pair_with_ref_cuz_trumped"]
        )

    def test_pairs_removed(self):
        self.assertTrue(self.pair1 not in self.bg.alleles[AlleleState.NORMAL])
        self.assertTrue(
            self.pair1 in self.bg.filtered_out["cant_pair_with_ref_cuz_trumped"]
        )


class TestABOCantPairWithRefCuzTrumped2(unittest.TestCase):
    def setUp(self):
        self.allele1 = Allele(
            genotype="FUT3*01",
            phenotype=".",
            genotype_alt=".",
            phenotype_alt=".",
            defining_variants=frozenset({}),
            weight_geno=1000,
            reference=True,
            sub_type="FUT3*01",
        )
        self.allele2 = Allele(
            genotype="FUT3*01.16",
            phenotype=".",
            genotype_alt=".",
            phenotype_alt=".",
            defining_variants=frozenset(
                {"19:5844043_C_T", "19:5844184_C_T", "19:5844367_C_T"}
            ),
            weight_geno=1000,
            reference=False,
            sub_type="FUT3*01",
        )
        self.allele3 = Allele(
            genotype="FUT3*01N.01.02",
            phenotype=".",
            genotype_alt=".",
            phenotype_alt=".",
            defining_variants=frozenset(
                {"19:5844184_C_T", "19:5844367_C_T", "19:5844838_C_T"}
            ),
            weight_geno=1,
            reference=False,
            sub_type="FUT3*01",
        )

        self.pair1 = Pair(allele1=self.allele1, allele2=self.allele2)  # 2x HET - ok
        self.pair2 = Pair(allele1=self.allele1, allele2=self.allele3)  # ok

        self.bg = BloodGroup(
            type="FUT3",
            alleles={AlleleState.NORMAL: [self.pair1, self.pair2]},
            sample="126",
            variant_pool={
                "19:5843883_C_G": "Heterozygous",
                "19:5844043_C_T": "Heterozygous",
                "19:5844184_C_T": "Heterozygous",
                "19:5844367_C_T": "Heterozygous",
                "19:5844838_C_T": "Homozygous",
            },
            filtered_out=defaultdict(list),
        )
        cant_pair_with_ref_cuz_trumped({1: self.bg})

    def test_pairs_not_removed(self):
        self.assertTrue(self.pair1 in self.bg.alleles[AlleleState.NORMAL])
        self.assertTrue(
            self.pair1 not in self.bg.filtered_out["cant_pair_with_ref_cuz_trumped"]
        )
        self.assertTrue(self.pair2 in self.bg.alleles[AlleleState.NORMAL])
        self.assertTrue(
            self.pair2 not in self.bg.filtered_out["cant_pair_with_ref_cuz_trumped"]
        )




class TestFilterHETPairsByWeight(unittest.TestCase):
    def setUp(self):
        self.allele1 = Allele(
            genotype="FUT2*01",
            genotype_alt=".",
            phenotype=".",
            phenotype_alt=".",
            defining_variants=frozenset({"19:49206250_ref"}),
            weight_geno=1000,
            reference=True,
            sub_type="FUT2*01",
        )
        self.allele2 = Allele(
            genotype="FUT2*01.03.01",
            genotype_alt=".",
            phenotype=".",
            phenotype_alt=".",
            defining_variants=frozenset({"19:49206286_A_G"}),
            weight_geno=1000,
            reference=False,
            sub_type="FUT2*01",
        )
        self.allele3 = Allele(
            genotype="FUT2*01N.02",
            genotype_alt=".",
            phenotype=".",
            phenotype_alt=".",
            defining_variants=frozenset({"19:49206674_G_A"}),
            weight_geno=1,
            reference=False,
            sub_type="FUT2*01",
        )
        self.allele4 = Allele(
            genotype="FUT2*01N.16",
            genotype_alt=".",
            phenotype=".",
            phenotype_alt=".",
            defining_variants=frozenset({"19:49206985_G_A"}),
            weight_geno=8,
            reference=False,
            sub_type="FUT2*01",
        )

        self.pair1 = Pair(allele1=self.allele1, allele2=self.allele2)  # not ok
        self.pair2 = Pair(allele1=self.allele1, allele2=self.allele4)  # not ok
        self.pair3 = Pair(allele1=self.allele1, allele2=self.allele3)  # ok
        self.bg = BloodGroup(
            type="FUT2",
            alleles={AlleleState.NORMAL: [self.pair1, self.pair2, self.pair3]},
            sample="001Kenya",
            variant_pool={
                "19:49206250_ref": "Homozygous",
                "19:49206286_A_G": "Heterozygous",
                "19:49206674_G_A": "Heterozygous",
                "19:49206985_G_A": "Heterozygous",
            },
            filtered_out=defaultdict(list),
        )
        filter_HET_pairs_by_weight({1: self.bg})

    def test_pairs_not_removed(self):
        self.assertTrue(self.pair3 in self.bg.alleles[AlleleState.NORMAL])
        self.assertTrue(
            self.pair3 not in self.bg.filtered_out["filter_HET_pairs_by_weight"]
        )

    def test_pairs_removed(self):
        self.assertTrue(self.pair1 not in self.bg.alleles[AlleleState.NORMAL])
        self.assertTrue(
            self.pair2 in self.bg.filtered_out["filter_HET_pairs_by_weight"]
        )


class TestFilterPairsByContext(unittest.TestCase):
    def setUp(self):
        self.allele1 = Allele(
            genotype="A4GALT*01",
            genotype_alt=".",
            phenotype=".",
            phenotype_alt=".",
            defining_variants=frozenset({"22:43113793_ref"}),
            weight_geno=1000,
            reference=True,
            sub_type="A4GALT*01",
        )
        self.allele2 = Allele(
            genotype="A4GALT*01.02",
            genotype_alt=".",
            phenotype=".",
            phenotype_alt=".",
            defining_variants=frozenset({"22:43089849_T_C"}),
            weight_geno=1000,
            reference=False,
            sub_type="A4GALT*01",
        )
        self.allele3 = Allele(
            genotype="A4GALT*02",
            genotype_alt=".",
            phenotype=".",
            phenotype_alt=".",
            defining_variants=frozenset({"22:43113793_C_A"}),
            weight_geno=1000,
            reference=False,
            sub_type="A4GALT*02",
        )
        self.allele4 = Allele(
            genotype="A4GALT*02.02",
            genotype_alt=".",
            phenotype=".",
            phenotype_alt=".",
            defining_variants=frozenset({"22:43113793_C_A", "22:43089849_T_C"}),
            weight_geno=1000,
            reference=False,
            sub_type="A4GALT*02",
        )

        self.pair1 = Pair(allele1=self.allele1, allele2=self.allele3)  # not ok
        self.pair2 = Pair(allele1=self.allele1, allele2=self.allele4)  # ok
        self.pair3 = Pair(allele1=self.allele2, allele2=self.allele3)
        # not ok (for different reason [antithetical is het])

        self.bg = BloodGroup(
            type="A4GALT",
            alleles={AlleleState.NORMAL: [self.pair1, self.pair2, self.pair3]},
            sample="Kenya",
            variant_pool={
                "22:43089849_T_C": "Heterozygous",
                "22:43113793_C_A": "Heterozygous",
                "22:43113793_ref": "Heterozygous",
            },
            filtered_out=defaultdict(list),
        )
        filter_pairs_by_context({1: self.bg})

    # def test_pairs_not_removed(self): TODO
    #     self.assertTrue(self.pair2 in self.bg.alleles[AlleleState.NORMAL])
    #     self.assertTrue(
    #         self.pair2 not in self.bg.filtered_out["filter_pairs_by_context"]
    #     )

    # def test_pairs_removed(self):
    #     self.assertTrue(self.pair1 not in self.bg.alleles[AlleleState.NORMAL])
    #     self.assertTrue(
    #         self.pair1 in self.bg.filtered_out["filter_pairs_by_context"]
    #     )


class TestFilterCoexistingPairsOnAntitheticalZygosity(unittest.TestCase):
    def test_no_co_alleles(self):
        """Test when bg.alleles[AlleleState.CO] is None."""
        bg = BloodGroup(
            type="ExampleType",
            alleles={AlleleState.CO: None},
            sample="sample1",
        )
        antitheticals = {"ExampleType": ["antigen1", "antigen2"]}

        # result_bg = filter_coexisting_pairs_on_antithetical_zygosity(bg, antitheticals)
        result_bg = list(
            filter_coexisting_pairs_on_antithetical_zygosity(
                {1: bg}, antitheticals
            ).values()
        )[0]

        self.assertIs(result_bg, bg)
        self.assertIsNone(result_bg.alleles[AlleleState.CO])

    def test_bg_type_not_in_antitheticals(self):
        """Test when bg.type is not in antitheticals."""
        allele1 = Allele(
            genotype="Allele1",
            phenotype="Phenotype1",
            genotype_alt=".",
            phenotype_alt=".",
            defining_variants=frozenset(),
            weight_geno=1,
            reference=False,
            sub_type="Subtype1",
        )
        allele2 = Allele(
            genotype="Allele2",
            phenotype="Phenotype2",
            genotype_alt=".",
            phenotype_alt=".",
            defining_variants=frozenset(),
            weight_geno=1,
            reference=False,
            sub_type="Subtype2",
        )
        pair = Pair(allele1=allele1, allele2=allele2)

        bg = BloodGroup(
            type="OtherType",
            alleles={AlleleState.CO: [pair], AlleleState.NORMAL: [allele1, allele2]},
            sample="sample2",
        )
        antitheticals = {"ExampleType": ["antigen1", "antigen2"]}

        # result_bg = filter_coexisting_pairs_on_antithetical_zygosity(bg, antitheticals)
        result_bg = list(
            filter_coexisting_pairs_on_antithetical_zygosity(
                {1: bg}, antitheticals
            ).values()
        )[0]
        self.assertIs(result_bg, bg)
        self.assertIn(pair, result_bg.alleles[AlleleState.CO])

    def test_flattened_sub_types_length_one(self):
        """Test when length of flattened_sub_types is 1."""
        allele1 = Allele(
            genotype="Allele1",
            phenotype="Phenotype1",
            genotype_alt=".",
            phenotype_alt=".",
            defining_variants=frozenset(),
            weight_geno=1,
            reference=False,
            sub_type="Subtype1",
        )
        allele2 = Allele(
            genotype="Allele2",
            phenotype="Phenotype2",
            genotype_alt=".",
            phenotype_alt=".",
            defining_variants=frozenset(),
            weight_geno=1,
            reference=False,
            sub_type="Subtype1",
        )
        pair = Pair(allele1=allele1, allele2=allele2)

        # For AlleleState.NORMAL, we use Pair objects
        pair_normal = Pair(allele1=allele1, allele2=allele2)

        bg = BloodGroup(
            type="ExampleType",
            alleles={
                AlleleState.CO: [pair],
                AlleleState.NORMAL: [pair_normal],  # Now contains Pair objects
            },
            sample="sample3",
        )
        antitheticals = {"ExampleType": ["antigen1", "antigen2"]}

        result_bg = list(
            filter_coexisting_pairs_on_antithetical_zygosity(
                {1: bg}, antitheticals
            ).values()
        )[0]

        self.assertIs(result_bg, bg)
        self.assertIn(pair, result_bg.alleles[AlleleState.CO])

    def test_pairs_with_flat_sub_types_equal(self):
        """Test when flat_sub_types == flattened_sub_types."""
        allele1 = Allele(
            genotype="Allele1",
            phenotype="Phenotype1",
            genotype_alt=".",
            phenotype_alt=".",
            defining_variants=frozenset(),
            weight_geno=1,
            reference=False,
            sub_type="Subtype1",
        )
        allele2 = Allele(
            genotype="Allele2",
            phenotype="Phenotype2",
            genotype_alt=".",
            phenotype_alt=".",
            defining_variants=frozenset(),
            weight_geno=1,
            reference=False,
            sub_type="Subtype2",
        )
        pair = Pair(allele1=allele1, allele2=allele2)

        # For AlleleState.NORMAL, we use Pair objects
        pair_normal = Pair(allele1=allele1, allele2=allele2)

        bg = BloodGroup(
            type="ExampleType",
            alleles={
                AlleleState.CO: [pair],
                AlleleState.NORMAL: [pair_normal],  # Now contains Pair objects
            },
            sample="sample4",
        )
        antitheticals = {"ExampleType": ["antigen1", "antigen2"]}

        result_bg = list(
            filter_coexisting_pairs_on_antithetical_zygosity(
                {1: bg}, antitheticals
            ).values()
        )[0]

        self.assertIs(result_bg, bg)
        self.assertIn(pair, result_bg.alleles[AlleleState.CO])

    def test_pairs_with_flat_sub_types_not_equal(self):
        """Test when flat_sub_types != flattened_sub_types."""
        allele1 = Allele(
            genotype="Allele1",
            phenotype="Phenotype1",
            genotype_alt=".",
            phenotype_alt=".",
            defining_variants=frozenset(),
            weight_geno=1,
            reference=False,
            sub_type="Subtype1",
        )
        allele2 = Allele(
            genotype="Allele2",
            phenotype="Phenotype2",
            genotype_alt=".",
            phenotype_alt=".",
            defining_variants=frozenset(),
            weight_geno=1,
            reference=False,
            sub_type="Subtype1",
        )
        allele3 = Allele(
            genotype="Allele3",
            phenotype="Phenotype3",
            genotype_alt=".",
            phenotype_alt=".",
            defining_variants=frozenset(),
            weight_geno=1,
            reference=False,
            sub_type="Subtype2",
        )
        # Create pairs for NORMAL state
        pair_normal1 = Pair(allele1=allele1, allele2=allele2)
        pair_normal2 = Pair(allele1=allele1, allele2=allele3)  # Different subtypes

        pair_co = Pair(allele1=allele1, allele2=allele2)

        bg = BloodGroup(
            type="ExampleType",
            alleles={
                AlleleState.CO: [pair_co],
                AlleleState.NORMAL: [pair_normal1, pair_normal2],
            },
            sample="sample5",
            filtered_out=defaultdict(list),
        )
        antitheticals = {"ExampleType": ["antigen1", "antigen2"]}

        # Adjusted function call to match the decorator
        result_bg = list(
            filter_coexisting_pairs_on_antithetical_zygosity(
                {1: bg}, antitheticals
            ).values()
        )[0]

        self.assertIs(result_bg, bg)
        self.assertNotIn(pair_co, result_bg.alleles[AlleleState.CO])
        self.assertIn(
            pair_co, result_bg.filtered_out["filter_co_pairs_on_antithetical_zygosity"]
        )


class TestEnsureCoExistingHetSnpUsed(unittest.TestCase):
    def test_no_co_alleles(self):
        """Test when bg.alleles[AlleleState.CO] is None."""
        bg = BloodGroup(
            type="ExampleType",
            alleles={AlleleState.CO: None},
            sample="sample1",
            variant_pool={"variant1": Zygosity.HET},
            filtered_out=defaultdict(list),
        )
        result_bg = list(ensure_co_existing_HET_SNP_used({1: bg}).values())[0]
        self.assertIs(result_bg, bg)
        self.assertIsNone(result_bg.alleles[AlleleState.CO])

    def test_no_het_variants(self):
        """Test when there are no heterozygous variants in variant_pool."""
        allele1 = Allele(
            genotype="Allele1",
            phenotype="Phenotype1",
            genotype_alt=".",
            phenotype_alt=".",
            defining_variants=frozenset({"var1"}),
            weight_geno=1,
            reference=False,
            sub_type="Subtype1",
        )
        allele2 = Allele(
            genotype="Allele2",
            phenotype="Phenotype2",
            genotype_alt=".",
            phenotype_alt=".",
            defining_variants=frozenset({"var2"}),
            weight_geno=1,
            reference=False,
            sub_type="Subtype1",
        )
        pair = Pair(allele1=allele1, allele2=allele2)
        bg = BloodGroup(
            type="ExampleType",
            alleles={AlleleState.CO: [pair]},
            sample="sample2",
            variant_pool={"variant3": Zygosity.HOM},
            filtered_out=defaultdict(list),
        )
        result_bg = list(ensure_co_existing_HET_SNP_used({1: bg}).values())[0]
        self.assertIs(result_bg, bg)
        self.assertIn(pair, result_bg.alleles[AlleleState.CO])
        self.assertNotIn(pair, result_bg.filtered_out["ensure_HET_SNP_used_CO"])

    def test_het_variants_no_matching_alleles(self):
        """Test when there are heterozygous variants but no matching alleles are found."""
        allele1 = Allele(
            genotype="Allele1",
            phenotype="Phenotype1",
            genotype_alt=".",
            phenotype_alt=".",
            defining_variants=frozenset({"var1"}),
            weight_geno=1,
            reference=False,
            sub_type="Subtype1",
        )
        allele2 = Allele(
            genotype="Allele2",
            phenotype="Phenotype2",
            genotype_alt=".",
            phenotype_alt=".",
            defining_variants=frozenset({"var2"}),
            weight_geno=1,
            reference=False,
            sub_type="Subtype1",
        )
        pair = Pair(allele1=allele1, allele2=allele2)
        allele3 = Allele(
            genotype="Allele3",
            phenotype="Phenotype3",
            genotype_alt=".",
            phenotype_alt=".",
            defining_variants=frozenset({"var3"}),
            weight_geno=1,
            reference=False,
            sub_type="Subtype1",
        )
        pair_other = Pair(allele1=allele3, allele2=allele3)
        bg = BloodGroup(
            type="ExampleType",
            alleles={AlleleState.CO: [pair, pair_other]},
            sample="sample3",
            variant_pool={"variant1": Zygosity.HET},
            filtered_out=defaultdict(list),
        )
        result_bg = list(ensure_co_existing_HET_SNP_used({1: bg}).values())[0]
        self.assertIs(result_bg, bg)
        self.assertIn(pair, result_bg.alleles[AlleleState.CO])
        self.assertIn(pair_other, result_bg.alleles[AlleleState.CO])
        self.assertNotIn(pair, result_bg.filtered_out["ensure_HET_SNP_used_CO"])

    def test_hits_not_greater_than_one(self):
        """Test when hits are not greater than 1, pairs are not removed."""
        allele1 = Allele(
            genotype="Allele1",
            phenotype="Phenotype1",
            genotype_alt=".",
            phenotype_alt=".",
            defining_variants=frozenset({"var1"}),
            weight_geno=1,
            reference=False,
            sub_type="Subtype1",
        )
        allele2 = Allele(
            genotype="Allele2",
            phenotype="Phenotype2",
            genotype_alt=".",
            phenotype_alt=".",
            defining_variants=frozenset({"var2"}),
            weight_geno=1,
            reference=False,
            sub_type="Subtype1",
        )
        pair = Pair(allele1=allele1, allele2=allele2)
        allele3 = Allele(
            genotype="Allele3",
            phenotype="Phenotype3",
            genotype_alt=".",
            phenotype_alt=".",
            defining_variants=frozenset({"var1", "variant1"}),
            weight_geno=1,
            reference=False,
            sub_type="Subtype1",
        )
        pair_other = Pair(allele1=allele3, allele2=allele3)
        bg = BloodGroup(
            type="ExampleType",
            alleles={AlleleState.CO: [pair, pair_other]},
            sample="sample4",
            variant_pool={"variant1": Zygosity.HET},
            filtered_out=defaultdict(list),
        )
        result_bg = list(ensure_co_existing_HET_SNP_used({1: bg}).values())[0]
        self.assertIs(result_bg, bg)
        self.assertIn(pair, result_bg.alleles[AlleleState.CO])
        self.assertNotIn(pair, result_bg.filtered_out["ensure_HET_SNP_used_CO"])

    def test_hits_greater_than_one_pair_removed(self):
        """Test when hits > 1, pairs are removed."""
        allele1 = Allele(
            genotype="Allele1",
            phenotype="Phenotype1",
            genotype_alt=".",
            phenotype_alt=".",
            defining_variants=frozenset({"var1"}),
            weight_geno=1,
            reference=False,
            sub_type="Subtype1",
        )
        allele2 = Allele(
            genotype="Allele2",
            phenotype="Phenotype2",
            genotype_alt=".",
            phenotype_alt=".",
            defining_variants=frozenset({"var2"}),
            weight_geno=1,
            reference=False,
            sub_type="Subtype1",
        )
        pair = Pair(allele1=allele1, allele2=allele2)
        allele3 = Allele(
            genotype="Allele3",
            phenotype="Phenotype3",
            genotype_alt=".",
            phenotype_alt=".",
            defining_variants=frozenset({"var1", "variant1"}),
            weight_geno=1,
            reference=False,
            sub_type="Subtype1",
        )
        allele4 = Allele(
            genotype="Allele4",
            phenotype="Phenotype4",
            genotype_alt=".",
            phenotype_alt=".",
            defining_variants=frozenset({"var2", "variant1"}),
            weight_geno=1,
            reference=False,
            sub_type="Subtype1",
        )
        pair_other1 = Pair(allele1=allele3, allele2=allele3)
        pair_other2 = Pair(allele1=allele4, allele2=allele4)
        bg = BloodGroup(
            type="ExampleType",
            alleles={AlleleState.CO: [pair, pair_other1, pair_other2]},
            sample="sample5",
            variant_pool={"variant1": Zygosity.HET},
            filtered_out=defaultdict(list),
        )
        result_bg = list(ensure_co_existing_HET_SNP_used({1: bg}).values())[0]
        self.assertIs(result_bg, bg)
        self.assertNotIn(pair, result_bg.alleles[AlleleState.CO])
        self.assertIn(pair, result_bg.filtered_out["ensure_co_existing_HET_SNP_used"])
        self.assertIn(pair_other1, result_bg.alleles[AlleleState.CO])
        self.assertIn(pair_other2, result_bg.alleles[AlleleState.CO])


# class TestFilterPairsByContext(unittest.TestCase):
#     def test_no_normal_alleles(self):
#         """Test when bg.alleles[AlleleState.NORMAL] is None."""
#         bg = BloodGroup(
#             type="ExampleType",
#             alleles={AlleleState.NORMAL: None},
#             sample="sample1",
#             variant_pool={},
#             filtered_out=defaultdict(list),
#         )
#         result_bg = list(filter_pairs_by_context({1: bg}).values())[0]
#         self.assertIs(result_bg, bg)
#         self.assertIsNone(result_bg.alleles[AlleleState.NORMAL])

#     def test_allele_reference(self):
#         """Test when an allele in the pair is a reference allele."""
#         # Create reference allele
#         ref_allele = Allele(
#             genotype="RefAllele",
#             phenotype="RefPhenotype",
#             genotype_alt=".",
#             phenotype_alt=".",
#             defining_variants=frozenset({"var_ref"}),
#             weight_geno=1,
#             reference=True,
#             sub_type="Subtype1",
#         )
#         # Create another allele
#         allele = Allele(
#             genotype="Allele1",
#             phenotype="Phenotype1",
#             genotype_alt=".",
#             phenotype_alt=".",
#             defining_variants=frozenset({"var1"}),
#             weight_geno=1,
#             weight_pheno=1,
#             reference=False,
#             sub_type="Subtype1",
#         )
#         pair = Pair(allele1=ref_allele, allele2=allele)
#         bg = BloodGroup(
#             type="ExampleType",
#             alleles={AlleleState.NORMAL: [pair]},
#             sample="sample2",
#             variant_pool={"var_ref": Zygosity.HOM, "var1": Zygosity.HET},
#             filtered_out=defaultdict(list),
#         )
#         result_bg = list(filter_pairs_by_context({1: bg}).values())[0]
#         # Since allele.reference is True, it should skip processing that allele
#         # and not remove the pair
#         self.assertIn(pair, result_bg.alleles[AlleleState.NORMAL])
#         self.assertNotIn(pair, result_bg.filtered_out["filter_pairs_by_context"])

#     def test_left_over_vars_empty_or_def_vars_small(self):
#         """Test when len(left_over_vars) == 0 or len(def_vars) < 2."""
#         # Create two alleles
#         allele1 = Allele(
#             genotype="Allele1",
#             phenotype="Phenotype1",
#             genotype_alt=".",
#             phenotype_alt=".",
#             defining_variants=frozenset({"var1"}),
#             weight_geno=1,
#             weight_pheno=1,
#             reference=False,
#             sub_type="Subtype1",
#         )
#         allele2 = Allele(
#             genotype="Allele2",
#             phenotype="Phenotype2",
#             genotype_alt=".",
#             phenotype_alt=".",
#             defining_variants=frozenset({"var2"}),
#             weight_geno=1,
#             weight_pheno=1,
#             reference=False,
#             sub_type="Subtype1",
#         )
#         pair = Pair(allele1=allele1, allele2=allele2)
#         # def_vars will be empty since there are no other pairs
#         bg = BloodGroup(
#             type="ExampleType",
#             alleles={AlleleState.NORMAL: [pair]},
#             sample="sample3",
#             variant_pool={"var1": Zygosity.HOM, "var2": Zygosity.HOM},
#             filtered_out=defaultdict(list),
#         )
#         result_bg = list(filter_pairs_by_context({1: bg}).values())[0]
#         # Since len(def_vars) < 2, it should skip and not remove the pair
#         self.assertIn(pair, result_bg.alleles[AlleleState.NORMAL])
#         self.assertNotIn(pair, result_bg.filtered_out["filter_pairs_by_context"])

    
#     def test_pair_not_removed(self):
#         """Test when the condition 'if all(variants in def_vars for variants in remaining)' is False."""
#         # Create alleles
#         allele1 = Allele(
#             genotype="Allele1",
#             phenotype="Phenotype1",
#             genotype_alt=".",
#             phenotype_alt=".",
#             defining_variants=frozenset({"var1"}),
#             weight_geno=1,
#             weight_pheno=1,
#             reference=False,
#             sub_type="Subtype1",
#         )
#         allele2 = Allele(
#             genotype="Allele2",
#             phenotype="Phenotype2",
#             genotype_alt=".",
#             phenotype_alt=".",
#             defining_variants=frozenset({"var2"}),
#             weight_geno=1,
#             weight_pheno=1,
#             reference=False,
#             sub_type="Subtype1",
#         )
#         # Other allele with different variants
#         other_allele = Allele(
#             genotype="OtherAllele",
#             phenotype="OtherPhenotype",
#             genotype_alt=".",
#             phenotype_alt=".",
#             defining_variants=frozenset({"var3"}),
#             weight_geno=1,
#             weight_pheno=1,
#             reference=False,
#             sub_type="Subtype1",
#         )
#         # Pairs
#         pair = Pair(allele1=allele1, allele2=allele2)
#         other_pair = Pair(allele1=other_allele, allele2=allele2)
#         bg = BloodGroup(
#             type="ExampleType",
#             alleles={AlleleState.NORMAL: [pair, other_pair]},
#             sample="sample5",
#             variant_pool={
#                 "var1": Zygosity.HET,
#                 "var2": Zygosity.HET,
#                 "var3": Zygosity.HET,
#             },
#             filtered_out=defaultdict(list),
#         )
#         result_bg = list(filter_pairs_by_context({1: bg}).values())[0]
#         # The pair should not be removed
#         self.assertIn(pair, result_bg.alleles[AlleleState.NORMAL])
#         self.assertNotIn(pair, result_bg.filtered_out["filter_pairs_by_context"])


class TestFilterPairsByPhase(unittest.TestCase):
    def setUp(self):
        # Define a reference allele for testing
        self.reference_allele = Allele(
            genotype="RefAllele",
            phenotype="RefPhenotype",
            genotype_alt=".",
            phenotype_alt=".",
            defining_variants=frozenset(),
            weight_geno=1,
            reference=True,
            sub_type="SubtypeRef",
            phases=(".",),
        )
        self.reference_alleles = {"ExampleType": self.reference_allele}

    def test_phased_false(self):
        """Test when phased is False."""
        allele1 = Allele(
            genotype="Allele1",
            phenotype="Phenotype1",
            genotype_alt=".",
            phenotype_alt=".",
            defining_variants=frozenset(),
            weight_geno=1,
            reference=False,
            sub_type="Subtype1",
            phases=("1",),
        )
        allele2 = Allele(
            genotype="Allele2",
            phenotype="Phenotype2",
            genotype_alt=".",
            phenotype_alt=".",
            defining_variants=frozenset(),
            weight_geno=1,
            reference=False,
            sub_type="Subtype2",
            phases=("2",),
        )
        pair = Pair(allele1=allele1, allele2=allele2)
        bg = BloodGroup(
            type="ExampleType",
            alleles={AlleleState.NORMAL: [pair]},
            sample="sample1",
            filtered_out=defaultdict(list),
        )
        result_bg = list(
            filter_pairs_by_phase({1: bg}, False, self.reference_alleles).values()
        )[0]
        self.assertIs(result_bg, bg)
        self.assertIn(pair, result_bg.alleles[AlleleState.NORMAL])
        self.assertEqual(len(result_bg.filtered_out), 0)

    def test_pair_contains_reference(self):
        """Test when pair contains a reference allele."""
        allele1 = self.reference_allele
        allele2 = Allele(
            genotype="Allele2",
            phenotype="Phenotype2",
            genotype_alt=".",
            phenotype_alt=".",
            defining_variants=frozenset(),
            weight_geno=1,
            reference=False,
            sub_type="Subtype2",
            phases=("1",),
        )
        pair = Pair(allele1=allele1, allele2=allele2)
        bg = BloodGroup(
            type="ExampleType",
            alleles={AlleleState.NORMAL: [pair]},
            sample="sample2",
            filtered_out=defaultdict(list),
        )
        result_bg = list(
            filter_pairs_by_phase({1: bg}, True, self.reference_alleles).values()
        )[0]
        self.assertIn(pair, result_bg.alleles[AlleleState.NORMAL])
        self.assertEqual(len(result_bg.filtered_out), 0)

    def test_both_homozygous(self):
        """Test when both alleles are homozygous (phase sets are {'.'})."""
        allele1 = Allele(
            genotype="Allele1",
            phenotype="Phenotype1",
            genotype_alt=".",
            phenotype_alt=".",
            defining_variants=frozenset(),
            weight_geno=1,
            reference=False,
            sub_type="Subtype1",
            phases=(".",),
        )
        allele2 = Allele(
            genotype="Allele2",
            phenotype="Phenotype2",
            genotype_alt=".",
            phenotype_alt=".",
            defining_variants=frozenset(),
            weight_geno=1,
            reference=False,
            sub_type="Subtype2",
            phases=(".",),
        )
        pair = Pair(allele1=allele1, allele2=allele2)
        bg = BloodGroup(
            type="ExampleType",
            alleles={AlleleState.NORMAL: [pair]},
            sample="sample3",
            filtered_out=defaultdict(list),
        )
        result_bg = list(
            filter_pairs_by_phase({1: bg}, True, self.reference_alleles).values()
        )[0]
        self.assertIn(pair, result_bg.alleles[AlleleState.NORMAL])
        self.assertEqual(len(result_bg.filtered_out), 0)

    def test_same_phase_sets_removed(self):
        """Test when both alleles have the same phase sets and are removed."""
        allele1 = Allele(
            genotype="Allele1",
            phenotype="Phenotype1",
            genotype_alt=".",
            phenotype_alt=".",
            defining_variants=frozenset(),
            weight_geno=1,
            reference=False,
            sub_type="Subtype1",
            phases=("1", "2"),
        )
        allele2 = Allele(
            genotype="Allele2",
            phenotype="Phenotype2",
            genotype_alt=".",
            phenotype_alt=".",
            defining_variants=frozenset(),
            weight_geno=1,
            reference=False,
            sub_type="Subtype2",
            phases=("1", "2"),
        )
        pair = Pair(allele1=allele1, allele2=allele2)
        bg = BloodGroup(
            type="ExampleType",
            alleles={AlleleState.NORMAL: [pair]},
            sample="sample4",
            filtered_out=defaultdict(list),
        )
        result_bg = list(
            filter_pairs_by_phase({1: bg}, True, self.reference_alleles).values()
        )[0]
        self.assertNotIn(pair, result_bg.alleles[AlleleState.NORMAL])
        self.assertIn(pair, result_bg.filtered_out["filter_pairs_by_phase"])

    def test_different_phase_sets_retained(self):
        """Test when alleles have different phase sets and are retained."""
        allele1 = Allele(
            genotype="Allele1",
            phenotype="Phenotype1",
            genotype_alt=".",
            phenotype_alt=".",
            defining_variants=frozenset(),
            weight_geno=1,
            reference=False,
            sub_type="Subtype1",
            phases=("1",),
        )
        allele2 = Allele(
            genotype="Allele2",
            phenotype="Phenotype2",
            genotype_alt=".",
            phenotype_alt=".",
            defining_variants=frozenset(),
            weight_geno=1,
            reference=False,
            sub_type="Subtype2",
            phases=("2",),
        )
        pair = Pair(allele1=allele1, allele2=allele2)
        bg = BloodGroup(
            type="ExampleType",
            alleles={AlleleState.NORMAL: [pair]},
            sample="sample5",
            filtered_out=defaultdict(list),
        )
        result_bg = list(
            filter_pairs_by_phase({1: bg}, True, self.reference_alleles).values()
        )[0]
        self.assertIn(pair, result_bg.alleles[AlleleState.NORMAL])
        self.assertEqual(len(result_bg.filtered_out), 0)

    def test_remove_homs_different(self):
        """Test when remove_homs(p1) != remove_homs(p2) and pair is retained."""
        allele1 = Allele(
            genotype="Allele1",
            phenotype="Phenotype1",
            genotype_alt=".",
            phenotype_alt=".",
            defining_variants=frozenset(),
            weight_geno=1,
            reference=False,
            sub_type="Subtype1",
            phases=("1", "."),
        )
        allele2 = Allele(
            genotype="Allele2",
            phenotype="Phenotype2",
            genotype_alt=".",
            phenotype_alt=".",
            defining_variants=frozenset(),
            weight_geno=1,
            reference=False,
            sub_type="Subtype2",
            phases=("2", "."),
        )
        pair = Pair(allele1=allele1, allele2=allele2)
        bg = BloodGroup(
            type="ExampleType",
            alleles={AlleleState.NORMAL: [pair]},
            sample="sample6",
            filtered_out=defaultdict(list),
        )
        result_bg = list(
            filter_pairs_by_phase({1: bg}, True, self.reference_alleles).values()
        )[0]
        self.assertIn(pair, result_bg.alleles[AlleleState.NORMAL])
        self.assertEqual(len(result_bg.filtered_out), 0)

    def test_remove_homs_same_removed(self):
        """Test when remove_homs(p1) == remove_homs(p2) and pair is removed."""
        allele1 = Allele(
            genotype="Allele1",
            phenotype="Phenotype1",
            genotype_alt=".",
            phenotype_alt=".",
            defining_variants=frozenset(),
            weight_geno=1,
            reference=False,
            sub_type="Subtype1",
            phases=("1", "."),
        )
        allele2 = Allele(
            genotype="Allele2",
            phenotype="Phenotype2",
            genotype_alt=".",
            phenotype_alt=".",
            defining_variants=frozenset(),
            weight_geno=1,
            reference=False,
            sub_type="Subtype2",
            phases=("1", "."),
        )
        pair = Pair(allele1=allele1, allele2=allele2)
        bg = BloodGroup(
            type="ExampleType",
            alleles={AlleleState.NORMAL: [pair]},
            sample="sample7",
            filtered_out=defaultdict(list),
        )
        result_bg = list(
            filter_pairs_by_phase({1: bg}, True, self.reference_alleles).values()
        )[0]
        self.assertNotIn(pair, result_bg.alleles[AlleleState.NORMAL])
        self.assertIn(pair, result_bg.filtered_out["filter_pairs_by_phase"])

    def test_all_pairs_removed_add_reference(self):
        """Test when all pairs are removed and new pairs with reference allele are added."""
        allele1 = Allele(
            genotype="Allele1",
            phenotype="Phenotype1",
            genotype_alt=".",
            phenotype_alt=".",
            defining_variants=frozenset(),
            weight_geno=1,
            reference=False,
            sub_type="Subtype1",
            phases=("1",),
        )
        allele2 = Allele(
            genotype="Allele2",
            phenotype="Phenotype2",
            genotype_alt=".",
            phenotype_alt=".",
            defining_variants=frozenset(),
            weight_geno=1,
            reference=False,
            sub_type="Subtype2",
            phases=("1",),
        )
        pair = Pair(allele1=allele1, allele2=allele2)
        bg = BloodGroup(
            type="ExampleType",
            alleles={AlleleState.NORMAL: [pair]},
            sample="sample8",
            filtered_out=defaultdict(list),
        )
        result_bg = list(
            filter_pairs_by_phase({1: bg}, True, self.reference_alleles).values()
        )[0]
        # Original pair should be removed
        self.assertNotIn(pair, result_bg.alleles[AlleleState.NORMAL])
        self.assertIn(pair, result_bg.filtered_out["filter_pairs_by_phase"])
        # New pairs with reference allele should be added
        new_pair1 = Pair(self.reference_allele, allele1)
        new_pair2 = Pair(self.reference_allele, allele2)
        self.assertIn(new_pair1, result_bg.alleles[AlleleState.NORMAL])
        self.assertIn(new_pair2, result_bg.alleles[AlleleState.NORMAL])


if __name__ == "__main__":
    unittest.main()
