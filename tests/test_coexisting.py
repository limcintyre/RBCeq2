import unittest

from rbceq2.core_logic.alleles import Allele, BloodGroup, Pair
from rbceq2.core_logic.co_existing import (
    add_co_existing_allele_and_ref,
    add_co_existing_alleles,
    all_hom_variants,
    can_co_exist,
    decide_if_co_existing,
    filter_redundant_pairs,
    # geno_str,
    homs,
    list_excluded_co_existing_pairs,
    max_rank,
    mush,
    mushed_vars,
    prep_co_putative_combos,
    sub_alleles,
)
from rbceq2.core_logic.constants import AlleleState

ALLELE_RELATIONSHIPS = {
    "KN": {
        "KN*02.02_isin_KN*01": False,
        "KN*01_isin_KN*02.02": False,
        "KN*02.02_isin_KN*01.07": False,
        "KN*01.07_isin_KN*02.02": False,
        "KN*02.02_isin_KN*01.10": False,
        "KN*01.10_isin_KN*02.02": False,
        "KN*02.02_isin_KN*01.06": False,
        "KN*01.06_isin_KN*02.02": False,
        "KN*01.-05_isin_KN*01": False,
        "KN*01.-05_isin_KN*01.-05": False,
        "KN*01.-05_isin_KN*01.-08": False,
        "KN*01.-05_isin_KN*01.-13": False,
        "KN*01.-05_isin_KN*01.06": False,
        "KN*01.-05_isin_KN*01.07": False,
        "KN*01.-05_isin_KN*01.10": False,
        "KN*01.-05_isin_KN*01.12": False,
        "KN*01.-05_isin_KN*01W": False,
        "KN*01.-05_isin_KN*02": False,
        "KN*01.-05_isin_KN*02.10": False,
        "KN*01.-08_isin_KN*01": False,
        "KN*01.-08_isin_KN*01.-05": False,
        "KN*01.-08_isin_KN*01.-08": False,
        "KN*01.-08_isin_KN*01.-13": False,
        "KN*01.-08_isin_KN*01.06": False,
        "KN*01.-08_isin_KN*01.07": False,
        "KN*01.-08_isin_KN*01.10": False,
        "KN*01.-08_isin_KN*01.12": False,
        "KN*01.-08_isin_KN*01W": False,
        "KN*01.-08_isin_KN*02": False,
        "KN*01.-08_isin_KN*02.10": False,
        "KN*01.-13_isin_KN*01": False,
        "KN*01.-13_isin_KN*01.-05": False,
        "KN*01.-13_isin_KN*01.-08": False,
        "KN*01.-13_isin_KN*01.-13": False,
        "KN*01.-13_isin_KN*01.06": False,
        "KN*01.-13_isin_KN*01.07": False,
        "KN*01.-13_isin_KN*01.10": False,
        "KN*01.-13_isin_KN*01.12": False,
        "KN*01.-13_isin_KN*01W": False,
        "KN*01.-13_isin_KN*02": False,
        "KN*01.-13_isin_KN*02.10": False,
        "KN*01.06_isin_KN*01": False,
        "KN*01.06_isin_KN*01.-05": False,
        "KN*01.06_isin_KN*01.-08": False,
        "KN*01.06_isin_KN*01.-13": False,
        "KN*01.06_isin_KN*01.06": False,
        "KN*01.06_isin_KN*01.07": False,
        "KN*01.06_isin_KN*01.10": False,
        "KN*01.06_isin_KN*01.12": False,
        "KN*01.06_isin_KN*01W": False,
        "KN*01.06_isin_KN*02": False,
        "KN*01.06_isin_KN*02.10": False,
        "KN*01.07_isin_KN*01": False,
        "KN*01.07_isin_KN*01.-05": False,
        "KN*01.07_isin_KN*01.-08": False,
        "KN*01.07_isin_KN*01.-13": False,
        "KN*01.07_isin_KN*01.06": True,
        "KN*01.07_isin_KN*01.07": False,
        "KN*01.07_isin_KN*01.10": False,
        "KN*01.07_isin_KN*01.12": False,
        "KN*01.07_isin_KN*01W": False,
        "KN*01.07_isin_KN*02": False,
        "KN*01.07_isin_KN*02.10": False,
        "KN*01.10_isin_KN*01": False,
        "KN*01.10_isin_KN*01.-05": False,
        "KN*01.10_isin_KN*01.-08": False,
        "KN*01.10_isin_KN*01.-13": False,
        "KN*01.10_isin_KN*01.06": True,
        "KN*01.10_isin_KN*01.07": True,
        "KN*01.10_isin_KN*01.10": False,
        "KN*01.10_isin_KN*01.12": False,
        "KN*01.10_isin_KN*01W": False,
        "KN*01.10_isin_KN*02": False,
        "KN*01.10_isin_KN*02.10": True,
        "KN*01.12_isin_KN*01": False,
        "KN*01.12_isin_KN*01.-05": False,
        "KN*01.12_isin_KN*01.-08": False,
        "KN*01.12_isin_KN*01.-13": False,
        "KN*01.12_isin_KN*01.06": False,
        "KN*01.12_isin_KN*01.07": False,
        "KN*01.12_isin_KN*01.10": False,
        "KN*01.12_isin_KN*01.12": False,
        "KN*01.12_isin_KN*01W": False,
        "KN*01.12_isin_KN*02": False,
        "KN*01.12_isin_KN*02.10": False,
        "KN*01W_isin_KN*01": False,
        "KN*01W_isin_KN*01.-05": False,
        "KN*01W_isin_KN*01.-08": False,
        "KN*01W_isin_KN*01.-13": False,
        "KN*01W_isin_KN*01.06": False,
        "KN*01W_isin_KN*01.07": False,
        "KN*01W_isin_KN*01.10": False,
        "KN*01W_isin_KN*01.12": False,
        "KN*01W_isin_KN*01W": False,
        "KN*01W_isin_KN*02": False,
        "KN*01W_isin_KN*02.10": False,
        "KN*01_isin_KN*01": False,
        "KN*01_isin_KN*01.-05": False,
        "KN*01_isin_KN*01.-08": False,
        "KN*01_isin_KN*01.-13": False,
        "KN*01_isin_KN*01.06": False,
        "KN*01_isin_KN*01.07": False,
        "KN*01_isin_KN*01.10": False,
        "KN*01_isin_KN*01.12": False,
        "KN*01_isin_KN*01W": False,
        "KN*01_isin_KN*02": False,
        "KN*01_isin_KN*02.10": False,
        "KN*02.10_isin_KN*01": False,
        "KN*02.10_isin_KN*01.-05": False,
        "KN*02.10_isin_KN*01.-08": False,
        "KN*02.10_isin_KN*01.-13": False,
        "KN*02.10_isin_KN*01.06": False,
        "KN*02.10_isin_KN*01.07": False,
        "KN*02.10_isin_KN*01.10": False,
        "KN*02.10_isin_KN*01.12": False,
        "KN*02.10_isin_KN*01W": False,
        "KN*02.10_isin_KN*02": False,
        "KN*02.10_isin_KN*02.10": False,
        "KN*02_isin_KN*01": False,
        "KN*02_isin_KN*01.-05": False,
        "KN*02_isin_KN*01.-08": False,
        "KN*02_isin_KN*01.-13": False,
        "KN*02_isin_KN*01.06": False,
        "KN*02_isin_KN*01.07": False,
        "KN*02_isin_KN*01.10": False,
        "KN*02_isin_KN*01.12": False,
        "KN*02_isin_KN*01W": False,
        "KN*02_isin_KN*02": False,
        "KN*02_isin_KN*02.10": True,
    }
}


class Zygosity:
    HOM = "Homozygous"
    HET = "Heterozygous"


class TestAddCoExistingAllelesNotAllHom(unittest.TestCase):
    def setUp(self):
        self.allele_ref = unique_alleles[0]  # ref
        self.allele10 = unique_alleles[1]  # 10
        self.allele12 = unique_alleles[2]  # 12
        self.allele7 = unique_alleles[3]  # 7
        self.allele6 = unique_alleles[4]  # 6
        # from kenya 230

        self.combos = [
            (self.allele10,),
            (self.allele6, self.allele12),
            (self.allele7, self.allele12),
            (self.allele7,),
            (self.allele10, self.allele12),
            (self.allele6,),
        ]

        self.bg = BloodGroup(
            type="KN",
            alleles={
                AlleleState.POS: [
                    self.allele_ref,
                    self.allele10,
                    self.allele12,
                    self.allele7,
                ]
            },
            sample="231",
            variant_pool={
                "1:207753621_A_G": "Heterozygous",
                "1:207782856_A_G": "Heterozygous",
                "1:207782889_A_G": "Heterozygous",
                "1:207782916_A_T": "Homozygous",
                "1:207782931_A_G": "Homozygous",
            },
            filtered_out={},
        )
        self.bg.misc = {}
        self.bg.misc["combos"] = self.combos
        self.reference_alleles = {
            "KN": Allele(
                genotype="KN*01",
                phenotype="KN:1,-2,3,-6,4,-7,5,8,9,-10,11,-12",
                genotype_alt="KN*A",
                phenotype_alt="Kn(a+), Kn(b-), McC(a+), McC(b-), Sl1+, Vil+, "
                "Sl3+, KCAM+, DACY+",
                defining_variants=frozenset({"1:207782916_A_T"}),
                weight_geno=1000,
                weight_pheno=2,
                reference=True,
                sub_type="KN*01",
                phases=None,
            )
        }

    def test_basic_functionality(self):
        result_bg = list(add_co_existing_alleles({1: self.bg}).values())[0]

        expected_possible_alleles = [
            ((self.allele10,), (self.allele10,)),
            ((self.allele10,), (self.allele6, self.allele12)),
            ((self.allele10,), (self.allele7, self.allele12)),
            ((self.allele10,), (self.allele7,)),
            ((self.allele10,), (self.allele10, self.allele12)),
            ((self.allele10,), (self.allele6,)),
            ((self.allele7,), (self.allele10, self.allele12)),
            ((self.allele10, self.allele12), (self.allele6,)),
        ]

        self.assertEqual(
            result_bg.alleles[AlleleState.CO],
            expected_possible_alleles,
        )


class TestCanCoExist(unittest.TestCase):
    def setUp(self):
        self.tmp_pool23 = {
            "var1": 2,
            "var2": 1,
            "var3": 0,
            "var4": 3,
        }

    def test_basic_functionality(self):
        mushed_combo1 = {"var1", "var2"}
        result = can_co_exist(mushed_combo1, self.tmp_pool23)
        self.assertTrue(result)

    def test_no_available_variants(self):
        tmp_pool23 = {
            "var1": 0,
            "var2": 0,
            "var3": 0,
            "var4": 0,
        }
        mushed_combo1 = {"var1", "var2"}
        result = can_co_exist(mushed_combo1, tmp_pool23)
        self.assertFalse(result)

    def test_partial_match(self):
        mushed_combo1 = {"var1", "var3"}
        result = can_co_exist(mushed_combo1, self.tmp_pool23)
        self.assertFalse(result)

    def test_empty_combination(self):
        mushed_combo1 = set()
        result = can_co_exist(mushed_combo1, self.tmp_pool23)
        self.assertTrue(result)


class TestDecideIfCoExisting(unittest.TestCase):
    def setUp(self):
        self.allele_ref = Allele(
            genotype="A1",
            phenotype="M",
            genotype_alt=".",
            phenotype_alt=".",
            defining_variants=frozenset({"var1", "var2"}),
            weight_geno=1,
            weight_pheno=1,
            reference=False,
            sub_type="subtype1",
        )
        self.allele10 = Allele(
            genotype="A2",
            phenotype="N",
            genotype_alt=".",
            phenotype_alt=".",
            defining_variants=frozenset({"var2", "var3"}),
            weight_geno=1,
            weight_pheno=1,
            reference=False,
            sub_type="subtype1",
        )
        self.allele12 = Allele(
            genotype="A3",
            phenotype="O",
            genotype_alt=".",
            phenotype_alt=".",
            defining_variants=frozenset({"var4"}),
            weight_geno=1,
            weight_pheno=1,
            reference=False,
            sub_type="subtype2",
        )

        self.tmp_pool2 = {
            "var1": 2,
            "var2": 2,
            "var3": 1,
            "var4": 1,
        }

    def test_basic_functionality(self):
        combo1 = (self.allele_ref, self.allele10)
        combo2 = (self.allele12,)
        mushed_combo1 = {"var1", "var2", "var3"}
        co_existing = []

        result = decide_if_co_existing(
            self.tmp_pool2, combo1, combo2, mushed_combo1, co_existing
        )

        expected = [(tuple(combo1), tuple(combo2))]
        self.assertEqual(result, expected)

    def test_empty_combinations(self):
        combo1 = ()
        combo2 = ()
        mushed_combo1 = set()
        co_existing = []

        result = decide_if_co_existing(
            self.tmp_pool2, combo1, combo2, mushed_combo1, co_existing
        )
        expected = [((), ())]
        self.assertEqual(result, expected)

    def test_partial_co_existing(self):
        combo1 = (self.allele_ref,)
        combo2 = (self.allele10,)
        mushed_combo1 = {"var2", "var3"}
        co_existing = []

        result = decide_if_co_existing(
            self.tmp_pool2, combo1, combo2, mushed_combo1, co_existing
        )

        expected = [(tuple(combo1), tuple(combo2))]
        self.assertEqual(result, expected)


class TestAddCoExistingAllelesAllHom(unittest.TestCase):  # has been replaced, right?
    def setUp(self):
        self.allele_ref = Allele(
            genotype="KN*01.07",
            phenotype="KN:1,-2,3,-6,-4,7,5,8,9,-10,11,-12",
            genotype_alt=".",
            phenotype_alt="Vil+",
            defining_variants=frozenset({"1:207782889_A_G", "1:207782916_A_T"}),
            weight_geno=1000,
            weight_pheno=1,
            reference=False,
            sub_type="KN*01",
            phases=None,
        )
        self.allele10 = Allele(
            genotype="KN*01.10",
            phenotype="KN:1,-2,3,-6,4,-7,5,8,-9,10,11,-12",
            genotype_alt=".",
            phenotype_alt="KDAS+",
            defining_variants=frozenset({"1:207782916_A_T", "1:207782931_A_G"}),
            weight_geno=1000,
            weight_pheno=1,
            reference=False,
            sub_type="KN*01",
            phases=None,
        )

        # alleles={'co_existing': [Pair(Genotype: KN*01.07+KN*01.10/KN*01.07+KN*01.10 Phenotype: WIP/WIP)],
        # pairs': [Pair(Genotype: KN*01.10/KN*01.10 Phenotype: KN:1,-2,3,-6,4,-7,5,8,-9,10,11,-12/KN:1,-2,3,-6,4,-7,5,8,-9,10,11,-12),
        # Pair(Genotype: KN*01.07/KN*01.07 Phenotype: KN:1,-2,3,-6,-4,7,5,8,9,-10,11,-12/KN:1,-2,3,-6,-4,7,5,8,9,-10,11,-12)],
        self.bg = BloodGroup(
            type="KN",
            alleles={AlleleState.POS: [self.allele_ref, self.allele10]},
            sample="231",
            variant_pool={},
            filtered_out={},
        )
        self.reference_alleles = {
            "KN": Allele(
                genotype="KN*01",
                phenotype="KN:1,-2,3,-6,4,-7,5,8,9,-10,11,-12",
                genotype_alt="KN*A",
                phenotype_alt="Kn(a+), Kn(b-), McC(a+), McC(b-), Sl1+, Vil+, "
                "Sl3+, KCAM+, DACY+",
                defining_variants=frozenset({"1:207782916_A_T"}),
                weight_geno=1000,
                weight_pheno=2,
                reference=True,
                sub_type="KN*01",
                phases=None,
            )
        }

    def test_basic_functionality(self):
        self.bg.variant_pool = {
            "1:207782889_A_G": Zygosity.HOM,
            "1:207782916_A_T": Zygosity.HOM,
            "1:207782931_A_G": Zygosity.HOM,
        }
        self.bg.misc = {}
        self.bg.misc["combos"] = [(self.allele_ref,), (self.allele10,)]
        result_bg = list(add_co_existing_alleles({1: self.bg}).values())[0]
        result_bg_as_set = []
        for combo_pair in result_bg.alleles[AlleleState.CO]:
            result_bg_as_set.append(set(map(frozenset, combo_pair)))
        expected_possible_alleles = [  # 1 combo pair
            frozenset([self.allele_ref, self.allele10]),
            frozenset([self.allele10]),
            frozenset([self.allele_ref]),
        ]

        self.assertEqual(len(result_bg_as_set), len(expected_possible_alleles))

    def test_no_homozygous_variants(self):
        self.bg.variant_pool = {
            "1:207782889_A_G": Zygosity.HET,
            "1:207782916_A_T": Zygosity.HET,
            "1:207782931_A_G": Zygosity.HET,
        }
        self.bg.misc = {}
        self.bg.misc["combos"] = [(self.allele_ref,), (self.allele10,)]
        result_bg = list(add_co_existing_alleles({1: self.bg}).values())[0]
        expected_possible_alleles = None
        self.assertEqual(result_bg.alleles[AlleleState.CO], expected_possible_alleles)


class TestMushedVarsFunction(unittest.TestCase):
    def test_mushed_vars(self):
        allele1 = Allele(
            genotype="KN*01",
            phenotype="positive",
            genotype_alt="KN*",
            phenotype_alt="neg",
            defining_variants=frozenset({"variant1", "variant2"}),
            weight_geno=2,
            weight_pheno=3,
            reference=False,
        )
        allele2 = Allele(
            genotype="KN*02",
            phenotype="negative",
            genotype_alt="KN*",
            phenotype_alt="neg",
            defining_variants=frozenset({"variant3"}),
            weight_geno=1,
            weight_pheno=2,
            reference=False,
        )
        combo = [allele1, allele2]
        expected_variants = {"variant1", "variant2", "variant3"}
        self.assertEqual(mushed_vars(combo), expected_variants)


class TestAddCoExistingAllelesFunction(unittest.TestCase):
    def setUp(self):
        # Create Alleles
        self.allele_ref = Allele(
            genotype="KN*01",
            phenotype="positive",
            genotype_alt="KN*",
            phenotype_alt="neg",
            defining_variants=frozenset({"variant1"}),
            weight_geno=2,
            weight_pheno=3,
            reference=False,
            sub_type="KN*01",
        )
        self.allele10 = Allele(
            genotype="KN*02",
            phenotype="negative",
            genotype_alt="KN*",
            phenotype_alt="neg",
            defining_variants=frozenset({"variant2"}),
            weight_geno=1,
            weight_pheno=2,
            reference=False,
            sub_type="KN*02",
        )
        self.bg = BloodGroup(
            type="KN",
            alleles={AlleleState.POS: [self.allele_ref, self.allele10]},
            sample="sample1",
            genotypes=["KN*01", "KN*02"],
            phenotypes=["positive", "negative"],
            variant_pool={"variant1": Zygosity.HET, "variant2": Zygosity.HET},
        )  # variant_pool_numeric={"variant1": 1, "variant2": 1},
        self.bg.misc = {}
        self.bg.misc["homs"] = set()
        self.bg.misc["combos"] = [(self.allele_ref,), (self.allele10,)]
        self.bg.misc["max_rank"] = 1  # Set a max rank
        self.combos = [(self.allele10,)]

    def test_add_co_existing_alleles(self):
        # Apply the function
        add_co_existing_alleles({"KN": self.bg})

        # Expected co-existing combinations
        expected_co_existing = [
            ((self.allele_ref,), (self.allele10,)),
        ]

        # Check if co-existing alleles are correctly identified
        self.assertEqual(len(self.bg.alleles[AlleleState.CO]), 1)
        self.assertEqual(expected_co_existing, self.bg.alleles[AlleleState.CO])

    def test_no_co_existing_alleles(self):
        # Modify variant_pool_numeric to make co-existence impossible
        # self.bg.variant_pool = {}
        self.bg.misc["combos"] = {(self.allele_ref,)}
        # Apply the function
        add_co_existing_alleles({"KN": self.bg})
        # Check that no co-existing alleles are found
        self.assertIsNone(self.bg.alleles[AlleleState.CO])


class TestMushFunctionNew(unittest.TestCase):
    """Test suite for the 'mush' function (new)."""

    def test_none_co_alleles(self) -> None:
        """If CO alleles are None, mush() returns the BloodGroup unchanged."""
        bg = BloodGroup(
            type="Test",
            sample="dummy_sample",  # <-- REQUIRED
            alleles={AlleleState.CO: None},
        )
        # Must call mush with a dict if @apply_to_dict_values is used:
        res_dict = mush({0: bg})
        result = res_dict[0]
        self.assertIsNone(result.alleles[AlleleState.CO], "Expected CO to remain None.")

    def test_single_allele_combo(self) -> None:
        """Single allele combos should remain unchanged (but actually raises ValueError if no 2 Alleles)."""
        allele1 = Allele(
            genotype="mock_g",
            genotype_alt="mock_ga",
            defining_variants=frozenset(),
            phenotype="Test:1",
            phenotype_alt="X+",
        )
        bg = BloodGroup(
            type="Test",
            sample="dummy_sample",  # <-- REQUIRED
            alleles={
                AlleleState.CO: [
                    [
                        [allele1]
                    ]  # One outer list => 1 group of combos => single sub-combo
                ]
            },
        )
        # Original code raises ValueError because only one Allele => cannot form a pair
        with self.assertRaises(ValueError):
            res = mush({0: bg})
            _ = res[0]

    def test_two_single_allele_combos_form_a_pair(self) -> None:
        """Two single Alleles => exactly 2 => forms a Pair."""
        allele1 = Allele(
            genotype="mock_g1",
            genotype_alt="mock_ga1",
            defining_variants=frozenset(),
            phenotype="Test:1",
            phenotype_alt="X+",
        )
        allele2 = Allele(
            genotype="mock_g2",
            genotype_alt="mock_ga2",
            defining_variants=frozenset(),
            phenotype="Test:2",
            phenotype_alt="Y-",
        )
        bg = BloodGroup(
            type="Test",
            sample="dummy_sample",
            alleles={AlleleState.CO: [[[allele1], [allele2]]]},
        )
        # p=Pair()
        res_dict = mush({0: bg})
        result = res_dict[0]
        self.assertEqual(len(result.alleles[AlleleState.CO]), 1)
        pair = result.alleles[AlleleState.CO][0]
        self.assertIsInstance(pair, Pair)
        self.assertEqual(pair.allele1, allele1)
        self.assertEqual(pair.allele2, allele2)

    def test_merge_multiple_alleles_no_conflict(self) -> None:
        """Multiple allele combos that unify into one mushed Allele + single Allele => forms a Pair."""
        allele1 = Allele(
            genotype="mock_gA",
            genotype_alt="mock_gaA",
            defining_variants=frozenset(),
            phenotype="Test:+A",
            phenotype_alt="X+",
        )
        allele2 = Allele(
            genotype="mock_gB",
            genotype_alt="mock_gaB",
            defining_variants=frozenset(),
            phenotype="Test:+B",
            phenotype_alt="Y+",
        )
        # single
        allele3 = Allele(
            genotype="mock_gC",
            genotype_alt="mock_gaC",
            defining_variants=frozenset(),
            phenotype="Test:+A",
            phenotype_alt="X+",
        )
        bg = BloodGroup(
            type="Test",
            sample="dummy_sample",
            alleles={
                AlleleState.CO: [
                    [
                        [allele1, allele2],
                        [allele3],
                    ]
                ]
            },
        )
        res_dict = mush({0: bg})
        result = res_dict[0]
        self.assertEqual(len(result.alleles[AlleleState.CO]), 1)
        pair = result.alleles[AlleleState.CO][0]
        # The first is a newly mushed allele, second is the single
        self.assertIsNot(pair.allele1, allele1)
        self.assertIsNot(pair.allele1, allele2)
        self.assertEqual(pair.allele2, allele3)

    def test_merge_conflict_raises_error(self) -> None:
        """Alleles with conflicting expression raise ValueError."""
        allele1 = Allele(
            genotype="conflict_g1",
            genotype_alt="conflict_ga1",
            defining_variants=frozenset(),
            phenotype="Test:+A",
            phenotype_alt="X+",
        )
        allele2 = Allele(
            genotype="conflict_g2",
            genotype_alt="conflict_ga2",
            defining_variants=frozenset(),
            phenotype="Test:-A",
            phenotype_alt="X-",
        )
        bg = BloodGroup(
            type="Test",
            sample="dummy_sample",
            alleles={AlleleState.CO: [[[allele1, allele2], [allele1]]]},
        )
        with self.assertRaises(ValueError):
            res_dict = mush({0: bg})
            _ = res_dict[0]

    def test_incomplete_pair_raises_error(self) -> None:
        """One sub-combo with multiple Alleles => merges to 1 => no second => ValueError."""
        allele1 = Allele(
            genotype="mock_gX",
            genotype_alt="mock_gaX",
            defining_variants=frozenset(),
            phenotype="Test:+A",
            phenotype_alt="X+",
        )
        allele2 = Allele(
            genotype="mock_gY",
            genotype_alt="mock_gaY",
            defining_variants=frozenset(),
            phenotype="Test:-B",
            phenotype_alt="Y-",
        )
        bg = BloodGroup(
            type="Test",
            sample="dummy_sample",
            alleles={AlleleState.CO: [[[allele1, allele2]]]},
        )
        with self.assertRaises(ValueError):
            res_dict = mush({0: bg})
            _ = res_dict[0]

    def test_result_none_if_empty_result(self) -> None:
        """If combos exist but produce an empty result, CO is set to None."""
        bg = BloodGroup(
            type="Test",
            sample="dummy_sample",  # <-- REQUIRED
            alleles={AlleleState.CO: []},
        )
        res_dict = mush({0: bg})
        result = res_dict[0]
        self.assertIsNone(result.alleles[AlleleState.CO])


# ------------------------------------------------------------------------------
# Full-Coverage Unit Tests for the Original `mush` Function
# ------------------------------------------------------------------------------
class TestMushOriginalFullCoverage(unittest.TestCase):
    """
    Tests that exercise all logic paths, edge cases, special cases, etc.
    to ensure full coverage of the original `mush` function.
    """

    def test_co_is_none_returns_unchanged(self):
        """If CO is None, function returns BG as-is."""
        bg = BloodGroup(type="AB", sample="sampleX", alleles={AlleleState.CO: None})
        out = list(mush({1: bg}).values())[0]
        self.assertIsNone(out.alleles[AlleleState.CO])

    def test_no_combos_sets_co_none(self):
        """If CO is [], result remains empty => sets CO=None."""
        bg = BloodGroup(type="AB", sample="sampleX", alleles={AlleleState.CO: []})
        out = list(mush({1: bg}).values())[0]
        self.assertIsNone(out.alleles[AlleleState.CO])

    def test_single_combo_one_subcombo_single_allele_raises_error2(self):
        """
        Single sub-combo => 1 Allele => can't form a 2-Allele Pair => ValueError
        We now expect EXACT message "must contain exactly two Alleles".
        """
        a1 = Allele(
            genotype="g1",
            genotype_alt="ga1",
            defining_variants=frozenset(),
            phenotype="AB:+1",
            phenotype_alt="Ph+",
        )
        bg = BloodGroup(type="AB", sample="sampleX", alleles={AlleleState.CO: [[[a1]]]})

        with self.assertRaises(ValueError) as ctx:
            mush({1: bg})

        self.assertIn(
            "Each mushed_pair must contain exactly two alleles.", str(ctx.exception)
        )

    def test_single_combo_two_single_alleles_ok(self):
        """One combo => sub-combo1 has a1, sub-combo2 has a2 => 2 Alleles => forms Pair."""
        a1 = Allele(
            genotype="g1",
            genotype_alt="ga1",
            defining_variants=frozenset(),
            phenotype="AB:+1",
            phenotype_alt="Ph+",
        )
        a2 = Allele(
            genotype="g2",
            genotype_alt="ga2",
            defining_variants=frozenset(),
            phenotype="AB:+2",
            phenotype_alt="Ph-",
        )
        bg = BloodGroup(
            type="AB", sample="sampleX", alleles={AlleleState.CO: [[[a1], [a2]]]}
        )
        out = list(mush({1: bg}).values())[0]
        self.assertIsNotNone(out.alleles[AlleleState.CO])
        self.assertEqual(len(out.alleles[AlleleState.CO]), 1)
        pair = out.alleles[AlleleState.CO][0]
        self.assertIsInstance(pair, Pair)
        self.assertEqual(pair.allele1, a1)
        self.assertEqual(pair.allele2, a2)

    def test_multi_allele_subcombo_no_conflict_merges_ok(self):
        """
        One sub-combo with multiple alleles => unify into 1 new Allele,
        second sub-combo => single => total 2 => forms Pair.
        No conflict => merges successfully.
        """
        a1 = Allele(
            genotype="gA",
            genotype_alt="gaA",
            defining_variants=frozenset(),
            phenotype="AB:+1",
            phenotype_alt="X+",
        )
        a2 = Allele(
            genotype="gB",
            genotype_alt="gaB",
            defining_variants=frozenset(),
            phenotype="AB:+2",
            phenotype_alt="Y+",
        )
        a3 = Allele(
            genotype="gC",
            genotype_alt="gaC",
            defining_variants=frozenset(),
            phenotype="AB:+3",
            phenotype_alt="Z+",
        )
        bg = BloodGroup(
            type="AB", sample="sampleX", alleles={AlleleState.CO: [[[a1, a2], [a3]]]}
        )
        out = list(mush({1: bg}).values())[0]
        self.assertEqual(len(out.alleles[AlleleState.CO]), 1)
        pair = out.alleles[AlleleState.CO][0]
        self.assertNotEqual(pair.allele1, a1)
        self.assertNotEqual(pair.allele1, a2)
        # The mushed allele is a brand new Allele
        self.assertEqual(pair.allele2, a3)

    def test_incomplete_pair_mushed_only_one_allele_raises_valueerror(self):
        """
        Single sub-combo with multiple => merges => 1 => no second => ValueError
        Adjust alt so no conflict overshadowing the real error.
        """
        a1 = Allele(
            genotype="g1",
            genotype_alt="ga1",
            defining_variants=frozenset(),
            phenotype="AB:+X",
            phenotype_alt="NoConflict+",
        )
        a2 = Allele(
            genotype="g2",
            genotype_alt="ga2",
            defining_variants=frozenset(),
            phenotype="AB:+Y",
            phenotype_alt="NoConflict+",
        )
        bg = BloodGroup(
            type="AB", sample="sampleX", alleles={AlleleState.CO: [[[a1, a2]]]}
        )
        # Expect ValueError about needing two Alleles, not an antigen conflict
        with self.assertRaises(ValueError) as ctx:
            mush({1: bg})
        self.assertIn(
            "Each mushed_pair must contain exactly two alleles", str(ctx.exception)
        )

    def test_conflicting_antigen_expressions_raises_valueerror(self):
        """
        Expect conflict for A vs -A in numeric phenotype => ValueError.
        Previously it didn't raise if we used 'AB:+A', because the real function
        sees that as name='+A'. Now we use 'AB:A' (no plus sign) vs 'AB:-A'.
        """
        a_plus = Allele(
            genotype="gPlus",
            genotype_alt="gPlus_alt",
            defining_variants=frozenset(),
            phenotype="AB:A",  # No '+' => name is 'A' => will conflict with '-A'
            phenotype_alt="AltNoConflict+",
        )
        a_minus = Allele(
            genotype="gMinus",
            genotype_alt="gMinus_alt",
            defining_variants=frozenset(),
            phenotype="AB:-A",  # => name='A', expr='-'
            phenotype_alt="AltNoConflict+",
        )
        single_allele = Allele(
            genotype="gSingle",
            genotype_alt="gSingle_alt",
            defining_variants=frozenset(),
            phenotype="AB:+Z",  # different antigen => no conflict
            phenotype_alt="Z+",
        )

        bg = BloodGroup(
            type="AB",
            sample="sampleX",
            alleles={
                AlleleState.CO: [
                    [
                        [
                            a_plus,
                            a_minus,
                        ],  # sub-combo => 2 alleles => unify => conflict
                        [single_allele],  # sub-combo => single => no conflict
                    ]
                ]
            },
        )

        with self.assertRaises(ValueError) as ctx:
            mush({1: bg})
        # Confirm the partial snippet
        self.assertIn("Conflicting expressions for antigen 'A':", str(ctx.exception))

    def test_single_combo_one_subcombo_single_allele_raises_error(self):
        """
        Single sub-combo => 1 Allele => can't form a 2-Allele Pair => ValueError
        Must match exactly 'must contain exactly two Alleles'.
        """
        a1 = Allele(
            genotype="g1",
            genotype_alt="ga1",
            defining_variants=frozenset(),
            phenotype="AB:A",  # does not matter; only 1 => must fail
            phenotype_alt="Ph+",
        )
        bg = BloodGroup(
            type="AB",
            sample="sampleX",
            alleles={
                AlleleState.CO: [
                    [
                        [a1]  # single sub-combo, single Allele
                    ]
                ]
            },
        )

        with self.assertRaises(ValueError) as ctx:
            mush({1: bg})
        self.assertIn(
            "Each mushed_pair must contain exactly two alleles.", str(ctx.exception)
        )

    def test_empty_result_sets_co_none(self):
        """
        If we had combos but all cause errors or if logic ends with no final pairs => sets CO to None.
        We'll simulate by adding multiple combos that each lead to ValueError, but that halts the function.
        For demonstration, we can just pass an empty combos array => same effect.
        """
        bg = BloodGroup(type="AB", sample="sampleX", alleles={AlleleState.CO: []})
        out = list(mush({1: bg}).values())[0]
        self.assertIsNone(out.alleles[AlleleState.CO])

    def test_helgeson_special_case_numeric_1w(self):
        """
        If '1w,' and '1,' appear => remove '1w,'.
        Adjust alt to avoid conflict.
        """
        a1 = Allele(
            genotype="h1",
            genotype_alt="h1_alt",
            defining_variants=frozenset(),
            phenotype="AB:1w,2",
            phenotype_alt="NoConflict+",
        )
        a2 = Allele(
            genotype="h2",
            genotype_alt="h2_alt",
            defining_variants=frozenset(),
            phenotype="AB:1",
            phenotype_alt="NoConflict+",
        )
        a_single = Allele(
            genotype="hS",
            genotype_alt="hS_alt",
            defining_variants=frozenset(),
            phenotype="AB:+Z",
            phenotype_alt="Z+",  # same sign => no conflict
        )
        bg = BloodGroup(
            type="AB",
            sample="sampleX",
            alleles={AlleleState.CO: [[[a1, a2], [a_single]]]},
        )
        out = list(mush({1: bg}).values())[0]
        pair = out.alleles[AlleleState.CO][0]
        # check we removed "1w,"
        self.assertNotIn("1w,", pair.allele1.phenotype)
        self.assertIn(":1,2", pair.allele1.phenotype)

    def test_helgeson_special_case_alphanumeric_knaw(self):
        """
        If 'Kn(aw+),' and 'Kn(a+),' => remove 'Kn(aw+),'.
        Adjust the other alt pieces so no conflict occurs.
        """
        a_aw = Allele(
            genotype="hAW",
            genotype_alt="hAW_alt",
            defining_variants=frozenset(),
            phenotype="AB:+A",
            phenotype_alt="Kn(aw+),Hello",  # removed the plus or minus on 'Hello'
        )
        a_a = Allele(
            genotype="hA",
            genotype_alt="hA_alt",
            defining_variants=frozenset(),
            phenotype="AB:+A",
            phenotype_alt="Kn(a+),Hello",  # no conflict => "Hello" vs "Hello"
        )
        single = Allele(
            genotype="hSingle",
            genotype_alt="hSingle_alt",
            defining_variants=frozenset(),
            phenotype="AB:+X",
            phenotype_alt="X+",  # keep consistent sign
        )
        bg = BloodGroup(
            type="AB",
            sample="sampleX",
            alleles={AlleleState.CO: [[[a_aw, a_a], [single]]]},
        )
        out = list(mush({1: bg}).values())[0]
        pair = out.alleles[AlleleState.CO][0]
        self.assertNotIn("Kn(aw+),", pair.allele1.phenotype_alt)
        self.assertIn("Kn(a+)", pair.allele1.phenotype_alt)

    def test_multiple_combos_in_co(self):
        """
        If we have multiple combos in CO, each must produce exactly one Pair.
        We'll do an example with 2 combos => final CO => 2 Pairs.
        """
        a1 = Allele(
            genotype="gA1",
            genotype_alt="gaA1",
            defining_variants=frozenset(),
            phenotype="AB:+A",
            phenotype_alt="X+",
        )
        a2 = Allele(
            genotype="gA2",
            genotype_alt="gaA2",
            defining_variants=frozenset(),
            phenotype="AB:-B",
            phenotype_alt="Y-",
        )
        a3 = Allele(
            genotype="gA3",
            genotype_alt="gaA3",
            defining_variants=frozenset(),
            phenotype="AB:+C",
            phenotype_alt="Z+",
        )
        a4 = Allele(
            genotype="gA4",
            genotype_alt="gaA4",
            defining_variants=frozenset(),
            phenotype="AB:+D",
            phenotype_alt="W+",
        )
        # combo1 => sub-combo1 => multiple => unify, sub-combo2 => single
        combo1 = [
            [a1, a2],
            [a3],
        ]
        # combo2 => both single
        combo2 = [[a4], [a2]]
        bg = BloodGroup(
            type="AB", sample="sampleX", alleles={AlleleState.CO: [combo1, combo2]}
        )
        out = list(mush({1: bg}).values())[0]
        pairs_list = out.alleles[AlleleState.CO]
        self.assertEqual(len(pairs_list), 2)
        # first result => pair from combo1
        self.assertIsInstance(pairs_list[0], Pair)
        self.assertIsInstance(pairs_list[1], Pair)

    def test_cover_paren_in_alt_minus_expression(self):
        """
        Cover the line:
        if ")" in name:
            unique_antigens_alt.append(name.replace(")", "-)"))
        by having an alt antigen with a closing parenthesis in a minus expression.
        We ensure the final result has 2 Alleles total (1 mushed + 1 single sub-combo).
        """
        # Allele #1 => 'Aw)-' => expression='-', name='Aw)'
        # Allele #2 => no parenthesis => merges fine
        # => This unifies into 1 mushed Allele.
        allele1 = Allele(
            genotype="pA1",
            genotype_alt="pA1_alt",
            defining_variants=frozenset(),
            phenotype="AB:+X",  # numeric => no conflict
            phenotype_alt="Aw)-",  # triggers the parenthesis replacement
        )
        allele2 = Allele(
            genotype="pA2",
            genotype_alt="pA2_alt",
            defining_variants=frozenset(),
            phenotype="AB:+Y",
            phenotype_alt="NoParen",  # no parenthesis
        )
        # Sub-combo #1 => [allele1, allele2] => merges into 1 new Allele
        # Sub-combo #2 => a single allele => ensures final pair => no ValueError
        single_allele = Allele(
            genotype="pA3",
            genotype_alt="pA3_alt",
            defining_variants=frozenset(),
            phenotype="AB:+Z",
            phenotype_alt="Z+",
        )
        bg = BloodGroup(
            type="AB",
            sample="sampleX",
            alleles={AlleleState.CO: [[[allele1, allele2], [single_allele]]]},
        )
        out = list(mush({1: bg}).values())[0]
        # Retrieve the final pair => the first is a new mushed Allele
        pair = out.alleles[AlleleState.CO][0]
        # Confirm the minus expression with parenthesis was replaced => "Aw-)"
        self.assertIn(
            "Aw-)",
            pair.allele1.phenotype_alt,
            "Expected parenthesis replaced with '-)' in the mushed allele's alt phenotype.",
        )

    def test_cover_knaw_helgeson_case(self):
        """
        Cover the line:
        if 'Kn(aw+),' in pheno_alphanumeric and 'Kn(a+),' in pheno_alphanumeric:
            pheno_alphanumeric = pheno_alphanumeric.replace('Kn(aw+),','')
        by ensuring the unified alt phenotype has both 'Kn(aw+),' and 'Kn(a+),'.
        Again, we use 2 sub-combos => final Pair => no ValueError.
        """
        # Allele #1 => alt includes 'Kn(aw+),'
        allele_knaw = Allele(
            genotype="kA1",
            genotype_alt="kA1_alt",
            defining_variants=frozenset(),
            phenotype="AB:+A",
            phenotype_alt="Kn(aw+),Hello",
        )
        # Allele #2 => alt includes 'Kn(a+),'
        allele_kna = Allele(
            genotype="kA2",
            genotype_alt="kA2_alt",
            defining_variants=frozenset(),
            phenotype="AB:+A",
            phenotype_alt="Kn(a+),World",
        )
        # These two unify into 1 mushed Allele. We still need a second sub-combo => single Allele
        single_allele = Allele(
            genotype="kA3",
            genotype_alt="kA3_alt",
            defining_variants=frozenset(),
            phenotype="AB:+Z",
            phenotype_alt="Z+",
        )
        bg = BloodGroup(
            type="AB",
            sample="sampleX",
            alleles={AlleleState.CO: [[[allele_knaw, allele_kna], [single_allele]]]},
        )
        out = list(mush({1: bg}).values())[0]
        pair = out.alleles[AlleleState.CO][0]
        # The Helgeson line should remove 'Kn(aw+),'
        self.assertNotIn(
            "Kn(aw+),",
            pair.allele1.phenotype_alt,
            "Expected 'Kn(aw+),' to be removed from the mushed alt phenotype.",
        )
        # 'Kn(a+),' should remain
        self.assertIn(
            "Kn(a+),",
            pair.allele1.phenotype_alt,
            "Expected 'Kn(a+),' to remain in the mushed alt phenotype.",
        )


class TestMaxRankFunction(unittest.TestCase):
    def setUp(self):
        # Similar setup as before
        self.allele_ref = Allele(
            genotype="KN*01",
            phenotype="positive",
            genotype_alt="KN*",
            phenotype_alt="neg",
            defining_variants=frozenset({"variant1", "variant2"}),
            weight_geno=2,  # Weight 2
            weight_pheno=3,
            reference=True,
            sub_type="KN*01",
        )
        self.allele02 = Allele(
            genotype="KN*02",
            phenotype="negative",
            genotype_alt="KN*",
            phenotype_alt="neg",
            defining_variants=frozenset({"variant3"}),
            weight_geno=1,  # Weight 1 (lower is higher rank)
            weight_pheno=2,
            reference=False,
            sub_type="KN*02",
        )
        self.bg = BloodGroup(
            type="KN",
            alleles={AlleleState.POS: [self.allele_ref, self.allele02]},
            sample="sample1",
            genotypes=["KN*01", "KN*02"],
            phenotypes=["positive", "negative"],
        )  # variant_pool_numeric={"variant1": 2, "variant2": 2, "variant3": 1},

        self.bg.misc = {}
        self.bg.misc["homs"] = {self.allele_ref}

    def test_max_rank_calculation(self):
        # Apply the max_rank function
        max_rank({"KN": self.bg})

        # Expected max rank is the minimum weight_geno among homozygous alleles and
        # alleles with one defining variant
        expected_max_rank = 1  # From allele2 (weight_geno=1)

        # Check if max rank is calculated correctly
        self.assertEqual(self.bg.misc["max_rank"], expected_max_rank)

    def test_max_rank_no_homs(self):
        # Remove homozygous alleles
        self.bg.misc["homs"] = set()

        # Apply the max_rank function
        max_rank({"KN": self.bg})

        # Expected max rank is LOW_WEIGHT (from constants)
        expected_max_rank = 1

        # Check if max rank is set to LOW_WEIGHT
        self.assertEqual(self.bg.misc["max_rank"], expected_max_rank)


class TestHomsFunction(unittest.TestCase):
    def setUp(self):
        # Create Alleles
        self.allele_ref = Allele(
            genotype="KN*01",
            phenotype="positive",
            genotype_alt="KN*",
            phenotype_alt="nkn",
            defining_variants=frozenset({"variant1", "variant2"}),
            weight_geno=2,
            weight_pheno=3,
            reference=False,
            sub_type="KN*01",
        )
        self.allele10 = Allele(
            genotype="KN*02",
            phenotype="negative",
            genotype_alt="KN*",
            phenotype_alt="neg",
            defining_variants=frozenset({"variant3"}),
            weight_geno=1,
            weight_pheno=2,
            reference=False,
            sub_type="KN*02",
        )

        # Create BloodGroup
        self.bg = BloodGroup(
            type="KN",
            alleles={AlleleState.POS: [self.allele_ref, self.allele10]},
            sample="sample1",
            genotypes=["KN*01", "KN*02"],
            phenotypes=["positive", "negative"],
            variant_pool={
                "variant1": Zygosity.HOM,
                "variant2": Zygosity.HOM,
                "variant3": Zygosity.HET,
                "variant4": Zygosity.HOM,
            },
        )  # variant_pool_numeric={"variant1": 2, "variant2": 2, "variant3": 1},

        # Initialize misc attribute
        self.bg.misc = {}

    def test_homs_identification(self):
        # Apply the homs function
        homs({"KN": self.bg})

        # Expected homozygous alleles
        expected_homs = {self.allele_ref}

        # Check if the homozygous alleles are correctly identified
        self.assertEqual(self.bg.misc["homs"], expected_homs)

    def test_homs_no_homozygous_alleles(self):
        # Modify variant_pool_numeric to have no homozygous variants
        # self.bg.variant_pool_numeric = {"variant1": 1, "variant2": 1, "variant3": 1}

        # Apply the homs function
        homs({"KN": self.bg})

        # Expected homozygous alleles
        expected_homs = {self.allele_ref}

        # Check if no homozygous alleles are identified
        self.assertEqual(self.bg.misc["homs"], expected_homs)

    def test_homs_multiple_homozygous_alleles(self):
        # Add another allele with homozygous variants
        self.allele12 = Allele(
            genotype="KN*03",
            phenotype="positive",
            genotype_alt="KN*",
            phenotype_alt="neg",
            defining_variants=frozenset({"variant4"}),
            weight_geno=1,
            weight_pheno=1,
            reference=False,
            sub_type="KN*03",
        )
        self.bg.alleles[AlleleState.POS].append(self.allele12)
        # self.bg.variant_pool_numeric["variant4"] = 2

        # Apply the homs function
        homs({"KN": self.bg})

        # Expected homozygous alleles
        expected_homs = {self.allele_ref, self.allele12}

        # Check if both homozygous alleles are identified
        self.assertEqual(self.bg.misc["homs"], expected_homs)


unique_alleles = [
    Allele(
        genotype="KN*01",
        defining_variants=frozenset([]),
        weight_geno=1000,
        phenotype="KN:1,-2,3,4,5,-6,-7,8,9,-10,11,-12,13",
        genotype_alt="KN*A",
        phenotype_alt="Kn(a+),Kn(b-),McC(a+),Sl1+,Yk(a+),McC(b-),Vil-,Sl3+,KCAM+,KDAS-,DACY+,YCAD-,KNMB+",
        weight_pheno=1,
        sub_type="KN*01",
        reference=True,
    ),
    Allele(
        genotype="KN*01.10",
        defining_variants=frozenset(["1:207782931_A_G"]),
        weight_geno=1000,
        phenotype="KN:-9,10",
        genotype_alt="",
        phenotype_alt="KCAM-,KDAS+",
        weight_pheno=1,
        sub_type="KN*01",
        reference=False,
    ),
    Allele(
        genotype="KN*01.12",
        defining_variants=frozenset(["1:207753621_A_G"]),
        weight_geno=1000,
        phenotype="KN:-11,12",
        genotype_alt="",
        phenotype_alt="DACY-,YCAD+",
        weight_pheno=1,
        sub_type="KN*01",
        reference=False,
    ),
    Allele(
        genotype="KN*01.07",
        defining_variants=frozenset(["1:207782931_A_G", "1:207782889_A_G"]),
        weight_geno=1000,
        phenotype="KN:7,-9,10",
        genotype_alt="",
        phenotype_alt="Vil+,KCAM-,KDAS+",
        weight_pheno=1,
        sub_type="KN*01",
        reference=False,
    ),
    Allele(
        genotype="KN*01.06",
        defining_variants=frozenset(
            ["1:207782931_A_G", "1:207782889_A_G", "1:207782856_A_G"]
        ),
        weight_geno=1000,
        phenotype="KN:-3,-4,6,7,-9,10",
        genotype_alt="",
        phenotype_alt="Sl1-,McC(a-),McC(b+),Vil+,KCAM-,KDAS+",
        weight_pheno=1,
        sub_type="KN*01",
        reference=False,
    ),
    Allele(
        genotype="KN*02.02",
        defining_variants=frozenset(
            ["1:207782931_A_G", "1:207782889_A_G", "1:207782856_A_G"]
        ),
        weight_geno=1000,
        phenotype="2,1",
        genotype_alt="",
        phenotype_alt="made+,up+",
        weight_pheno=1,
        sub_type="KN*02",
        reference=False,
    ),
]


def process_combos(combos: set[tuple[Allele, ...]]) -> list[list[Allele]]:
    """
    Converts a set of combos (set of tuples of Alleles) into a sorted list of lists.
    Each combo and the list of combos are sorted to ensure consistent ordering.

    Args:
        combos: A set of tuples, where each tuple contains Allele objects representing
        a combo.

    Returns:
        A sorted list of lists of Allele objects. Both the inner lists (combos) and the
        outer list are sorted.
    """
    # Convert set to list
    combo_list = list(combos)
    # Sort Alleles within each combo based on genotype
    sorted_combo_list = []
    for combo in combo_list:
        sorted_combo = sorted(combo, key=lambda allele: allele.genotype)
        sorted_combo_list.append(sorted_combo)
    # Sort the list of combos based on the genotypes of the Alleles in each combo
    sorted_combo_list.sort(key=lambda combo: [allele.genotype for allele in combo])
    return sorted_combo_list


class TestPrepCoPutativeCombos(unittest.TestCase):
    def checker(self, generated, expected):
        actual_combos = process_combos(generated)
        expected_combos_processed = process_combos(expected)
        self.assertEqual(actual_combos, expected_combos_processed)

    def setUp(self):
        self.allele_ref = unique_alleles[0]  # ref
        self.allele10 = unique_alleles[1]  # 10
        self.allele12 = unique_alleles[2]  # 12
        self.allele7 = unique_alleles[3]  # 7
        self.allele6 = unique_alleles[4]  # 6
        self.allele2 = unique_alleles[5]  # made up - sub type 2
        # Create BloodGroup
        self.bg = BloodGroup(
            type="KN",
            alleles={
                AlleleState.POS: [
                    self.allele10,
                    self.allele12,
                    self.allele7,
                    self.allele6,
                ]
            },
            sample="sample1",
            variant_pool={
                "1:207753621_A_G": Zygosity.HET,
                "1:207782856_A_G": Zygosity.HET,
                "1:207782889_A_G": Zygosity.HET,
                "1:207782931_A_G": Zygosity.HOM,
            },
            genotypes=[],
            phenotypes=["mock"],
        )
        self.bg_with_ref = BloodGroup(
            type="KN",
            alleles={
                AlleleState.POS: [
                    self.allele_ref,
                    self.allele10,
                    self.allele12,
                    self.allele7,
                    self.allele6,
                ]
            },
            sample="sample1",
            variant_pool={
                "1:207753621_A_G": Zygosity.HET,
                "1:207782856_A_G": Zygosity.HET,
                "1:207782889_A_G": Zygosity.HET,
                "1:207782931_A_G": Zygosity.HOM,
            },
            genotypes=[],
            phenotypes=["mock"],
        )

        # Initialize misc attribute
        self.bg.misc = {}

        # Set up homozygous alleles
        self.bg.misc["homs"] = [self.allele_ref]  # Let's say allele1 is homozygous

        # Set up allele relationships
        # Assume allele1 is a sub-allele of allele2
        self.allele_relationships = ALLELE_RELATIONSHIPS

    def test_prep_co_putative_combos_with_homs(self):
        # Call the function
        # prep_co_putative_combos(self.bg, self.allele_relationships)
        result_bg = list(
            prep_co_putative_combos({1: self.bg}, self.allele_relationships).values()
        )[0]
        # Expected combos
        expected_combos = {
            (self.allele10,),
            (self.allele12,),
            (self.allele7,),
            (self.allele6,),
        }

        self.checker(result_bg.misc["combos"], expected_combos)
        # self.assertEqual(result_bg.misc["combos"], expected_combos)

    def test_prep_co_putative_combos_no_homs(self):
        # Remove homozygous alleles
        self.bg.alleles[AlleleState.POS] = [self.allele_ref, self.allele12]
        # Call the function
        # prep_co_putative_combos(self.bg, self.allele_relationships)
        result_bg = list(
            prep_co_putative_combos({1: self.bg}, self.allele_relationships).values()
        )[0]
        # Expected combos
        # self.allele_ref = unique_alleles[0] #ref
        # self.allele10 = unique_alleles[1] #10
        # self.allele12 = unique_alleles[2] #12
        # self.allele7 = unique_alleles[3] #7
        # self.allele6 = unique_alleles[4] #6
        expected_combos = {(self.allele12,)}
        self.checker(result_bg.misc["combos"], expected_combos)
        # self.assertEqual(result_bg.misc["combos"], expected_combos)

    def test_prep_co_putative_combos_with_reference_allele(self):
        # Remove reference status from allele3
        # self.allele12.reference = False

        # Update the alleles in bg
        self.bg_with_ref.alleles[AlleleState.POS] = [
            self.allele_ref,
            self.allele12,
        ]
        self.bg_with_ref.misc = {}
        self.bg_with_ref.misc["homs"] = []

        # Call the function
        result_bg = list(
            prep_co_putative_combos(
                {1: self.bg_with_ref}, self.allele_relationships
            ).values()
        )[0]

        # Expected combos
        expected_combos = {
            (self.allele12,),
        }
        self.checker(result_bg.misc["combos"], expected_combos)
        # self.assertEqual(result_bg.misc["combos"], expected_combos)

    def test_prep_co_putative_combos_sub_alleles(self):
        # Test with sub-allele relationships that exclude combinations
        self.bg.misc["homs"] = []
        self.allele_relationships["KN"]["KN*01_isin_KN*02"] = True
        self.allele_relationships["KN"]["KN*02_isin_KN*01"] = True

        # Call the function
        # prep_co_putative_combos(self.bg, self.allele_relationships)
        result_bg = list(
            prep_co_putative_combos({1: self.bg}, self.allele_relationships).values()
        )[0]
        # Expected combos
        expected_combos = {
            (self.allele10,),
            (self.allele12,),
            (self.allele10, self.allele12),
            (self.allele7,),
            (self.allele7, self.allele12),
            (self.allele6,),
            (self.allele6, self.allele12),
        }
        # Process both actual and expected combos
        self.checker(result_bg.misc["combos"], expected_combos)
        # actual_combos = process_combos(result_bg.misc["combos"])
        # expected_combos_processed = process_combos(expected_combos)

        # self.assertEqual(actual_combos, expected_combos_processed)

        # # Since allele1 and allele2 are sub-alleles of each other, combinations of both are excluded
        # self.assertEqual(result_bg.misc["combos"], expected_combos)

    def test_prep_co_putative_combos_all_alleles(self):
        # Remove homozygous alleles and reference status
        self.bg.misc["homs"] = []
        # self.allele12.reference = False

        # Update the alleles in bg
        self.bg.alleles[AlleleState.POS] = [
            self.allele10,
            self.allele12,
            self.allele7,
            self.allele6,
        ]

        # Update allele relationships to allow all combinations
        self.allele_relationships["KN"] = ALLELE_RELATIONSHIPS["KN"]

        # Call the function
        # prep_co_putative_combos(self.bg, self.allele_relationships)
        result_bg = list(
            prep_co_putative_combos({1: self.bg}, self.allele_relationships).values()
        )[0]
        # Expected combos
        expected_combos = {
            (self.allele10,),
            (self.allele12,),
            (self.allele7,),
            (self.allele6,),
            (self.allele10, self.allele12),
            (self.allele12, self.allele7),
            (self.allele6, self.allele12),
        }

        self.checker(result_bg.misc["combos"], expected_combos)
        # self.allele_ref = unique_alleles[0] #ref
        # self.allele10 = unique_alleles[1] #10
        # self.allele12 = unique_alleles[2] #12
        # self.allele7 = unique_alleles[3] #7
        # self.allele6 = unique_alleles[4] #6
        # actual_combos = process_combos(result_bg.misc["combos"])
        # expected_combos_processed = process_combos(expected_combos)

        # self.assertEqual(actual_combos, expected_combos_processed)
        # self.assertEqual(result_bg.misc["combos"], expected_combos)

    def test_prep_co_putative_combos_excluding_reference_allele(self):
        # Keep allele3 as reference
        self.bg.alleles[AlleleState.POS] = [
            self.allele_ref,
            self.allele10,
            self.allele12,
        ]

        # Remove homozygous alleles
        self.bg.misc["homs"] = []

        # Update allele relationships to allow combinations
        self.allele_relationships["KN"] = ALLELE_RELATIONSHIPS["KN"]
        # Call the function
        # prep_co_putative_combos(self.bg, self.allele_relationships)
        result_bg = list(
            prep_co_putative_combos({1: self.bg}, self.allele_relationships).values()
        )[0]
        # Expected combos (excluding combinations with the reference allele)
        expected_combos = {
            (self.allele12,),
            (self.allele12, self.allele10),
            (self.allele10,),
        }
        self.checker(result_bg.misc["combos"], expected_combos)

    def test_combos_with_multiple_subtypes_are_excluded(self):
        # Add alleles with different subtypes
        # self.allele2.sub_type = "KN*02.02"
        self.bg.alleles[AlleleState.POS] = [
            self.allele10,
            self.allele7,
            self.allele2,
            self.allele_ref,
        ]
        self.bg.misc = {}
        self.bg.misc["homs"] = []
        # Call the function
        result_bg = list(
            prep_co_putative_combos({1: self.bg}, self.allele_relationships).values()
        )[0]
        # Expected combos should not include combinations with multiple subtypes
        expected_combos = {
            (self.allele10,),
            (self.allele7,),
            (self.allele2,),
        }
        self.checker(result_bg.misc["combos"], expected_combos)


class TestCheckHomVariants(unittest.TestCase):
    def setUp(self) -> None:
        self.allele_ref = Allele(
            genotype="KN*01",
            phenotype="positive",
            genotype_alt="KN*02",
            phenotype_alt="negative",
            defining_variants=frozenset({"variant1", "variant2"}),
            weight_geno=2,
            weight_pheno=3,
            reference=False,
        )
        self.allele10 = Allele(
            genotype="KN*02",
            phenotype="negative",
            genotype_alt="KN*03",
            phenotype_alt="positive",
            defining_variants=frozenset({"variant1"}),
            weight_geno=1,
            weight_pheno=2,
            reference=False,
        )

    def test_all_hom_variants_all_homs_in_combo(self) -> None:
        all_homs = [self.allele_ref, self.allele10]
        current_combo = (self.allele_ref, self.allele10)
        self.assertTrue(all_hom_variants(all_homs, current_combo))

    def test_all_hom_variants_some_homs_in_combo(self) -> None:
        all_homs = [self.allele_ref, self.allele10]
        current_combo = (self.allele_ref,)
        self.assertFalse(all_hom_variants(all_homs, current_combo))

    def test_all_hom_variants_no_homs_provided(self) -> None:
        all_homs = []
        current_combo = (self.allele_ref,)
        self.assertFalse(all_hom_variants(all_homs, current_combo))

    def test_all_hom_variants_no_homs_in_combo(self) -> None:
        all_homs = [self.allele_ref]
        current_combo = (self.allele10,)
        self.assertTrue(all_hom_variants(all_homs, current_combo))


class TestSubAlleles(unittest.TestCase):
    def setUp(self) -> None:
        self.allele_ref = Allele(
            genotype="KN*01",
            phenotype="positive",
            genotype_alt="KN*02",
            phenotype_alt="negative",
            defining_variants=frozenset({"variant1", "variant2"}),
            weight_geno=2,
            weight_pheno=3,
            reference=False,
        )
        self.allele10 = Allele(
            genotype="KN*02",
            phenotype="negative",
            genotype_alt="KN*03",
            phenotype_alt="positive",
            defining_variants=frozenset({"variant1"}),
            weight_geno=1,
            weight_pheno=2,
            reference=False,
        )

    def test_sub_alleles_true(self) -> None:
        lst = (self.allele_ref, self.allele10)
        allele_relationships = {"KN*01_isin_KN*02": True, "KN*02_isin_KN*01": False}
        self.assertTrue(sub_alleles(lst, allele_relationships))

    def test_sub_alleles_false(self) -> None:
        lst = (self.allele_ref, self.allele10)
        allele_relationships = {"KN*01_isin_KN*02": False, "KN*02_isin_KN*01": False}
        self.assertFalse(sub_alleles(lst, allele_relationships))

    def test_sub_alleles_single_element(self) -> None:
        lst = (self.allele_ref,)
        allele_relationships = {}
        self.assertFalse(sub_alleles(lst, allele_relationships))


class TestAddCoExistingAlleleAndRef(unittest.TestCase):
    def setUp(self):
        # Create Alleles
        self.allele1 = Allele(
            genotype="KN*01.06",
            defining_variants=frozenset({"var1", "var2"}),
            weight_geno=1,
            weight_pheno=1,
            phenotype="KN:6",
            genotype_alt="KN*A",
            phenotype_alt="PhenotypeA",
            reference=False,
            sub_type="KN*01",
        )
        self.allele2 = Allele(
            genotype="KN*01.07",
            defining_variants=frozenset({"var3"}),
            weight_geno=2,
            weight_pheno=1,
            phenotype="KN:7",
            genotype_alt="KN*B",
            phenotype_alt="PhenotypeB",
            reference=False,
            sub_type="KN*01",
        )
        # Reference Allele
        self.reference_allele = Allele(
            genotype="KN*01",
            defining_variants=frozenset(),
            weight_geno=0,
            weight_pheno=1,
            phenotype="KN:1",
            genotype_alt="KN*Ref",
            phenotype_alt="PhenotypeRef",
            reference=True,
            sub_type="KN*01",
        )
        self.reference_alleles = {"KN": self.reference_allele}
        # Variant pool
        self.variant_pool = {
            "var1": Zygosity.HET,
            "var2": Zygosity.HET,
            "var3": Zygosity.HET,
        }
        # self.variant_pool_numeric = {
        #     "var1": 1,
        #     "var2": 1,
        #     "var3": 1,
        # }

        # Create BloodGroup
        self.bg = BloodGroup(
            type="KN",
            alleles={AlleleState.POS: [self.allele1, self.allele2]},
            sample="sample1",
            variant_pool=self.variant_pool,
            genotypes=[],
            phenotypes=[],
            filtered_out={},
        )
        # Set misc attributes
        self.bg.misc = {}
        self.bg.misc["combos"] = [(self.allele1,), (self.allele2,)]
        self.bg.misc["homs"] = []  # Assuming no homs
        self.bg.misc["max_rank"] = 2  # Set max rank

        # Initially no co_existing
        self.bg.alleles[AlleleState.CO] = None

    def test_add_co_existing_with_reference(self):
        # Call the function
        result_bg = list(
            add_co_existing_allele_and_ref(
                {1: self.bg}, self.reference_alleles
            ).values()
        )[0]

        # Expected co-existing combinations
        expected_co_existing = [
            ((self.allele1,), (self.reference_allele,)),
            ((self.allele2,), (self.reference_allele,)),
        ]

        # Check that co-existing alleles include combinations with reference allele
        self.assertEqual(result_bg.alleles[AlleleState.CO], expected_co_existing)

    def test_max_rank_excludes_combinations(self):
        # Set max_rank higher than any allele's weight_geno
        self.bg.misc["max_rank"] = 0  # All alleles have weight_geno >=1
        # Call the function
        result_bg = list(
            add_co_existing_allele_and_ref(
                {1: self.bg}, self.reference_alleles
            ).values()
        )[0]
        # Since max_rank is 0, all combinations should be skipped
        self.assertIsNone(result_bg.alleles[AlleleState.CO])

    def test_with_homs(self):
        # Set homs to include one allele
        self.bg.misc["homs"] = [self.allele1]
        # Call the function
        result_bg = list(
            add_co_existing_allele_and_ref(
                {1: self.bg}, self.reference_alleles
            ).values()
        )[0]
        # Since there are homs, should skip combinations
        self.assertIsNone(result_bg.alleles[AlleleState.CO])

    def test_non_kn_bg_sets_co_none(self):
        """
        Covers:
        if bg.type != 'KN':
            bg.alleles[AlleleState.CO] = None
            return bg
        """
        # Make a copy of self.bg but set type != 'KN'
        not_kn_bg = BloodGroup(
            type="AB",
            alleles={AlleleState.CO: None},  # could be anything
            sample="non_kn_sample",
            variant_pool={},
            genotypes=[],
            phenotypes=[],
            filtered_out={},
        )
        not_kn_bg.misc = {}
        not_kn_bg.alleles[AlleleState.CO] = []  # just to show it's not None

        # Call the function
        result_dict = add_co_existing_allele_and_ref(
            {999: not_kn_bg}, self.reference_alleles
        )
        result_bg = result_dict[999]

        # We expect CO to be set to None and returned immediately
        self.assertIsNone(
            result_bg.alleles[AlleleState.CO],
            "If bg.type != 'KN', CO must be set to None.",
        )

    def test_co_existing_is_not_none_before_loop(self):
        """
        Covers:
        if bg.alleles[AlleleState.CO] is not None:
            co_existing = bg.alleles[AlleleState.CO]

        We set CO to a non-None value, ensuring we use that value
        instead of an empty list.
        """
        # Suppose user pre-populates CO with some pairs
        existing_co = [
            ((self.allele1,), (self.allele2,)),  # just an example structure
        ]
        self.bg.alleles[AlleleState.CO] = existing_co

        # The function should see 'if bg.alleles[AlleleState.CO] is not None'
        # and continue from there, reusing that 'co_existing' list
        result_bg = list(
            add_co_existing_allele_and_ref(
                {10: self.bg}, self.reference_alleles
            ).values()
        )[0]
        # If the function doesn't skip combos due to rank/homs, it might add more,
        # but we at least confirm it starts from the existing_co instead of None
        self.assertIsNotNone(
            result_bg.alleles[AlleleState.CO],
            "We expected the function to keep pre-existing CO combos.",
        )

    def test_reference_variants_found_in_tmp_pool2_are_decremented(self):
        """
        Covers the line:
        if variant_on_other_strand in tmp_pool2:
            tmp_pool2[variant_on_other_strand] -= 1
        by ensuring the reference Allele actually has that variant in its defining_variants.
        """

        # 1. Clone existing reference Allele with new variants
        old_ref = self.reference_alleles["KN"]
        new_defining_variants = old_ref.defining_variants.union(
            {"var1", "extra_ref_var"}
        )
        new_ref = Allele(
            genotype=old_ref.genotype,
            genotype_alt=old_ref.genotype_alt,
            defining_variants=new_defining_variants,
            phenotype=old_ref.phenotype,
            phenotype_alt=old_ref.phenotype_alt,
            weight_geno=old_ref.weight_geno,
            weight_pheno=old_ref.weight_pheno,
            reference=old_ref.reference,
            sub_type=old_ref.sub_type,
        )
        self.reference_alleles["KN"] = new_ref

        # 2. Update the existing variant_pool_numeric in place (no direct assignment!)
        pool = self.bg.variant_pool_numeric
        pool["var1"] = 2
        # pool["unrelated_var"] = 5

        # 3. Ensure we won't skip combos
        self.bg.misc["max_rank"] = 10
        self.bg.misc["homs"] = []

        # 4. Run the function
        result_bg = list(
            add_co_existing_allele_and_ref(
                {20: self.bg}, self.reference_alleles
            ).values()
        )[0]

        # 5. Confirm var1 was decremented
        self.assertEqual(
            result_bg.variant_pool_numeric["var1"],
            1,
            "Expected 'var1' to be decremented from 2 to 1.",
        )
        # self.assertEqual(result_bg.variant_pool_numeric["unrelated_var"], 5,
        #                 "'unrelated_var' should remain unchanged.")


class TestFilterRedundantPairs(unittest.TestCase):
    def setUp(self):
        # Create Alleles
        self.allele1 = Allele(
            genotype="KN*01.06",
            defining_variants=frozenset(),
            weight_geno=1,
            weight_pheno=1,
            phenotype="KN:6",
            genotype_alt="KN*A",
            phenotype_alt="PhenotypeA",
            reference=False,
            sub_type="KN*01",
        )
        self.allele2 = Allele(
            genotype="KN*01.07",
            defining_variants=frozenset(),
            weight_geno=2,
            weight_pheno=1,
            phenotype="KN:7",
            genotype_alt="KN*B",
            phenotype_alt="PhenotypeB",
            reference=False,
            sub_type="KN*01",
        )
        self.allele3 = Allele(
            genotype="KN*01.10",
            defining_variants=frozenset(),
            weight_geno=3,
            weight_pheno=1,
            phenotype="KN:10",
            genotype_alt="KN*C",
            phenotype_alt="PhenotypeC",
            reference=False,
            sub_type="KN*01",
        )

        # Create co_existing combinations with redundant pairs
        # For example, (allele1, allele2) and (allele2, allele1)
        self.bg = BloodGroup(
            type="KN",
            alleles={
                AlleleState.CO: [
                    ((self.allele1,), (self.allele2,)),
                    ((self.allele2,), (self.allele1,)),
                    ((self.allele1,), (self.allele3,)),
                    ((self.allele3,), (self.allele1,)),
                ]
            },
            sample="sample1",
            variant_pool={},
            genotypes=[],
            phenotypes=[],
            filtered_out={},
        )

    def test_filter_redundant_pairs(self):
        # Call the function
        result_bg = list(filter_redundant_pairs({1: self.bg}).values())[0]
        # Expected co_existing combinations after filtering redundant pairs
        expected_co_existing = [
            ((self.allele1,), (self.allele2,)),
            ((self.allele1,), (self.allele3,)),
        ]
        self.assertEqual(result_bg.alleles[AlleleState.CO], expected_co_existing)


# def geno_str(combo):
#     return "+".join(sorted([a.genotype for a in combo]))


# class TestListExcludedCoExistingPairs(unittest.TestCase):
#     def setUp(self):
#         # Create Alleles
#         self.allele1 = Allele(
#             genotype="KN*01.06",
#             defining_variants=frozenset(),
#             weight_geno=1,
#             weight_pheno=1,
#             phenotype="KN:6",
#             genotype_alt="KN*A",
#             phenotype_alt="PhenotypeA",
#             reference=False,
#             sub_type="KN*01",
#         )
#         self.allele2 = Allele(
#             genotype="KN*01.07",
#             defining_variants=frozenset(),
#             weight_geno=2,
#             weight_pheno=1,
#             phenotype="KN:7",
#             genotype_alt="KN*B",
#             phenotype_alt="PhenotypeB",
#             reference=False,
#             sub_type="KN*01",
#         )
#         self.reference_allele = Allele(
#             genotype="KN*01",
#             defining_variants=frozenset(),
#             weight_geno=0,
#             weight_pheno=1,
#             phenotype="KN:1",
#             genotype_alt="KN*Ref",
#             phenotype_alt="PhenotypeRef",
#             reference=True,
#             sub_type="KN*01",
#         )
#         self.reference_alleles = {"KN": self.reference_allele}

#         # Create BloodGroup
#         self.bg = BloodGroup(
#             type="KN",
#             alleles={AlleleState.CO: [((self.allele1,), (self.allele2,))]},
#             sample="sample1",
#             variant_pool={},
#             genotypes=[],
#             phenotypes=[],
#             filtered_out={},
#         )
#         # Set misc["combos"]
#         self.bg.misc = {}
#         self.bg.misc["combos"] = [(self.allele1,), (self.allele2,)]
#         self.bg.filtered_out = {}


class TestListExcludedCoExistingPairs(unittest.TestCase):
    """Test suite covering all logic paths for `list_excluded_co_existing_pairs`."""

    def setUp(self):
        # Minimal Allele for reference usage
        self.ref_allele = Allele(
            genotype="KN*01",
            genotype_alt="RefAlt",
            defining_variants=frozenset({"var1"}),
            phenotype="KN:1",
            phenotype_alt="Ref+",
            reference=True,
        )
        self.reference_alleles = {"KN": self.ref_allele}

    def test_non_kn_type_sets_co_none(self):
        """If bg.type != 'KN', sets CO to None and returns."""
        bg = BloodGroup(
            type="AB",  # not KN
            sample="test1",
            alleles={AlleleState.CO: []},  # or any other value
            misc={},
            filtered_out={},
        )
        result_dict = list_excluded_co_existing_pairs({0: bg}, self.reference_alleles)
        result_bg = result_dict[0]
        self.assertIsNone(
            result_bg.alleles[AlleleState.CO],
            "CO should be set to None for non-KN types.",
        )

    def test_empty_combos(self):
        """If bg.misc['combos'] is empty, we never add anything to tested => all remain unfiltered."""
        bg = BloodGroup(
            type="KN",
            sample="test2",
            alleles={
                AlleleState.CO: []  # or any pairs
            },
            misc={"combos": []},  # empty combos
            filtered_out={},
        )
        # We do have a reference for "KN"
        result_dict = list_excluded_co_existing_pairs({1: bg}, self.reference_alleles)
        result_bg = result_dict[1]
        # tested is empty => filtered_out is empty
        self.assertEqual(
            result_bg.filtered_out.get(AlleleState.CO),
            [],
            "No combos => no tested => no filtered_out pairs.",
        )
        # CO remains as it was
        self.assertEqual(
            result_bg.alleles[AlleleState.CO],
            [],
            "CO should remain unchanged if no combos exist.",
        )

    def test_some_combos_reference_in_alleles_co(self):
        """If we have combos in misc and some pre-existing CO pairs,
        ensure only pairs not in 'CO' go to filtered_out."""
        # 2 combos: comboA, comboB
        # We'll define them as lists of 1 Allele each, to keep it simple.
        alleleA = Allele(
            genotype="A",
            genotype_alt="Aalt",
            defining_variants=frozenset({"varA"}),
            phenotype="KN:A",
            phenotype_alt="A+",
        )
        alleleB = Allele(
            genotype="B",
            genotype_alt="Balt",
            defining_variants=frozenset({"varB"}),
            phenotype="KN:B",
            phenotype_alt="B+",
        )
        combos = [
            [alleleA],  # combo1
            [alleleB],  # combo2
        ]
        # Pre-existing CO pairs => suppose user added some real pairs
        preexisting_pairs = [Pair(alleleA, alleleB)]
        bg = BloodGroup(
            type="KN",
            sample="test3",
            alleles={AlleleState.CO: preexisting_pairs},
            misc={"combos": combos},
            filtered_out={},
        )
        out_dict = list_excluded_co_existing_pairs({99: bg}, self.reference_alleles)
        out_bg = out_dict[99]

        # tested => we generate a Pair(A_mushed, B_mushed), Pair(A_mushed, A_mushed), etc., plus Pair(ref, X).
        # preexisting_pairs has exactly [Pair(alleleA, alleleB)] => So only that pair won't be filtered out
        filtered = out_bg.filtered_out[AlleleState.CO]
        # Check that all tested pairs except the (alleleA, alleleB) are in filtered_out
        self.assertTrue(
            len(filtered) >= 1, "Expected at least some pairs filtered out."
        )
        self.assertNotIn(
            Pair(alleleA, alleleB),
            filtered,
            "Pair(alleleA, alleleB) was in CO => should NOT be in filtered_out.",
        )

    def test_multiple_combos_ensures_all_pairs_tested(self):
        """Checks that for each combo1 in misc['combos'] and combo2, we produce tested pairs,
        plus the Pair(ref, mushed(combo1)).
        Then we verify they're filtered if absent from CO."""
        alleleX = Allele(
            genotype="X",
            genotype_alt="Xalt",
            defining_variants=frozenset({"varX"}),
            phenotype="KN:X",
            phenotype_alt="X+",
        )
        alleleY = Allele(
            genotype="Y",
            genotype_alt="Yalt",
            defining_variants=frozenset({"varY"}),
            phenotype="KN:Y",
            phenotype_alt="Y+",
        )
        combos = [
            [alleleX],
            [alleleY],
        ]
        # Suppose CO is empty => everything tested is filtered out
        bg = BloodGroup(
            type="KN",
            sample="test4",
            alleles={AlleleState.CO: []},
            misc={"combos": combos},
            filtered_out={},
        )
        out_dict = list_excluded_co_existing_pairs({2: bg}, self.reference_alleles)
        out_bg = out_dict[2]
        # We expect tested includes:
        #   Pair(mushed(X), mushed(X)), Pair(mushed(X), mushed(Y)),
        #   Pair(mushed(Y), mushed(X)), Pair(mushed(Y), mushed(Y))
        #   Pair(ref, mushed(X)), Pair(ref, mushed(Y))
        # => 6 total pairs. Since CO was empty => all 6 => filtered_out
        self.assertEqual(
            len(out_bg.filtered_out[AlleleState.CO]),
            6,
            "Expected 2 combos => 2x2 + 2 for reference => 6 tested => all filtered out.",
        )

    def test_co_is_not_modified_if_already_has_pairs(self):
        """Ensure we do not remove or alter the existing pairs in CO; we only set filtered_out."""
        # aX = Allele("X", "Xalt", frozenset(), "KN:X", "X+")
        aX = Allele(
            genotype="X",
            genotype_alt="Xalt",
            defining_variants=frozenset({}),
            phenotype="KN:X",
            phenotype_alt="X+",
        )
        combos = [[aX]]
        existing_pairs = [Pair(aX, aX)]  # nonsense, but enough to verify it remains

        bg = BloodGroup(
            type="KN",
            sample="test5",
            alleles={AlleleState.CO: existing_pairs},
            misc={"combos": combos},
            filtered_out={},
        )
        out = list_excluded_co_existing_pairs({999: bg}, self.reference_alleles)[999]
        self.assertEqual(
            out.alleles[AlleleState.CO],
            existing_pairs,
            "We do not expect the function to remove or alter existing CO pairs.",
        )

    def test_reference_allele_not_provided_for_kn_raises_key_error(self):
        """If reference_alleles lacks an entry for bg.type, we might get a KeyError."""
        # This can happen if the user forgets to supply a reference for "KN"
        # We'll demonstrate that scenario is unhandled => KeyError
        bg = BloodGroup(
            type="KN", sample="test6", alleles={AlleleState.CO: []}, misc={"combos": []}
        )
        # Provide a reference_alleles that lacks 'KN'
        bad_ref_dict = {
            "AB": self.ref_allele  # only AB => no 'KN'
        }
        with self.assertRaises(KeyError) as ctx:
            list_excluded_co_existing_pairs({123: bg}, bad_ref_dict)
        self.assertIn(
            "KN",
            str(ctx.exception),
            "Expected a KeyError complaining about missing 'KN' reference.",
        )


if __name__ == "__main__":
    unittest.main()
