import unittest
from collections import defaultdict
from rbceq2.core_logic.utils import Zygosity
from rbceq2.core_logic.alleles import Allele, BloodGroup, Pair
from rbceq2.core_logic.constants import AlleleState
from rbceq2.filters.geno import (
    ABO_cant_pair_with_ref_cuz_261delG_HET,
    cant_pair_with_ref_cuz_SNPs_must_be_on_other_side,
    cant_pair_with_ref_cuz_trumped,
    filter_HET_pairs_by_weight,
    filter_pairs_by_context,
    filter_pairs_on_antithetical_zygosity,
    flatten_alleles,
    split_pair_by_ref,
    antithetical_modifying_SNP_is_HOM
)



class TestFlattenAlleles(unittest.TestCase):
    def test_flatten_alleles(self):
        allele1 = Allele(
            genotype="A1",
            phenotype="M",
            genotype_alt=".",
            phenotype_alt=".",
            defining_variants=frozenset({"var1"}),
            null=False,
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
            null=False,
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
            null=False,
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
            null=False,
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
            null=False,
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
            null=False,
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
            null=False,
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
            null=False,
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
            null=False,
            weight_geno=2,
            reference=False,
            sub_type="subtype2",
        )
        pair = Pair(allele1=allele1, allele2=allele2)
        with self.assertRaises(ValueError):
            split_pair_by_ref(pair)







class TestFilterPairsOnAntitheticalZygosity(unittest.TestCase):
    def setUp(self):
        allele2 = Allele(
            genotype="FY*02",
            phenotype="FY:2",
            genotype_alt="FY*B",
            phenotype_alt="Fy(b+)",
            defining_variants=frozenset({"1:159175354_G_A"}),
            null=False,
            weight_geno=1000,
            reference=False,
            sub_type="FY*02",
        )
        allele3 = Allele(
            genotype="FY*01",
            phenotype="FY:1",
            genotype_alt="FY*A",
            phenotype_alt="Fy(a+)",
            defining_variants=frozenset({"1:159175354_ref"}),
            null=False,
            weight_geno=1000,
            reference=True,
            sub_type="FY*01",
        )
        allele4 = Allele(
            genotype="FY*01N.01",
            phenotype=".",
            genotype_alt=".",
            phenotype_alt="Fy(a-b-)",
            defining_variants=frozenset({"1:159175354_ref", "1:159174683_T_C"}),
            null=False,
            weight_geno=7,
            reference=False,
            sub_type="FY*01",
        )
        # have to have both subtypes in pair
        self.pair1 = Pair(allele1=allele2, allele2=allele4)  # ok
        self.pair2 = Pair(allele1=allele3, allele2=allele4)  # not ok
        self.bg = BloodGroup(
            type="FY",
            alleles={AlleleState.NORMAL: [self.pair1, self.pair2]},
            sample="013Kenya",
            variant_pool={"1:159175354_G_A": Zygosity.HET},
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
            variant_pool={"1:159175354_G_A": Zygosity.HET},
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
            null=False,
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
            null=False,
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
            null=False,
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
            antithetical_modifying_SNP_is_HOM(
                {1: bg}, self.antitheticals
            ).values()
        )[0]
        self.assertTrue(self.pair2 not in filtered_bg.alleles[AlleleState.NORMAL])
        self.assertTrue(
            self.pair2
            in filtered_bg.filtered_out["antithetical_modifying_SNP_is_HOM"]
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
            antithetical_modifying_SNP_is_HOM(
                {1: bg}, self.antitheticals
            ).values()
        )[0]
        self.assertTrue(self.pair1 in filtered_bg.alleles[AlleleState.NORMAL])
        self.assertTrue(
            self.pair1
            not in filtered_bg.filtered_out[
                "antithetical_modifying_SNP_is_HOM"
            ]
        )

    def test_empty_normal_list(self):
        bg = BloodGroup(
            type="LU",
            alleles={AlleleState.NORMAL: []},
            sample="128",
            variant_pool={"19:45315445_G_A": Zygosity.HET},
            filtered_out=defaultdict(list),
        )
        filtered_bg = list(
            antithetical_modifying_SNP_is_HOM(
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
            null=False,
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
            null=False,
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
            null=False,
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
            null=False,
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
            null=False,
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
            null=False,
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
            null=False,
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
            null=False,
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
            null=False,
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
            null=False,
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
            null=False,
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
            null=False,
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
            null=False,
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
            null=False,
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
            null=False,
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
            null=False,
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
            null=False,
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
            null=False,
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
            null=False,
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
            null=False,
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
            null=False,
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




if __name__ == "__main__":
    unittest.main()
