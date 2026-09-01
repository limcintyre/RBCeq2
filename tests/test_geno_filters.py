import unittest
from collections import defaultdict
from rbceq2.core_logic.utils import Zygosity
from rbceq2.core_logic.alleles import Allele, BloodGroup, Pair
from rbceq2.core_logic.constants import AlleleState, UNDETERMINED_SLOT
from rbceq2.core_logic.data_procesing import get_genotypes
from rbceq2.filters.geno import (
    ABO_cant_pair_with_ref_cuz_261delG_HET,
    cant_name_second_slot_cuz_ref_impossible,
    cant_name_second_slot_cuz_shared_variant_has_too_few_copies,
    cant_pair_with_ref_cuz_shared_variant_has_too_few_copies,
    pool_cant_supply_both,
    ref_slot_is_impossible,
    cant_pair_with_ref_cuz_SNPs_must_be_on_other_side,
    cant_pair_with_ref_cuz_trumped,
    filter_HET_pairs_by_weight,
    filter_pairs_by_context,
    filter_pairs_on_antithetical_zygosity,
    flatten_alleles,
    split_pair_by_ref,
    antithetical_modifying_SNP_is_HOM,
    cant_have_2_non_ref_alleles_cuz_only_1_gene_copy,
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
            variant_pool={
                "1:159175354_G_A": Zygosity.HET,
                "1:159175354_ref": Zygosity.HET,
                "1:159174683_T_C": Zygosity.HET,
            },
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
            antithetical_modifying_SNP_is_HOM({1: bg}, self.antitheticals).values()
        )[0]
        self.assertTrue(self.pair2 not in filtered_bg.alleles[AlleleState.NORMAL])
        self.assertTrue(
            self.pair2 in filtered_bg.filtered_out["antithetical_modifying_SNP_is_HOM"]
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
            antithetical_modifying_SNP_is_HOM({1: bg}, self.antitheticals).values()
        )[0]
        self.assertTrue(self.pair1 in filtered_bg.alleles[AlleleState.NORMAL])
        self.assertTrue(
            self.pair1
            not in filtered_bg.filtered_out["antithetical_modifying_SNP_is_HOM"]
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
            antithetical_modifying_SNP_is_HOM({1: bg}, self.antitheticals).values()
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


class TestPoolCantSupplyBoth(unittest.TestCase):
    """The arithmetic shared by ABO_cant_pair_with_ref_cuz_261delG_HET and
    cant_pair_with_ref_cuz_shared_variant_has_too_few_copies."""

    def setUp(self):
        self.ref = Allele(
            genotype="GYPA*02",
            phenotype=".",
            genotype_alt=".",
            phenotype_alt=".",
            defining_variants=frozenset({"4:144120554_ref", "4:144120555_ref", "4:144120567_ref"}),
            null=False,
            weight_geno=1000,
            reference=True,
            sub_type="GYPA*02",
        )
        self.other = Allele(
            genotype="GYPA*08",
            phenotype=".",
            genotype_alt=".",
            phenotype_alt=".",
            defining_variants=frozenset({"4:144120554_ref", "4:144120555_ref", "4:144120567_A_G"}),
            null=False,
            weight_geno=1000,
            reference=False,
            sub_type="GYPA*08",
        )

    def test_shared_variant_with_one_copy_cannot_supply_both(self):
        pool = {
            "4:144120554_ref": 1,
            "4:144120555_ref": 2,
            "4:144120567_A_G": 1,
            "4:144120567_ref": 1,
        }
        self.assertTrue(pool_cant_supply_both(self.ref, self.other, pool))

    def test_shared_variant_with_two_copies_can_supply_both(self):
        pool = {
            "4:144120554_ref": 2,
            "4:144120555_ref": 2,
            "4:144120567_A_G": 1,
            "4:144120567_ref": 1,
        }
        self.assertFalse(pool_cant_supply_both(self.ref, self.other, pool))

    def test_pool_is_not_spent_by_the_call(self):
        pool = {
            "4:144120554_ref": 1,
            "4:144120555_ref": 2,
            "4:144120567_A_G": 1,
            "4:144120567_ref": 1,
        }
        pool_cant_supply_both(self.ref, self.other, pool)
        self.assertEqual(pool["4:144120554_ref"], 1)

    def test_reference_variant_absent_from_pool_is_not_spent(self):
        """A reference reaches here from the database, so its variants need not be in
        the pool. Absence is cant_name_second_slot_cuz_ref_impossible's claim, not
        this one's, and must not raise."""
        pool = {"4:144120554_ref": 1, "4:144120555_ref": 2, "4:144120567_A_G": 1}
        ref_needing_absent = Allele(
            genotype="GYPA*02",
            phenotype=".",
            genotype_alt=".",
            phenotype_alt=".",
            defining_variants=frozenset({"4:144120554_ref", "4:144120555_ref", "4:144120567_ref"}),
            null=False,
            weight_geno=1000,
            reference=True,
            sub_type="GYPA*02",
        )
        self.assertTrue(
            pool_cant_supply_both(ref_needing_absent, self.other, pool)
        )


class TestCantPairWithRefCuzSharedVariantHasTooFewCopies(unittest.TestCase):
    """HG03097 RHCE, short read. RHCE*01 is a reference defined partly by an alternate
    at a lane locus, and RHCE*01.20.01 is RHCE*01 plus one more variant, so both need
    1:25420739_G_C - which is heterozygous."""

    def setUp(self):
        self.ref = Allele(
            genotype="RHCE*01",
            phenotype=".",
            genotype_alt=".",
            phenotype_alt=".",
            defining_variants=frozenset({"1:25390874_ref", "1:25408711_ref", "1:25420739_G_C"}),
            null=False,
            weight_geno=1000,
            reference=True,
            sub_type="RHCE*01",
        )
        self.allele_01_01 = Allele(
            genotype="RHCE*01.01",
            phenotype=".",
            genotype_alt=".",
            phenotype_alt=".",
            defining_variants=frozenset({"1:25390874_ref", "1:25408711_ref", "1:25420739_ref"}),
            null=False,
            weight_geno=1000,
            reference=False,
            sub_type="RHCE*01",
        )
        self.allele_20_01 = Allele(
            genotype="RHCE*01.20.01",
            phenotype=".",
            genotype_alt=".",
            phenotype_alt=".",
            defining_variants=frozenset({"1:25390817_G_C", "1:25390874_ref", "1:25408711_ref", "1:25420739_G_C"}),
            null=False,
            weight_geno=1000,
            reference=False,
            sub_type="RHCE*01",
        )

        self.impossible = Pair(allele1=self.ref, allele2=self.allele_20_01)
        self.possible = Pair(allele1=self.ref, allele2=self.allele_01_01)
        self.no_reference = Pair(allele1=self.allele_01_01, allele2=self.allele_20_01)

        self.bg = BloodGroup(
            type="RHCE",
            alleles={
                AlleleState.NORMAL: [
                    self.impossible,
                    self.possible,
                    self.no_reference,
                ]
            },
            sample="HG03097",
            variant_pool={
                "1:25390817_G_C": Zygosity.HET,
                "1:25390874_ref": Zygosity.HOM,
                "1:25408711_ref": Zygosity.HOM,
                "1:25420739_G_C": Zygosity.HET,
                "1:25420739_ref": Zygosity.HET,
            },
            filtered_out=defaultdict(list),
        )
        cant_pair_with_ref_cuz_shared_variant_has_too_few_copies({1: self.bg})

    def test_pair_needing_two_copies_of_one_variant_is_removed(self):
        self.assertTrue(self.impossible not in self.bg.alleles[AlleleState.NORMAL])
        self.assertTrue(
            self.impossible
            in self.bg.filtered_out[
                "cant_pair_with_ref_cuz_shared_variant_has_too_few_copies"
            ]
        )

    def test_pair_whose_shared_variants_are_hom_is_kept(self):
        self.assertTrue(self.possible in self.bg.alleles[AlleleState.NORMAL])

    def test_pair_without_a_reference_is_left_to_pair_can_exist(self):
        self.assertTrue(self.no_reference in self.bg.alleles[AlleleState.NORMAL])


class TestCantPairWithRefCuzSharedVariantIsRefToken(unittest.TestCase):
    """HG00436 GYPA, long read. Neither allele is defined by an alternate, and the
    shared variant is a '_ref' token - which says the reference base is on one
    chromosome and 4:144120554_C_A is on the other."""

    def setUp(self):
        self.ref = Allele(
            genotype="GYPA*02",
            phenotype=".",
            genotype_alt=".",
            phenotype_alt=".",
            defining_variants=frozenset({"4:144120554_ref", "4:144120555_ref", "4:144120567_ref"}),
            null=False,
            weight_geno=1000,
            reference=True,
            sub_type="GYPA*02",
        )
        self.other = Allele(
            genotype="GYPA*08",
            phenotype=".",
            genotype_alt=".",
            phenotype_alt=".",
            defining_variants=frozenset({"4:144120554_ref", "4:144120555_ref", "4:144120567_A_G"}),
            null=False,
            weight_geno=1000,
            reference=False,
            sub_type="GYPA*08",
        )
        self.pair = Pair(allele1=self.ref, allele2=self.other)

        self.bg = BloodGroup(
            type="GYPA",
            alleles={AlleleState.NORMAL: [self.pair]},
            sample="HG00436",
            variant_pool={
                "4:144120554_ref": Zygosity.HET,
                "4:144120555_ref": Zygosity.HOM,
                "4:144120567_A_G": Zygosity.HET,
                "4:144120567_ref": Zygosity.HET,
            },
            filtered_out=defaultdict(list),
        )
        cant_pair_with_ref_cuz_shared_variant_has_too_few_copies({1: self.bg})

    def test_pair_removed(self):
        self.assertEqual(self.bg.alleles[AlleleState.NORMAL], [])
        self.assertTrue(
            self.pair
            in self.bg.filtered_out[
                "cant_pair_with_ref_cuz_shared_variant_has_too_few_copies"
            ]
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


class TestCantHave2NonRefAllelesCuzOnly1GeneCopy(unittest.TestCase):
    """One gene copy on two chromosomes leaves one place for a real allele.

    The pair naming two different ones is claiming both are on it. Until v2.4.6 this was
    warned about at the reporting layer and then written out as called - the one
    impossible result the pipeline knowingly reported.
    """

    FILTER = "cant_have_2_non_ref_alleles_cuz_only_1_gene_copy"

    def _allele(self, genotype, reference=False, variants=("1:1_A_G",)):
        return Allele(
            genotype=genotype,
            phenotype=".",
            genotype_alt=".",
            phenotype_alt=".",
            defining_variants=frozenset(variants),
            null=False,
            weight_geno=1000,
            reference=reference,
            sub_type="RHD*01",
        )

    def _bg(self, pairs, locus_copies=1, chrom_copies=2, co=None):
        alleles = {AlleleState.NORMAL: pairs, AlleleState.CO: co}
        return BloodGroup(
            type="RHD",
            alleles=alleles,
            sample="s1",
            locus_copies=locus_copies,
            chrom_copies=chrom_copies,
            filtered_out=defaultdict(list),
        )

    def setUp(self):
        self.ref = self._allele("RHD*01", reference=True, variants=("1:1_ref",))
        self.alt1 = self._allele("RHD*09.01", variants=("1:1_A_G",))
        self.alt2 = self._allele("RHD*15", variants=("1:2_C_T",))
        self.two_non_ref = Pair(allele1=self.alt1, allele2=self.alt2)
        self.one_non_ref = Pair(allele1=self.ref, allele2=self.alt1)
        self.duplicate = Pair(allele1=self.alt1, allele2=self.alt1)

    def test_two_distinct_non_reference_alleles_are_excluded(self):
        bg = self._bg([self.two_non_ref, self.one_non_ref])

        cant_have_2_non_ref_alleles_cuz_only_1_gene_copy({1: bg})

        self.assertNotIn(self.two_non_ref, bg.alleles[AlleleState.NORMAL])

    def test_the_exclusion_is_recorded_under_the_filter_name(self):
        """The whole point - a warning is not an audit trail."""
        bg = self._bg([self.two_non_ref, self.one_non_ref])

        cant_have_2_non_ref_alleles_cuz_only_1_gene_copy({1: bg})

        self.assertIn(self.two_non_ref, bg.filtered_out[self.FILTER])

    def test_a_pair_with_the_reference_is_kept(self):
        """One real allele plus the slot the missing copy displaces. Ordinary."""
        bg = self._bg([self.two_non_ref, self.one_non_ref])

        cant_have_2_non_ref_alleles_cuz_only_1_gene_copy({1: bg})

        self.assertIn(self.one_non_ref, bg.alleles[AlleleState.NORMAL])

    def test_the_same_allele_twice_is_kept(self):
        """The duplicate the pairing machinery writes. One allele, one copy."""
        bg = self._bg([self.duplicate])

        cant_have_2_non_ref_alleles_cuz_only_1_gene_copy({1: bg})

        self.assertIn(self.duplicate, bg.alleles[AlleleState.NORMAL])

    def test_two_gene_copies_are_untouched(self):
        """The ordinary sample. Two alleles, two places to put them."""
        bg = self._bg([self.two_non_ref], locus_copies=None)

        cant_have_2_non_ref_alleles_cuz_only_1_gene_copy({1: bg})

        self.assertIn(self.two_non_ref, bg.alleles[AlleleState.NORMAL])
        self.assertEqual(bg.filtered_out[self.FILTER], [])

    def test_one_chromosome_is_untouched(self):
        """XK in a male is one slot, not two, and get_genotypes renders it with '-'."""
        bg = self._bg([self.two_non_ref], locus_copies=1, chrom_copies=1)

        cant_have_2_non_ref_alleles_cuz_only_1_gene_copy({1: bg})

        self.assertIn(self.two_non_ref, bg.alleles[AlleleState.NORMAL])

    def test_co_existing_alleles_are_left_alone(self):
        """Two alleles on one chromosome is what co-existing means."""
        bg = self._bg([self.two_non_ref], co=[self.two_non_ref])

        cant_have_2_non_ref_alleles_cuz_only_1_gene_copy({1: bg})

        self.assertIn(self.two_non_ref, bg.alleles[AlleleState.NORMAL])
        self.assertEqual(bg.filtered_out[self.FILTER], [])

    def test_an_input_with_no_copy_number_channel_is_a_no_op(self):
        """locus_copies is None on every input that does not encode it."""
        bg = self._bg([self.two_non_ref], locus_copies=None)

        cant_have_2_non_ref_alleles_cuz_only_1_gene_copy({1: bg})

        self.assertEqual(bg.filtered_out[self.FILTER], [])

    def test_an_empty_blood_group_does_not_raise(self):
        bg = self._bg([])

        cant_have_2_non_ref_alleles_cuz_only_1_gene_copy({1: bg})

        self.assertEqual(bg.alleles[AlleleState.NORMAL], [])


class TestCantNameSecondSlotCuzRefImpossible(unittest.TestCase):
    """Rule 4 - one slot named, the other refused.

    Built on the real GYPA case. GYPA is a lane blood group, so the genome reference is
    GYPA*02 and its defining variants are all '_ref'. A sample heterozygous at 144120554,
    with the caller doubting 144120555 so GYPA*01 is dropped, and homozygous alternate at
    144120567, leaves GYPA*08 supported and GYPA*02 impossible - there is no reference
    copy at 144120567 for it to sit on.
    """

    FILTER = "cant_name_second_slot_cuz_ref_impossible"

    REF_554 = "4:144120554_ref"
    REF_555 = "4:144120555_ref"
    REF_567 = "4:144120567_ref"
    ALT_567 = "4:144120567_A_G"

    def _allele(self, genotype, variants, reference=False, sub_type="GYPA*01"):
        return Allele(
            genotype=genotype,
            phenotype=".",
            genotype_alt=".",
            phenotype_alt=".",
            defining_variants=frozenset(variants),
            null=False,
            weight_geno=1000,
            reference=reference,
            sub_type=sub_type,
        )

    def _bg(self, pairs, pool=None, co=None):
        return BloodGroup(
            type="GYPA",
            alleles={AlleleState.NORMAL: pairs, AlleleState.CO: co},
            sample="s1",
            variant_pool=self.pool if pool is None else pool,
            filtered_out=defaultdict(list),
        )

    def setUp(self):
        self.ref = self._allele(
            "GYPA*02",
            (self.REF_554, self.REF_555, self.REF_567),
            reference=True,
            sub_type="GYPA*02",
        )
        self.mc = self._allele(
            "GYPA*08",
            (self.REF_554, self.REF_555, self.ALT_567),
            sub_type="GYPA*08",
        )
        self.other = self._allele("GYPA*11", ("4:144120558_G_T",), sub_type="GYPA*11")
        # 144120567_ref is missing because the alternate there is homozygous.
        self.pool = {
            self.REF_554: Zygosity.HET,
            self.REF_555: Zygosity.HOM,
            self.ALT_567: Zygosity.HOM,
        }
        self.impossible = Pair(allele1=self.ref, allele2=self.mc)

    def test_the_supported_half_is_named_and_the_other_refused(self):
        bg = self._bg([self.impossible])

        cant_name_second_slot_cuz_ref_impossible({1: bg})

        self.assertEqual(
            bg.single_slot_genotypes, [f"GYPA*08/{UNDETERMINED_SLOT}"]
        )

    def test_the_pair_it_came_from_is_removed(self):
        """Reporting GYPA*02 asserts wildtype on a chromosome there is evidence against."""
        bg = self._bg([self.impossible])

        cant_name_second_slot_cuz_ref_impossible({1: bg})

        self.assertEqual(bg.alleles[AlleleState.NORMAL], [])

    def test_the_exclusion_is_recorded_under_the_filter_name(self):
        """Hard rule 3 - a result is not an audit trail."""
        bg = self._bg([self.impossible])

        cant_name_second_slot_cuz_ref_impossible({1: bg})

        self.assertIn(self.impossible, bg.filtered_out[self.FILTER])

    def test_get_genotypes_writes_the_named_slot_out(self):
        """The half call has to survive to the TSV or none of this is visible."""
        bg = self._bg([self.impossible])

        cant_name_second_slot_cuz_ref_impossible({1: bg})
        out = get_genotypes({"GYPA": bg})["GYPA"]

        self.assertEqual(out.genotypes, [f"GYPA*08/{UNDETERMINED_SLOT}"])

    def test_a_supported_reference_is_untouched(self):
        """The ordinary case - the pool holds every '_ref' the reference needs."""
        pool = dict(self.pool)
        pool[self.REF_567] = Zygosity.HET
        bg = self._bg([self.impossible], pool=pool)

        cant_name_second_slot_cuz_ref_impossible({1: bg})

        self.assertEqual(bg.alleles[AlleleState.NORMAL], [self.impossible])
        self.assertEqual(bg.single_slot_genotypes, [])

    def test_one_callable_pair_stops_it(self):
        """no_defining_variant's case. A blood group with an answer keeps it."""
        callable_pair = Pair(allele1=self.mc, allele2=self.other)
        bg = self._bg([self.impossible, callable_pair])

        cant_name_second_slot_cuz_ref_impossible({1: bg})

        self.assertEqual(
            bg.alleles[AlleleState.NORMAL], [self.impossible, callable_pair]
        )
        self.assertEqual(bg.single_slot_genotypes, [])

    def test_two_reference_alleles_are_left_alone(self):
        """The reference fallback of a blood group whose alleles were all filtered."""
        bg = self._bg([Pair(allele1=self.ref, allele2=self.ref)])

        cant_name_second_slot_cuz_ref_impossible({1: bg})

        self.assertEqual(bg.single_slot_genotypes, [])

    def test_an_empty_pool_is_absence_of_evidence(self):
        """Nothing was measured, so nothing is contradicted. Same reading as
        no_defining_variant."""
        bg = self._bg([self.impossible], pool={})

        cant_name_second_slot_cuz_ref_impossible({1: bg})

        self.assertEqual(bg.alleles[AlleleState.NORMAL], [self.impossible])
        self.assertEqual(bg.single_slot_genotypes, [])

    def test_a_co_existing_result_is_left_alone(self):
        bg = self._bg([self.impossible], co=[self.impossible])

        cant_name_second_slot_cuz_ref_impossible({1: bg})

        self.assertEqual(bg.single_slot_genotypes, [])

    def test_an_empty_blood_group_does_not_raise(self):
        bg = self._bg([])

        cant_name_second_slot_cuz_ref_impossible({1: bg})

        self.assertEqual(bg.single_slot_genotypes, [])

    def test_duplicate_partners_are_named_once(self):
        """Two impossible pairs naming the same allele is one answer, not two."""
        twin = self._allele(
            "GYPA*02",
            (self.REF_554, self.REF_567),
            reference=True,
            sub_type="GYPA*02",
        )
        bg = self._bg([self.impossible, Pair(allele1=twin, allele2=self.mc)])

        cant_name_second_slot_cuz_ref_impossible({1: bg})

        self.assertEqual(
            bg.single_slot_genotypes, [f"GYPA*08/{UNDETERMINED_SLOT}"]
        )


class TestCantNameSecondSlotCuzSharedVariantHasTooFewCopies(unittest.TestCase):
    """Both alleles possible on their own, the pair impossible, so both are candidates.

    Built on the real HG00436 GYPA case, long read. The caller emitted no row at
    144120555, so GYPA*01 was never built and the only pair left is GYPA*02/GYPA*08 -
    which both need 144120554_ref, and the sample has one copy of it. One of them is on
    that chromosome and the other chromosome carries 144120554_C_A, which completes
    nothing. Without phase there is no saying which, so both are reported.
    """

    FILTER = "cant_name_second_slot_cuz_shared_variant_has_too_few_copies"

    REF_554 = "4:144120554_ref"
    REF_555 = "4:144120555_ref"
    REF_567 = "4:144120567_ref"
    ALT_567 = "4:144120567_A_G"

    def _allele(self, genotype, variants, reference=False, sub_type="GYPA*01"):
        return Allele(
            genotype=genotype,
            phenotype=".",
            genotype_alt=".",
            phenotype_alt=".",
            defining_variants=frozenset(variants),
            null=False,
            weight_geno=1000,
            reference=reference,
            sub_type=sub_type,
        )

    def _bg(self, pairs, pool=None, co=None):
        return BloodGroup(
            type="GYPA",
            alleles={AlleleState.NORMAL: pairs, AlleleState.CO: co},
            sample="HG00436",
            variant_pool=self.pool if pool is None else pool,
            filtered_out=defaultdict(list),
        )

    def setUp(self):
        self.ref = self._allele(
            "GYPA*02",
            (self.REF_554, self.REF_555, self.REF_567),
            reference=True,
            sub_type="GYPA*02",
        )
        self.mc = self._allele(
            "GYPA*08",
            (self.REF_554, self.REF_555, self.ALT_567),
            sub_type="GYPA*08",
        )
        self.elsewhere = self._allele(
            "GYPA*11", ("4:144120558_G_T",), sub_type="GYPA*11"
        )
        self.pool = {
            self.REF_554: Zygosity.HET,
            self.REF_555: Zygosity.HOM,
            self.ALT_567: Zygosity.HET,
            self.REF_567: Zygosity.HET,
        }
        self.impossible = Pair(allele1=self.ref, allele2=self.mc)

    def test_both_candidates_are_named(self):
        """Naming one would be a guess; naming neither throws away what is known."""
        bg = self._bg([self.impossible])

        cant_name_second_slot_cuz_shared_variant_has_too_few_copies({1: bg})

        self.assertEqual(
            bg.single_slot_genotypes,
            [f"GYPA*02/{UNDETERMINED_SLOT}", f"GYPA*08/{UNDETERMINED_SLOT}"],
        )

    def test_the_pair_it_came_from_is_removed(self):
        bg = self._bg([self.impossible])

        cant_name_second_slot_cuz_shared_variant_has_too_few_copies({1: bg})

        self.assertEqual(bg.alleles[AlleleState.NORMAL], [])

    def test_the_exclusion_is_recorded_under_the_filter_name(self):
        """Hard rule 3 - a result is not an audit trail."""
        bg = self._bg([self.impossible])

        cant_name_second_slot_cuz_shared_variant_has_too_few_copies({1: bg})

        self.assertIn(self.impossible, bg.filtered_out[self.FILTER])

    def test_get_genotypes_writes_both_named_slots_out(self):
        bg = self._bg([self.impossible])

        cant_name_second_slot_cuz_shared_variant_has_too_few_copies({1: bg})
        out = get_genotypes({"GYPA": bg})["GYPA"]

        self.assertEqual(
            out.genotypes,
            [f"GYPA*02/{UNDETERMINED_SLOT}", f"GYPA*08/{UNDETERMINED_SLOT}"],
        )

    def test_a_callable_pair_stops_it(self):
        """The blood group has an answer, so cant_pair_with_ref_cuz_shared_variant_
        has_too_few_copies removes the impossible pair and this does nothing."""
        pool = dict(self.pool)
        pool["4:144120558_G_T"] = Zygosity.HET
        callable_pair = Pair(allele1=self.ref, allele2=self.elsewhere)
        bg = self._bg([self.impossible, callable_pair], pool=pool)

        cant_name_second_slot_cuz_shared_variant_has_too_few_copies({1: bg})

        self.assertEqual(
            bg.alleles[AlleleState.NORMAL], [self.impossible, callable_pair]
        )
        self.assertEqual(bg.single_slot_genotypes, [])

    def test_a_pool_that_supplies_both_is_untouched(self):
        """The ordinary case - two copies of the shared variant, so no contradiction."""
        pool = dict(self.pool)
        pool[self.REF_554] = Zygosity.HOM
        bg = self._bg([self.impossible], pool=pool)

        cant_name_second_slot_cuz_shared_variant_has_too_few_copies({1: bg})

        self.assertEqual(bg.alleles[AlleleState.NORMAL], [self.impossible])
        self.assertEqual(bg.single_slot_genotypes, [])

    def test_a_reference_the_pool_contradicts_is_not_offered(self):
        """A reference needing a variant that is absent entirely cannot be a candidate,
        which is ref_slot_is_impossible's claim rather than this filter's."""
        pool = {
            self.REF_554: Zygosity.HET,
            self.REF_555: Zygosity.HOM,
            self.ALT_567: Zygosity.HOM,
        }
        bg = self._bg([self.impossible], pool=pool)

        cant_name_second_slot_cuz_shared_variant_has_too_few_copies({1: bg})

        self.assertEqual(
            bg.single_slot_genotypes, [f"GYPA*08/{UNDETERMINED_SLOT}"]
        )

    def test_a_co_existing_result_is_not_overridden(self):
        bg = self._bg([self.impossible], co=[self.impossible])

        cant_name_second_slot_cuz_shared_variant_has_too_few_copies({1: bg})

        self.assertEqual(bg.alleles[AlleleState.NORMAL], [self.impossible])
        self.assertEqual(bg.single_slot_genotypes, [])


class TestRefSlotIsImpossible(unittest.TestCase):
    """The two exemptions, which are no_defining_variant's and are here for its reasons."""

    def _allele(self, variants, genotype="X*01"):
        return Allele(
            genotype=genotype,
            phenotype=".",
            genotype_alt=".",
            phenotype_alt=".",
            defining_variants=frozenset(variants),
            null=False,
            weight_geno=1000,
            reference=True,
            sub_type=genotype,
        )

    def test_a_missing_ref_token_is_impossible(self):
        allele = self._allele(("4:144120567_ref",))

        self.assertTrue(ref_slot_is_impossible(allele, {"4:144120567_A_G": Zygosity.HOM}))

    def test_the_ABO_delG_locus_is_exempt(self):
        """The database treats the deletion as reference, so ABO*A1.01 is a reference
        allele defined by an alternate and its absence is the ordinary group O state."""
        allele = self._allele(("9:133257521_T_TC",), genotype="ABO*A1.01")

        self.assertFalse(ref_slot_is_impossible(allele, {"9:133257521_ref": Zygosity.HOM}))

    def test_an_allele_defined_only_by_absence_markers_is_exempt(self):
        """It makes no claim about a locus, so there is nothing to contradict."""
        allele = self._allele((".",))

        self.assertFalse(ref_slot_is_impossible(allele, {"4:1_A_G": Zygosity.HOM}))


if __name__ == "__main__":
    unittest.main()
