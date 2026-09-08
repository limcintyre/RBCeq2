import unittest
from collections import defaultdict
from unittest.mock import MagicMock, patch

# --- Import the actual components from your project ---
from rbceq2.core_logic.alleles import Allele, BloodGroup, Pair
from rbceq2.core_logic.constants import AlleleState, UNDETERMINED_SLOT
from rbceq2.core_logic.utils import Zygosity

# Import the functions to be tested
from rbceq2.filters.phased import (
    locus_has_a_copy_number,
    no_defining_variant,
    _get_allele_phase_info,
    cant_be_hom_ref_due_to_HET_SNP,
    cant_name_second_slot_cuz_hom_ref_impossible,
    cant_name_second_slot_cuz_ref_not_phased,
    check_phase,
    ref_not_phased,
    filter_if_all_HET_vars_on_same_side_and_phased,
    filter_on_in_relationship_if_HET_vars_on_dif_side_and_phased,
    filter_pairs_by_phase,
    impossible_alleles_phased,
    narrow_second_slot_candidates_by_phase,
    iterate_over_list,
    low_weight_hom,
    remove_unphased,
    rm_ref_if_2x_HET_phased,
)


# --- Helper Mock Class for Testing ---
class MockAllele(Allele):
    def __init__(self, **kwargs):
        defaults = {
            "phenotype": ".",
            "genotype_alt": ".",
            "phenotype_alt": ".",
            "defining_variants": frozenset(),
            "null": False,
            "weight_geno": 1000,
            "reference": False,
            "sub_type": "default",
        }
        
        defaults.update(kwargs)
        
        # Ensure defining_variants is a frozenset for set operations
        if "defining_variants" in defaults and not isinstance(defaults["defining_variants"], frozenset):
             defaults["defining_variants"] = frozenset(defaults["defining_variants"])
             
        # Call parent init (Allele is frozen, so we can't set attrs normally after this)
        super().__init__(**defaults)
        
        # Manually set .alleles to [self] using object.__setattr__ 
        # This mimics a Pair-like interface if utilities try to flatten/iterate it
        object.__setattr__(self, "alleles", [self])

    def __eq__(self, other):
        return isinstance(other, Allele) and self.genotype == other.genotype

    def __hash__(self):
        return hash(self.genotype)
    
    def __contains__(self, other):
        """
        Allows the check: if other in self:
        Returns True if 'other' variants are a subset of 'self' variants.
        """
        if not isinstance(other, Allele):
            return False
        # CRITICAL FIX: The real Allele class returns False if comparing to self
        # to prevent an allele from filtering itself out.
        if self == other:
            return False
            
        return other.defining_variants.issubset(self.defining_variants)


# --- Base Test Class ---
class TestPhasedFilters(unittest.TestCase):
    def setUp(self):
        """Set up a fresh mock BloodGroup for each test."""
        self.mock_bg = MagicMock(spec=BloodGroup)
        self.mock_bg.alleles = {
            AlleleState.FILT: [],
            AlleleState.NORMAL: [],
            AlleleState.CO: [],
        }
        self.mock_bg.filtered_out = defaultdict(list)
        self.mock_bg.variant_pool = {}
        self.mock_bg.variant_pool_phase = {}
        self.mock_bg.variant_pool_phase_set = {}
        self.mock_bg.type = "Undefined"

        def mock_remove(items_to_remove, filter_name, state=AlleleState.NORMAL):
            self.mock_bg.filtered_out[filter_name].extend(items_to_remove)
            current_items = self.mock_bg.alleles.get(state)
            if not current_items:  # Prevent error if list is already empty
                return

            is_allele_list = isinstance(current_items[0], Allele)
            if is_allele_list:
                self.mock_bg.alleles[state] = [
                    a for a in current_items if a not in items_to_remove
                ]
            else:
                self.mock_bg.alleles[state] = [
                    p for p in current_items if p not in items_to_remove
                ]

        self.mock_bg.remove_alleles.side_effect = mock_remove
        self.mock_bg.remove_pairs.side_effect = mock_remove


# --- Tests for Helper Functions ---
class TestPhasedHelperFunctions(unittest.TestCase):
    def test_get_allele_phase_info(self):
        allele = MockAllele(genotype="A1", defining_variants={"var1", "var2"})
        phase_dict = {"var1": "1|0", "var2": "0|1", "var3": "1/1"}
        result = _get_allele_phase_info(allele, phase_dict)
        self.assertCountEqual(result, ["1|0", "0|1"])

    def test_check_phase(self):
        allele = MockAllele(genotype="A1", defining_variants={"het1", "het2", "hom1"})
        variant_pool_true = {"het1": "setA", "het2": "setA", "hom1": "."}
        self.assertTrue(check_phase(variant_pool_true, allele, "."))
        variant_pool_false = {"het1": "setA", "het2": "setB", "hom1": "."}
        self.assertFalse(check_phase(variant_pool_false, allele, "."))

    def test_iterate_over_list(self):
        a1 = MockAllele(genotype="A", defining_variants={"v1"})
        a2 = MockAllele(genotype="B", defining_variants={"v1", "v2"})
        a3 = MockAllele(genotype="C", defining_variants={"v3"})
        self.assertCountEqual(iterate_over_list([a1, a2, a3]), [a1])


# --- Tests for Main Filter Functions ---


class TestRemoveUnphased(TestPhasedFilters):
    @patch("rbceq2.filters.phased.identify_unphased")
    def test_removes_unphased_alleles_when_phased(self, mock_identify_unphased):
        allele_to_remove = MockAllele(genotype="unphased")
        mock_identify_unphased.return_value = [allele_to_remove]
        self.mock_bg.alleles[AlleleState.FILT] = [allele_to_remove]
        remove_unphased({1: self.mock_bg}, phased=True)
        mock_identify_unphased.assert_called_once()
        self.mock_bg.remove_alleles.assert_called_once_with(
            [allele_to_remove], "remove_unphased", AlleleState.FILT
        )


class TestFilterIfAllHetVarsOnSameSide(TestPhasedFilters):
    def test_removes_pair_if_het_vars_share_phase(self):
        het1 = "1:100_A_G"
        het2 = "1:200_C_T"
        hom1 = "1:300_G_A"
        a1 = MockAllele(genotype="A1", defining_variants={het1, hom1})
        a2 = MockAllele(genotype="A2", defining_variants={het2})
        pair_to_remove = Pair(a1, a2)
        self.mock_bg.alleles[AlleleState.NORMAL] = [pair_to_remove]
        self.mock_bg.variant_pool = {
            het1: Zygosity.HET,
            hom1: Zygosity.HOM,
            het2: Zygosity.HET,
        }
        self.mock_bg.variant_pool_phase = {het1: "0|1", hom1: "1/1", het2: "0|1"}
        self.mock_bg.variant_pool_phase_set = {het1: "1000", hom1: ".", het2: "1000"}
        filter_if_all_HET_vars_on_same_side_and_phased({1: self.mock_bg}, phased=True)
        self.assertIn(
            pair_to_remove,
            self.mock_bg.filtered_out["filter_if_all_HET_vars_on_same_side_and_phased"],
        )


class TestFilterOnInRelationshipIfHET(TestPhasedFilters):
    def test_removes_hom_subset_pair_when_hets_are_opposite(self):
        """
        Mixed Hom/Het alleles are now skipped by find_phase (returns len=2),
        so the filter does NOT remove this pair.
        """
        hom_allele = MockAllele(
            genotype="LU*02", defining_variants={"hom1"}, reference=True
        )
        het_allele1 = MockAllele(
            genotype="LU*02.-13", defining_variants={"hom1", "het1", "het2"}
        )
        het_allele2 = MockAllele(
            genotype="LU*02.19", defining_variants={"hom1", "het3"}
        )
        pair_to_remove = Pair(hom_allele, het_allele1)
        pair_to_keep = Pair(het_allele1, het_allele2)
        self.mock_bg.alleles[AlleleState.NORMAL] = [pair_to_remove, pair_to_keep]
        self.mock_bg.variant_pool = {
            "hom1": Zygosity.HOM,
            "het1": Zygosity.HET,
            "het2": Zygosity.HET,
            "het3": Zygosity.HET,
        }
        self.mock_bg.variant_pool_phase = {
            "hom1": "1/1",
            "het1": "1|0",
            "het2": "1|0",
            "het3": "0|1",
        }
        self.mock_bg.variant_pool_phase_set = {
            "hom1": "1",
            "het1": "",
            "het2": "1",
            "het3": "1",
        }
        filter_on_in_relationship_if_HET_vars_on_dif_side_and_phased(
            {1: self.mock_bg}, phased=True
        )
        
        self.assertNotIn(
            pair_to_remove,
            self.mock_bg.filtered_out[
                "filter_on_in_relationship_if_HET_vars_on_dif_side_and_phased"
            ],
        )


class TestFilterPairsByPhase(TestPhasedFilters):
    def test_removes_pair_with_same_phase_set(self):
        # CORRECTED TEST
        het1 = "1:100_A_G"
        het2 = "1:200_C_T"
        het3 = "1:300_G_A"
        a1_same_phase = MockAllele(genotype="A1", defining_variants={het1})
        a2_same_phase = MockAllele(genotype="A2", defining_variants={het2})
        a3_diff_phase = MockAllele(genotype="A3", defining_variants={het3})

        pair_to_remove = Pair(a1_same_phase, a2_same_phase)
        # Add a valid pair to ensure the "replace all with ref" logic is not triggered
        pair_to_keep = Pair(a1_same_phase, a3_diff_phase)

        self.mock_bg.alleles[AlleleState.NORMAL] = [pair_to_remove, pair_to_keep]
        self.mock_bg.variant_pool = {
            het1: Zygosity.HET,
            het2: Zygosity.HET,
            het3: Zygosity.HET,
        }
        self.mock_bg.variant_pool_phase = {het1: "0|1", het2: "0|1", het3: "1|0"}
        self.mock_bg.variant_pool_phase_set = {
            het1: "setA",
            het2: "setA",
            het3: "setA",
        }

        filter_pairs_by_phase({1: self.mock_bg}, phased=True, reference_alleles={})

        self.assertIn(
            pair_to_remove, self.mock_bg.filtered_out["filter_pairs_by_phase"]
        )
        self.assertNotIn(pair_to_remove, self.mock_bg.alleles[AlleleState.NORMAL])
        self.assertIn(pair_to_keep, self.mock_bg.alleles[AlleleState.NORMAL])

    def test_replaces_with_ref_if_all_pairs_removed(self):
        self.mock_bg.type = "FUT2"
        ref_allele = MockAllele(genotype="FUT2*REF", reference=True)
        v1 = "19:48703728_G_A"
        v2 = "19:48703417_G_A"
        a1 = MockAllele(genotype="FUT2*01N.16", defining_variants={v1})
        a2 = MockAllele(genotype="FUT2*01N.02", defining_variants={v2})
        pair_to_remove = Pair(a1, a2)

        self.mock_bg.alleles[AlleleState.NORMAL] = [pair_to_remove]
        self.mock_bg.variant_pool = {v1: Zygosity.HET, v2: Zygosity.HET}
        self.mock_bg.variant_pool_phase = {v1: "0|1", v2: "0|1"}
        self.mock_bg.variant_pool_phase_set = {v1: "setA", v2: "setA"}

        filter_pairs_by_phase(
            {1: self.mock_bg}, phased=True, reference_alleles={"FUT2": ref_allele}
        )

        self.assertIn(
            pair_to_remove, self.mock_bg.filtered_out["filter_pairs_by_phase"]
        )
        self.assertCountEqual(
            self.mock_bg.alleles[AlleleState.NORMAL],
            [Pair(ref_allele, a1), Pair(ref_allele, a2)],
        )


class TestImpossibleAllelesPhased(TestPhasedFilters):
    def test_removes_phased_subset_allele_pair(self):
        allele_subset = MockAllele(
            genotype="GYPB*03N.03", defining_variants={"het1", "hom1"}
        )
        allele_superset = MockAllele(
            genotype="GYPB*03N.04", defining_variants={"het1", "hom1", "het2", "het3"}
        )
        other_allele = MockAllele(genotype="Other")
        pair_to_remove = Pair(allele_subset, other_allele)
        pair_to_keep = Pair(allele_superset, other_allele)
        self.mock_bg.alleles[AlleleState.NORMAL] = [pair_to_remove, pair_to_keep]
        self.mock_bg.variant_pool = {
            "het1": Zygosity.HET,
            "hom1": Zygosity.HOM,
            "het2": Zygosity.HET,
            "het3": Zygosity.HET,
        }
        self.mock_bg.variant_pool_phase = {
            "het1": "0|1",
            "hom1": "1/1",
            "het2": "0|1",
            "het3": "0|1",
        }
        self.mock_bg.variant_pool_phase_set = {
            "het1": "setA",
            "hom1": ".",
            "het2": "setA",
            "het3": "setA",
        }
        impossible_alleles_phased({1: self.mock_bg}, phased=True)
        self.assertIn(
            pair_to_remove,
            self.mock_bg.filtered_out["filter_impossible_alleles_phased"],
        )
        self.assertNotIn(pair_to_remove, self.mock_bg.alleles[AlleleState.NORMAL])

class TestNoDefiningVariantEmptyPool(TestPhasedFilters):
    """An empty variant pool is absence of evidence, not evidence of impossibility.

    The filter exists to drop a reference allele the data contradicts - '_ref' at a locus
    where the alternate is homozygous. When a filter has emptied the blood group the pool
    is empty too, so every defining variant is 'not in the pool' and the reference pair is
    removed, turning the conventional fall back to the reference allele into a no call.
    HG02308 KN and HG03673 RHD both reported no call where the answer was the reference.
    """

    def _ref_pair(self):
        ref = MockAllele(
            genotype="KN*01",
            defining_variants={"1:207609424_ref", "1:207609571_A_T"},
            reference=True,
        )
        return Pair(ref, ref)

    def test_empty_pool_keeps_the_reference_pair(self):
        self.mock_bg.alleles[AlleleState.NORMAL] = [self._ref_pair()]
        self.mock_bg.variant_pool = {}
        no_defining_variant({1: self.mock_bg})
        self.mock_bg.remove_pairs.assert_not_called()
        self.assertEqual(len(self.mock_bg.alleles[AlleleState.NORMAL]), 1)

    def test_populated_pool_still_removes_a_contradicted_reference(self):
        """The behaviour the filter is for is unchanged."""
        pair = self._ref_pair()
        self.mock_bg.alleles[AlleleState.NORMAL] = [pair]
        # the locus was seen, and the reference token is not what was found there.
        # 207609424, not 207609571 - the contradicting call has to be at the locus the
        # missing token names, which is what locus_has_a_copy_number asks.
        self.mock_bg.variant_pool = {"1:207609424_C_A": Zygosity.HOM}
        no_defining_variant({1: self.mock_bg})
        self.mock_bg.remove_pairs.assert_called_once_with([pair], "no_defining_variant")

    def test_it_runs_unphased_too(self):
        """The flag used to decide this, and it is about something else. Nothing in the
        filter reads phase, so the same pool has to give the same answer in both arms."""
        pair = self._ref_pair()
        self.mock_bg.alleles[AlleleState.NORMAL] = [pair]
        self.mock_bg.variant_pool = {"1:207609424_C_A": Zygosity.HOM}
        no_defining_variant({1: self.mock_bg})
        self.mock_bg.remove_pairs.assert_called_once_with([pair], "no_defining_variant")

    def test_a_locus_nobody_called_is_not_evidence(self):
        """HG00109 HPA3 on an array. The probe did not call the locus, so the reference
        token is missing for want of a measurement rather than because of one."""
        pair = self._ref_pair()
        self.mock_bg.alleles[AlleleState.NORMAL] = [pair]
        self.mock_bg.variant_pool = {"1:207609424_C_A": Zygosity.NO_DATA}
        no_defining_variant({1: self.mock_bg})
        self.mock_bg.remove_pairs.assert_not_called()

    def test_a_locus_with_no_copies_is_evidence(self):
        """Zero copies is a measurement of absence, and a locus with no copies has no
        reference copy either. Unobserved on any input to hand; pinned so the choice
        cannot drift."""
        pair = self._ref_pair()
        self.mock_bg.alleles[AlleleState.NORMAL] = [pair]
        self.mock_bg.variant_pool = {"1:207609424_C_A": Zygosity.NO_COPIES}
        no_defining_variant({1: self.mock_bg})
        self.mock_bg.remove_pairs.assert_called_once_with([pair], "no_defining_variant")


class TestLocusHasACopyNumber(unittest.TestCase):
    """Which absences are evidence."""

    def test_a_sibling_token_with_a_copy_number_counts(self):
        self.assertTrue(
            locus_has_a_copy_number(
                "1:25408711_ref", {"1:25408711_G_A": Zygosity.HOM}
            )
        )

    def test_no_data_does_not_count(self):
        self.assertFalse(
            locus_has_a_copy_number(
                "17:42453065_ref", {"17:42453065_A_C": Zygosity.NO_DATA}
            )
        )

    def test_an_empty_locus_does_not_count(self):
        self.assertFalse(
            locus_has_a_copy_number(
                "1:25408711_ref", {"1:25390874_C_G": Zygosity.HOM}
            )
        )

    def test_a_longer_position_at_the_same_prefix_is_not_the_same_locus(self):
        """'1:2540871' must not match '1:25408711'."""
        self.assertFalse(
            locus_has_a_copy_number(
                "1:2540871_ref", {"1:25408711_G_A": Zygosity.HOM}
            )
        )

    def test_another_chromosome_at_the_same_position_is_not_the_same_locus(self):
        """The substring hazard: '1:159205564' is inside '11:159205564'."""
        self.assertFalse(
            locus_has_a_copy_number(
                "1:159205564_ref", {"11:159205564_G_A": Zygosity.HOM}
            )
        )

    def test_the_alternate_form_works_too(self):
        """The missing token is not always a '_ref' one - four references are defined
        partly by a real alternate at a lane locus."""
        self.assertTrue(
            locus_has_a_copy_number(
                "1:25420739_G_C", {"1:25420739_ref": Zygosity.HOM}
            )
        )

class TestNoDefiningVariantAboDelG(TestPhasedFilters):
    """The ABO c.261delG insertion is exempt, and was only exempt on one build.

    The database treats the deletion as the reference sequence, so ABO*A1.01 is a reference
    allele defined by an alternate. Its absence from the pool is the ordinary state of a
    group O sample, not a contradiction, so removing the pair over it would make ABO
    uncallable for anyone who is not group A or B. The GRCh37 token used to be written
    without its chromosome prefix and so matched nothing.
    """

    def _abo_ref_pair(self, delg_token):
        ref = MockAllele(
            genotype="ABO*A1.01", defining_variants={delg_token}, reference=True
        )
        return Pair(ref, ref)

    def _run(self, delg_token):
        pair = self._abo_ref_pair(delg_token)
        self.mock_bg.type = "ABO"
        self.mock_bg.alleles[AlleleState.NORMAL] = [pair]
        # the locus was seen, and what was found there is the deletion, ie no insertion
        self.mock_bg.variant_pool = {delg_token.replace("_T_TC", "_ref"): Zygosity.HOM}
        no_defining_variant({1: self.mock_bg})

    def test_exempt_on_grch38(self):
        self._run("9:133257521_T_TC")
        self.mock_bg.remove_pairs.assert_not_called()

    def test_exempt_on_grch37(self):
        """This is the one that used to fail - the token was missing its '9:' prefix."""
        self._run("9:136132908_T_TC")
        self.mock_bg.remove_pairs.assert_not_called()

    def test_an_ordinary_reference_variant_is_still_removed(self):
        pair = self._abo_ref_pair("1:25390874_ref")
        self.mock_bg.alleles[AlleleState.NORMAL] = [pair]
        self.mock_bg.variant_pool = {"1:25390874_C_G": Zygosity.HOM}
        no_defining_variant({1: self.mock_bg})
        self.mock_bg.remove_pairs.assert_called_once_with([pair], "no_defining_variant")


class TestCantNameSecondSlotCuzHomRefImpossible(unittest.TestCase):
    """Rule 4, the half where the reference is the slot the data settles.

    Built on the real NA18571 RHCE case. FILTER drops RHCE*02 for two LowQual defining
    variants and remove_unphased drops the other two candidates, leaving only
    RHCE*01/RHCE*01, which cant_be_hom_ref_due_to_HET_SNP then removes because
    25408711 is heterozygous. That leaves nothing, and the answer was
    'Undetermined/Undetermined' when one chromosome is plainly RHCE*01.

    The sibling of TestCantNameSecondSlotCuzRefImpossible in test_geno_filters.py, which
    covers the mirror: reference in the slot that cannot be named.
    """

    REF_874 = "1:25390874_ref"
    REF_711 = "1:25408711_ref"
    ALT_711 = "1:25408711_G_A"
    ALT_739 = "1:25420739_G_C"

    def _allele(self, genotype, variants, reference=False):
        return Allele(
            genotype=genotype,
            phenotype=".",
            genotype_alt=".",
            phenotype_alt=".",
            defining_variants=frozenset(variants),
            null=False,
            weight_geno=1000,
            reference=reference,
            sub_type="RHCE*01",
        )

    def _bg(self, phase, pairs=None, co=None):
        return BloodGroup(
            type="RHCE",
            alleles={AlleleState.NORMAL: [] if pairs is None else pairs,
                     AlleleState.CO: co},
            sample="NA18571.vcf",
            variant_pool=dict(self.pool),
            variant_pool_phase=phase,
            filtered_out=defaultdict(list),
        )

    def setUp(self):
        self.ref = self._allele(
            "RHCE*01", (self.REF_874, self.REF_711, self.ALT_739), reference=True
        )
        self.hom_ref_pair = Pair(allele1=self.ref, allele2=self.ref)
        self.pool = {
            self.REF_874: Zygosity.HOM,
            self.ALT_711: Zygosity.HET,
            self.REF_711: Zygosity.HET,
            self.ALT_739: Zygosity.HET,
        }
        # The reference's two phased defining variants are both on the left.
        self.coherent = {
            self.REF_874: "1/1",
            self.ALT_711: "0|1",
            self.REF_711: "1|0",
            self.ALT_739: "1|0",
        }

    def _emptied(self, phase):
        """A blood group cant_be_hom_ref_due_to_HET_SNP has just emptied."""
        bg = self._bg(phase, pairs=[self.hom_ref_pair])
        cant_be_hom_ref_due_to_HET_SNP({1: bg}, phased=True)
        self.assertEqual(bg.alleles[AlleleState.NORMAL], [])
        return bg

    def test_the_reference_slot_is_named_and_the_other_refused(self):
        bg = self._emptied(self.coherent)

        cant_name_second_slot_cuz_hom_ref_impossible({1: bg}, phased=True)

        self.assertEqual(
            bg.single_slot_genotypes, [f"RHCE*01/{UNDETERMINED_SLOT}"]
        )

    def test_nothing_further_is_excluded(self):
        """The pair was removed and recorded by the filter this one reads."""
        bg = self._emptied(self.coherent)

        cant_name_second_slot_cuz_hom_ref_impossible({1: bg}, phased=True)

        self.assertEqual(
            list(bg.filtered_out), ["cant_be_hom_ref_due_to_HET_SNP"]
        )

    def test_unphased_is_left_alone(self):
        """Without phase there is nothing saying which chromosome carries the reference."""
        bg = self._emptied(self.coherent)

        cant_name_second_slot_cuz_hom_ref_impossible({1: bg}, phased=False)

        self.assertEqual(bg.single_slot_genotypes, [])

    def test_defining_variants_on_opposite_sides_are_left_alone(self):
        """Neither chromosome carries the whole reference, so refusing both is honest."""
        split = dict(self.coherent)
        split[self.ALT_739] = "0|1"
        bg = self._emptied(split)

        cant_name_second_slot_cuz_hom_ref_impossible({1: bg}, phased=True)

        self.assertEqual(bg.single_slot_genotypes, [])

    def test_a_heterozygote_without_a_bar_is_left_alone(self):
        """A partly phased file says nothing about sides, so this must not guess."""
        unbarred = dict(self.coherent)
        unbarred[self.REF_711] = "0/1"
        unbarred[self.ALT_739] = "0/1"
        bg = self._emptied(unbarred)

        cant_name_second_slot_cuz_hom_ref_impossible({1: bg}, phased=True)

        self.assertEqual(bg.single_slot_genotypes, [])

    def test_a_surviving_pair_stops_it_dead(self):
        """If anything paired the blood group has an answer and this must not touch it."""
        other = self._allele("RHCE*02", (self.REF_874, self.ALT_711))
        bg = self._bg(self.coherent, pairs=[self.hom_ref_pair])
        cant_be_hom_ref_due_to_HET_SNP({1: bg}, phased=True)
        bg.alleles[AlleleState.NORMAL] = [Pair(allele1=self.ref, allele2=other)]

        cant_name_second_slot_cuz_hom_ref_impossible({1: bg}, phased=True)

        self.assertEqual(bg.single_slot_genotypes, [])

    def test_a_coexisting_result_is_not_overridden(self):
        bg = self._bg(self.coherent, pairs=[self.hom_ref_pair])
        cant_be_hom_ref_due_to_HET_SNP({1: bg}, phased=True)
        bg.alleles[AlleleState.CO] = []

        cant_name_second_slot_cuz_hom_ref_impossible({1: bg}, phased=True)

        self.assertEqual(bg.single_slot_genotypes, [])

    def test_it_does_nothing_when_that_filter_did_not_fire(self):
        """An empty blood group emptied by something else is not this filter's case."""
        bg = self._bg(self.coherent)

        cant_name_second_slot_cuz_hom_ref_impossible({1: bg}, phased=True)

        self.assertEqual(bg.single_slot_genotypes, [])

    def test_the_removal_does_not_promise_a_revert_to_reference(self):
        """The pair removed *is* the reference pair, so the default warning would lie."""
        bg = self._bg(self.coherent, pairs=[self.hom_ref_pair])
        with patch("rbceq2.core_logic.alleles.logger") as mock_logger:
            cant_be_hom_ref_due_to_HET_SNP({1: bg}, phased=True)
        mock_logger.warning.assert_not_called()


class TestCantNameSecondSlotCuzRefNotPhased(unittest.TestCase):
    """HG00128 RHCE: the reference is split across both chromosomes and its partner is not.

    Every defining variant is heterozygous and the phase pool separates them cleanly.
    RHCE*01 wants 25390874_ref from one chromosome and 25408711_ref and 25420739_G_C
    from the other, so remove_unphased drops it; process_genetic_data re-adds it anyway
    and ref_not_phased removes the RHCE*01/RHCE*03 pair that made. That leaves nothing,
    and the answer was 'Undetermined/Undetermined' when one chromosome is plainly
    RHCE*03.

    The third of the family. TestCantNameSecondSlotCuzHomRefImpossible above covers the
    hom reference pair, and TestCantNameSecondSlotCuzRefImpossible in
    test_geno_filters.py covers the reference being the slot that cannot be named.
    """

    REF_874 = "1:25390874_ref"
    ALT_874 = "1:25390874_C_G"
    REF_711 = "1:25408711_ref"
    ALT_711 = "1:25408711_G_A"
    REF_739 = "1:25420739_ref"
    ALT_739 = "1:25420739_G_C"

    def _allele(self, genotype, variants, reference=False):
        return Allele(
            genotype=genotype,
            phenotype=".",
            genotype_alt=".",
            phenotype_alt=".",
            defining_variants=frozenset(variants),
            null=False,
            weight_geno=1000,
            reference=reference,
            sub_type="RHCE*01",
        )

    def _bg(self, phase, pairs=None, co=None, unphased_ref=True):
        filtered_out = defaultdict(list)
        if unphased_ref:
            filtered_out["remove_unphased"].append(self.ref)
        return BloodGroup(
            type="RHCE",
            alleles={AlleleState.NORMAL: [] if pairs is None else pairs,
                     AlleleState.CO: co},
            sample="HG00128.vcf",
            variant_pool=dict(self.pool),
            variant_pool_phase=phase,
            filtered_out=filtered_out,
        )

    def setUp(self):
        self.ref = self._allele(
            "RHCE*01", (self.REF_874, self.REF_711, self.ALT_739), reference=True
        )
        self.partner = self._allele(
            "RHCE*03", (self.ALT_874, self.REF_711, self.ALT_739)
        )
        self.pair = Pair(allele1=self.ref, allele2=self.partner)
        self.pool = {
            self.REF_874: Zygosity.HET,
            self.ALT_874: Zygosity.HET,
            self.REF_711: Zygosity.HET,
            self.ALT_711: Zygosity.HET,
            self.REF_739: Zygosity.HET,
            self.ALT_739: Zygosity.HET,
        }
        # The reference is split; every one of the partner's variants is on the left.
        self.coherent = {
            self.REF_874: "0|1",
            self.ALT_874: "1|0",
            self.REF_711: "1|0",
            self.ALT_711: "0|1",
            self.REF_739: "0|1",
            self.ALT_739: "1|0",
        }

    def _emptied(self, phase):
        """A blood group ref_not_phased has just emptied."""
        bg = self._bg(phase, pairs=[self.pair])
        ref_not_phased({1: bg}, phased=True)
        self.assertEqual(bg.alleles[AlleleState.NORMAL], [])
        return bg

    def test_the_partner_slot_is_named_and_the_other_refused(self):
        bg = self._emptied(self.coherent)

        cant_name_second_slot_cuz_ref_not_phased({1: bg}, phased=True)

        self.assertEqual(
            bg.single_slot_genotypes, [f"RHCE*03/{UNDETERMINED_SLOT}"]
        )

    def test_nothing_further_is_excluded(self):
        """The pair was removed and recorded by the filter this one reads."""
        bg = self._emptied(self.coherent)

        cant_name_second_slot_cuz_ref_not_phased({1: bg}, phased=True)

        self.assertEqual(
            sorted(bg.filtered_out), ["ref_not_phased", "remove_unphased"]
        )

    def test_unphased_is_left_alone(self):
        """Without phase there is nothing saying which chromosome carries the partner."""
        bg = self._emptied(self.coherent)

        cant_name_second_slot_cuz_ref_not_phased({1: bg}, phased=False)

        self.assertEqual(bg.single_slot_genotypes, [])

    def test_a_partner_split_across_both_sides_is_left_alone(self):
        """Neither chromosome carries the whole partner, so refusing both is honest."""
        split = dict(self.coherent)
        split[self.ALT_874] = "0|1"
        bg = self._emptied(split)

        cant_name_second_slot_cuz_ref_not_phased({1: bg}, phased=True)

        self.assertEqual(bg.single_slot_genotypes, [])

    def test_a_heterozygote_without_a_bar_is_left_alone(self):
        """A partly phased file says nothing about sides, so this must not guess."""
        unbarred = dict(self.coherent)
        for variant in (self.ALT_874, self.REF_711, self.ALT_739):
            unbarred[variant] = "0/1"
        bg = self._emptied(unbarred)

        cant_name_second_slot_cuz_ref_not_phased({1: bg}, phased=True)

        self.assertEqual(bg.single_slot_genotypes, [])

    def test_a_partner_located_by_nothing_is_left_alone(self):
        """Every defining variant homozygous locates the allele on neither chromosome."""
        homs = dict(self.coherent)
        for variant in (self.ALT_874, self.REF_711, self.ALT_739):
            homs[variant] = "1/1"
        bg = self._emptied(homs)

        cant_name_second_slot_cuz_ref_not_phased({1: bg}, phased=True)

        self.assertEqual(bg.single_slot_genotypes, [])

    def test_a_surviving_pair_stops_it_dead(self):
        """If anything paired the blood group has an answer and this must not touch it."""
        other = self._allele("RHCE*02", (self.REF_874, self.ALT_711))
        bg = self._bg(self.coherent, pairs=[self.pair])
        ref_not_phased({1: bg}, phased=True)
        bg.alleles[AlleleState.NORMAL] = [Pair(allele1=self.partner, allele2=other)]

        cant_name_second_slot_cuz_ref_not_phased({1: bg}, phased=True)

        self.assertEqual(bg.single_slot_genotypes, [])

    def test_a_slot_already_named_is_not_overwritten(self):
        """The filters that ran earlier own their answer."""
        bg = self._emptied(self.coherent)
        bg.single_slot_genotypes = [f"RHCE*01/{UNDETERMINED_SLOT}"]

        cant_name_second_slot_cuz_ref_not_phased({1: bg}, phased=True)

        self.assertEqual(
            bg.single_slot_genotypes, [f"RHCE*01/{UNDETERMINED_SLOT}"]
        )

    def test_a_coexisting_result_is_not_overridden(self):
        bg = self._bg(self.coherent, pairs=[self.pair])
        ref_not_phased({1: bg}, phased=True)
        bg.alleles[AlleleState.CO] = []

        cant_name_second_slot_cuz_ref_not_phased({1: bg}, phased=True)

        self.assertEqual(bg.single_slot_genotypes, [])

    def test_it_does_nothing_when_that_filter_did_not_fire(self):
        """An empty blood group emptied by something else is not this filter's case."""
        bg = self._bg(self.coherent)

        cant_name_second_slot_cuz_ref_not_phased({1: bg}, phased=True)

        self.assertEqual(bg.single_slot_genotypes, [])

    def test_every_settled_candidate_is_named(self):
        """Two partners, both settled, are two candidates rather than a choice."""
        second = self._allele("RHCE*03.18", (self.ALT_874, self.REF_711))
        bg = self._bg(
            self.coherent,
            pairs=[self.pair, Pair(allele1=self.ref, allele2=second)],
        )
        ref_not_phased({1: bg}, phased=True)

        cant_name_second_slot_cuz_ref_not_phased({1: bg}, phased=True)

        self.assertEqual(
            bg.single_slot_genotypes,
            [f"RHCE*03.18/{UNDETERMINED_SLOT}", f"RHCE*03/{UNDETERMINED_SLOT}"],
        )


class TestPhaseSetBoundaries(unittest.TestCase):
    """A GT side is comparable only within a named block on one chromosome."""

    LEFT = "1:100_A_G"
    EXTRA = "1:150_G_T"
    HOM = "1:300_G_A"
    SAME_SIDE = "filter_if_all_HET_vars_on_same_side_and_phased"
    REMOVE_REF = "rm_ref_if_2x_HET_phased"
    LOW_WEIGHT = "low_weight_hom"

    @staticmethod
    def _allele(name, variants=(), *, weight=1, reference=False):
        """Make a generic definition without changing curated allele biology."""
        return Allele(
            genotype=name,
            phenotype="TEST:1",
            genotype_alt=name,
            phenotype_alt="Test+",
            defining_variants=frozenset(variants),
            null=False,
            weight_geno=weight,
            reference=reference,
            sub_type="TEST*01",
        )

    def _group(
        self, *, right_gt="1|0", right_chrom="1", first_variants=None,
        state=AlleleState.NORMAL,
    ):
        """Return two weighted nonreference pairs and three reference choices."""
        first_variants = first_variants or (self.LEFT,)
        left = self._allele("TEST*01.01", first_variants, weight=1)
        right = self._allele("TEST*01.02", (f"{right_chrom}:200_C_T",), weight=2)
        other = self._allele("TEST*01.03", (f"{right_chrom}:250_G_A",), weight=8)
        reference = self._allele("TEST*01", reference=True)
        candidates = [Pair(left, right), Pair(left, other)]
        reference_pairs = [Pair(reference, allele) for allele in (left, right, other)]
        phase = {variant: "0|1" for variant in first_variants}
        phase.update({
            variant: right_gt
            for allele in (right, other)
            for variant in allele.defining_variants
        })
        states = {AlleleState.NORMAL: [], AlleleState.CO: []}
        states[state] = candidates + reference_pairs
        bg = BloodGroup(
            # The CO stage is restricted to KN; these definitions remain synthetic.
            type="KN" if state == AlleleState.CO else "TEST",
            sample="phase_set_boundary",
            alleles=states,
            variant_pool={variant: Zygosity.HET for variant in phase},
            variant_pool_phase=phase,
            variant_pool_phase_set={variant: "1000" for variant in phase},
        )
        return bg, candidates, reference_pairs

    @staticmethod
    def _stages():
        return (
            (filter_if_all_HET_vars_on_same_side_and_phased, AlleleState.NORMAL),
            (filter_if_all_HET_vars_on_same_side_and_phased, AlleleState.CO),
            (rm_ref_if_2x_HET_phased, AlleleState.NORMAL),
            (low_weight_hom, AlleleState.NORMAL),
        )

    def _assert_unchanged(self, stage, bg, state, *, phased=True):
        initial = list(bg.alleles[state])

        result = stage({1: bg}, phased=phased)

        self.assertIs(result[1], bg)
        self.assertEqual(bg.alleles[state], initial)
        self.assertEqual(dict(bg.filtered_out), {})

    def test_different_phase_sets_preserve_same_and_opposite_gt_pairs(self):
        for stage, state in self._stages():
            for right_gt in ("0|1", "1|0"):
                with self.subTest(stage=stage.__name__, state=state, right_gt=right_gt):
                    bg, candidates, _ = self._group(right_gt=right_gt, state=state)
                    for pair in candidates:
                        for variant in pair.allele2.defining_variants:
                            bg.variant_pool_phase_set[variant] = "2000"

                    self._assert_unchanged(stage, bg, state)

    def test_reused_phase_set_label_on_different_chromosomes_preserves_pairs(self):
        for stage, state in self._stages():
            for right_gt in ("0|1", "1|0"):
                with self.subTest(stage=stage.__name__, state=state, right_gt=right_gt):
                    bg, _, _ = self._group(
                        right_gt=right_gt, right_chrom="2", state=state
                    )

                    self._assert_unchanged(stage, bg, state)

    def test_missing_or_unknown_phase_sets_supply_no_block_evidence(self):
        for stage, state in self._stages():
            for phase_set in (None, "", ".", "unknown"):
                for right_gt in ("0|1", "1|0"):
                    with self.subTest(
                        stage=stage.__name__, state=state,
                        phase_set=phase_set, right_gt=right_gt,
                    ):
                        bg, candidates, _ = self._group(right_gt=right_gt, state=state)
                        for pair in candidates:
                            for variant in pair.allele2.defining_variants:
                                if phase_set is None:
                                    bg.variant_pool_phase_set.pop(variant)
                                else:
                                    bg.variant_pool_phase_set[variant] = phase_set

                        self._assert_unchanged(stage, bg, state)

    def test_homozygous_and_hemizygous_alleles_cannot_supply_a_side(self):
        for stage, state in self._stages():
            for copy_state, gt in (
                (Zygosity.HOM, "1/1"),
                (Zygosity.HOM, "1|1"),
                (Zygosity.HEM, "1"),
            ):
                for side in ("first", "second"):
                    with self.subTest(
                        stage=stage.__name__, state=state,
                        copy_state=copy_state, gt=gt, side=side,
                    ):
                        bg, candidates, _ = self._group(state=state)
                        alleles = (
                            [candidates[0].allele1] if side == "first"
                            else [pair.allele2 for pair in candidates]
                        )
                        for allele in alleles:
                            for variant in allele.defining_variants:
                                bg.variant_pool[variant] = copy_state
                                bg.variant_pool_phase[variant] = gt

                        self._assert_unchanged(stage, bg, state)

    def test_mixed_first_allele_het_sides_do_not_prove_trans(self):
        for stage in (rm_ref_if_2x_HET_phased, low_weight_hom):
            with self.subTest(stage=stage.__name__):
                bg, _, _ = self._group(
                    first_variants=(self.LEFT, self.EXTRA)
                )
                bg.variant_pool_phase[self.EXTRA] = "1|0"

                self._assert_unchanged(stage, bg, AlleleState.NORMAL)

    def test_phased_false_preserves_every_pair_and_exclusion_trail(self):
        for stage, state in self._stages():
            for right_gt in ("0|1", "1|0"):
                with self.subTest(stage=stage.__name__, state=state, right_gt=right_gt):
                    bg, _, _ = self._group(right_gt=right_gt, state=state)

                    self._assert_unchanged(stage, bg, state, phased=False)

    def test_same_block_cis_exclusions_keep_the_filter_name_in_both_states(self):
        for state in (AlleleState.NORMAL, AlleleState.CO):
            with self.subTest(state=state):
                bg, candidates, reference_pairs = self._group(right_gt="0|1", state=state)

                filter_if_all_HET_vars_on_same_side_and_phased({1: bg}, phased=True)

                self.assertEqual(bg.alleles[state], reference_pairs)
                self.assertEqual(set(bg.filtered_out), {self.SAME_SIDE})
                self.assertCountEqual(bg.filtered_out[self.SAME_SIDE], candidates)

    def test_same_block_trans_does_not_trigger_same_side_exclusion(self):
        for state in (AlleleState.NORMAL, AlleleState.CO):
            with self.subTest(state=state):
                bg, _, _ = self._group(state=state)

                self._assert_unchanged(
                    filter_if_all_HET_vars_on_same_side_and_phased, bg, state
                )

    def test_same_block_trans_removes_reference_pairs_by_name(self):
        bg, candidates, reference_pairs = self._group()

        rm_ref_if_2x_HET_phased({1: bg}, phased=True)

        self.assertEqual(bg.alleles[AlleleState.NORMAL], candidates)
        self.assertEqual(set(bg.filtered_out), {self.REMOVE_REF})
        self.assertCountEqual(bg.filtered_out[self.REMOVE_REF], reference_pairs)

    def test_same_block_trans_retains_the_best_weighted_pair_by_name(self):
        bg, candidates, reference_pairs = self._group()

        low_weight_hom({1: bg}, phased=True)

        self.assertEqual(bg.alleles[AlleleState.NORMAL], [candidates[0]])
        self.assertEqual(set(bg.filtered_out), {self.LOW_WEIGHT})
        self.assertCountEqual(
            bg.filtered_out[self.LOW_WEIGHT], candidates[1:] + reference_pairs
        )

    def test_cis_pairs_do_not_establish_trans_for_reference_or_weight_filters(self):
        for stage in (rm_ref_if_2x_HET_phased, low_weight_hom):
            with self.subTest(stage=stage.__name__):
                bg, _, _ = self._group(right_gt="0|1")

                self._assert_unchanged(stage, bg, AlleleState.NORMAL)

    def test_hom_tokens_do_not_obscure_valid_het_phase_evidence(self):
        for hom_gt in ("1/1", "1|1"):
            for stage, state in self._stages():
                with self.subTest(hom_gt=hom_gt, stage=stage.__name__, state=state):
                    same_side = stage is filter_if_all_HET_vars_on_same_side_and_phased
                    bg, candidates, reference_pairs = self._group(
                        right_gt="0|1" if same_side else "1|0",
                        first_variants=(self.LEFT, self.HOM), state=state,
                    )
                    bg.variant_pool[self.HOM] = Zygosity.HOM
                    bg.variant_pool_phase[self.HOM] = hom_gt
                    bg.variant_pool_phase_set[self.HOM] = "."

                    stage({1: bg}, phased=True)

                    if same_side:
                        reason, kept, removed = self.SAME_SIDE, reference_pairs, candidates
                    elif stage is rm_ref_if_2x_HET_phased:
                        reason, kept, removed = self.REMOVE_REF, candidates, reference_pairs
                    else:
                        reason = self.LOW_WEIGHT
                        kept = [candidates[0]]
                        removed = candidates[1:] + reference_pairs
                    self.assertEqual(bg.alleles[state], kept)
                    self.assertEqual(set(bg.filtered_out), {reason})
                    self.assertCountEqual(bg.filtered_out[reason], removed)

    def test_pair_phase_filter_keeps_pairs_without_a_known_shared_block(self):
        scenarios = (
            ("both_dot", ".", ".", "1"),
            ("both_unknown", "unknown", "unknown", "1"),
            ("both_absent", None, None, "1"),
            ("different_blocks", "1000", "2000", "1"),
            ("different_chromosomes", "1000", "1000", "2"),
        )
        for name, first_ps, second_ps, right_chrom in scenarios:
            with self.subTest(scenario=name):
                bg, candidates, _ = self._group(
                    right_gt="0|1", right_chrom=right_chrom
                )
                pair = candidates[0]
                bg.alleles[AlleleState.NORMAL] = [pair]
                bg.variant_pool_phase_set = {
                    variant: phase_set
                    for allele, phase_set in (
                        (pair.allele1, first_ps), (pair.allele2, second_ps)
                    )
                    if phase_set is not None
                    for variant in allele.defining_variants
                }

                result = filter_pairs_by_phase(
                    {1: bg}, phased=True, reference_alleles={}
                )

                self.assertIs(result[1], bg)
                self.assertEqual(bg.alleles[AlleleState.NORMAL], [pair])
                self.assertEqual(dict(bg.filtered_out), {})

    def test_pair_phase_filter_retains_known_shared_block_cis_exclusion(self):
        bg, candidates, reference_pairs = self._group(right_gt="0|1")

        filter_pairs_by_phase({1: bg}, phased=True, reference_alleles={})

        self.assertEqual(bg.alleles[AlleleState.NORMAL], reference_pairs)
        self.assertEqual(set(bg.filtered_out), {"filter_pairs_by_phase"})
        self.assertCountEqual(bg.filtered_out["filter_pairs_by_phase"], candidates)


if __name__ == "__main__":
    unittest.main(argv=["first-arg-is-ignored"], exit=False)


class TestNarrowSecondSlotCandidatesByPhase(unittest.TestCase):
    """HG01527 RHCE: six candidates, one chromosome, one of them holds every variant.

    cant_name_second_slot_cuz_ref_impossible names one chromosome and refuses the other,
    leaving one genotype per candidate in single_slot_genotypes. Nothing narrowed that,
    because it removed the pairs to build it and the phase based 'in' filters work on
    pairs, so they never see it.
    """

    @staticmethod
    def _bg(candidates, pool, phase):
        """A BloodGroup mid-pipeline, after the second slot was refused."""
        alleles = [
            MockAllele(genotype=name, defining_variants=variants)
            for name, variants in candidates.items()
        ]
        bg = BloodGroup(
            type="RHCE",
            alleles={AlleleState.RAW: alleles, AlleleState.NORMAL: []},
            sample="test_sample",
            variant_pool=dict(pool),
            variant_pool_phase=dict(phase),
        )
        bg.single_slot_genotypes = [
            f"{name}/{UNDETERMINED_SLOT}" for name in sorted(candidates)
        ]
        return bg

    # a subset chain, every heterozygous variant on side two - the HG01527 shape
    CANDIDATES = {
        "RHCE*01.01": {"c"},
        "RHCE*01.02.01": {"a", "c"},
        "RHCE*01.20.04.02": {"a", "b", "c"},
    }
    POOL = {
        "a": Zygosity.HET,
        "b": Zygosity.HET,
        "c": Zygosity.HOM,
    }
    ONE_SIDE = {"a": "0|1", "b": "0|1", "c": "1/1"}

    def test_the_superset_is_the_only_answer_left(self):
        bg = self._bg(self.CANDIDATES, self.POOL, self.ONE_SIDE)
        narrow_second_slot_candidates_by_phase({"RHCE": bg}, phased=True)
        self.assertEqual(
            bg.single_slot_genotypes, [f"RHCE*01.20.04.02/{UNDETERMINED_SLOT}"]
        )

    def test_what_was_dropped_is_recorded(self):
        """Every exclusion carries the filter's name - it is not a silent drop."""
        bg = self._bg(self.CANDIDATES, self.POOL, self.ONE_SIDE)
        narrow_second_slot_candidates_by_phase({"RHCE": bg}, phased=True)
        dropped = {
            allele.genotype
            for allele in bg.filtered_out["narrow_second_slot_candidates_by_phase"]
        }
        self.assertEqual(dropped, {"RHCE*01.01", "RHCE*01.02.01"})

    def test_two_sides_means_the_subset_may_be_the_right_one(self):
        """Candidates describing different chromosomes must both survive."""
        bg = self._bg(
            self.CANDIDATES, self.POOL, {"a": "0|1", "b": "1|0", "c": "1/1"}
        )
        narrow_second_slot_candidates_by_phase({"RHCE": bg}, phased=True)
        self.assertEqual(len(bg.single_slot_genotypes), 3)

    def test_unphased_heterozygotes_leave_it_alone(self):
        """Ambiguity is the correct output when nothing locates the variants."""
        bg = self._bg(
            self.CANDIDATES, self.POOL, {"a": "0/1", "b": "0/1", "c": "1/1"}
        )
        narrow_second_slot_candidates_by_phase({"RHCE": bg}, phased=True)
        self.assertEqual(len(bg.single_slot_genotypes), 3)

    def test_does_nothing_without_the_phased_flag(self):
        bg = self._bg(self.CANDIDATES, self.POOL, self.ONE_SIDE)
        narrow_second_slot_candidates_by_phase({"RHCE": bg}, phased=False)
        self.assertEqual(len(bg.single_slot_genotypes), 3)

    def test_homozygous_variants_locate_nothing_and_are_not_consulted(self):
        """A homozygous variant is on both chromosomes, so its phase says nothing.

        Here the only heterozygous variants still share a side, so the narrowing must
        still happen even though the homozygous one reads '1/1'.
        """
        bg = self._bg(
            {"RHCE*01.01": {"c"}, "RHCE*01.20.04.02": {"a", "c"}},
            {"a": Zygosity.HET, "c": Zygosity.HOM},
            {"a": "0|1", "c": "1/1"},
        )
        narrow_second_slot_candidates_by_phase({"RHCE": bg}, phased=True)
        self.assertEqual(
            bg.single_slot_genotypes, [f"RHCE*01.20.04.02/{UNDETERMINED_SLOT}"]
        )

    def test_candidates_that_are_not_subsets_all_survive(self):
        """Nothing to choose between alleles neither of which contains the other."""
        bg = self._bg(
            {"RHCE*01.01": {"a"}, "RHCE*02": {"b"}},
            {"a": Zygosity.HET, "b": Zygosity.HET},
            {"a": "0|1", "b": "0|1"},
        )
        narrow_second_slot_candidates_by_phase({"RHCE": bg}, phased=True)
        self.assertEqual(len(bg.single_slot_genotypes), 2)

    def test_a_single_candidate_is_untouched(self):
        bg = self._bg({"RHCE*01": {"a"}}, {"a": Zygosity.HET}, {"a": "0|1"})
        narrow_second_slot_candidates_by_phase({"RHCE": bg}, phased=True)
        self.assertEqual(len(bg.single_slot_genotypes), 1)
