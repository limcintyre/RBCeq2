import unittest
from collections import defaultdict
from unittest.mock import MagicMock, patch

# --- Import the actual components from your project ---
from rbceq2.core_logic.alleles import Allele, BloodGroup, Pair
from rbceq2.core_logic.constants import AlleleState
from rbceq2.core_logic.utils import Zygosity

# Import the functions to be tested
from rbceq2.filters.phased import (
    _get_allele_phase_info,
    check_phase,
    filter_if_all_HET_vars_on_same_side_and_phased,
    filter_on_in_relationship_if_HET_vars_on_dif_side_and_phased,
    filter_pairs_by_phase,
    impossible_alleles_phased,
    iterate_over_list,
    remove_unphased,
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
        a1 = MockAllele(genotype="A1", defining_variants={"het1", "hom1"})
        a2 = MockAllele(genotype="A2", defining_variants={"het2"})
        pair_to_remove = Pair(a1, a2)
        self.mock_bg.alleles[AlleleState.NORMAL] = [pair_to_remove]
        self.mock_bg.variant_pool = {
            "het1": Zygosity.HET,
            "hom1": Zygosity.HOM,
            "het2": Zygosity.HET,
        }
        self.mock_bg.variant_pool_phase = {"het1": "0|1", "hom1": "1/1", "het2": "0|1"}
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
        a1_same_phase = MockAllele(genotype="A1", defining_variants={"het1"})
        a2_same_phase = MockAllele(genotype="A2", defining_variants={"het2"})
        a3_diff_phase = MockAllele(genotype="A3", defining_variants={"het3"})

        pair_to_remove = Pair(a1_same_phase, a2_same_phase)
        # Add a valid pair to ensure the "replace all with ref" logic is not triggered
        pair_to_keep = Pair(a1_same_phase, a3_diff_phase)

        self.mock_bg.alleles[AlleleState.NORMAL] = [pair_to_remove, pair_to_keep]
        self.mock_bg.variant_pool = {
            "het1": Zygosity.HET,
            "het2": Zygosity.HET,
            "het3": Zygosity.HET,
        }
        self.mock_bg.variant_pool_phase = {"het1": "0|1", "het2": "0|1", "het3": "1|0"}
        self.mock_bg.variant_pool_phase_set = {
            "het1": "setA",
            "het2": "setA",
            "het3": "setA",
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
        a1 = MockAllele(genotype="FUT2*01N.16", defining_variants={"v1"})
        a2 = MockAllele(genotype="FUT2*01N.02", defining_variants={"v2"})
        pair_to_remove = Pair(a1, a2)

        self.mock_bg.alleles[AlleleState.NORMAL] = [pair_to_remove]
        self.mock_bg.variant_pool = {"v1": Zygosity.HET, "v2": Zygosity.HET}
        self.mock_bg.variant_pool_phase = {"v1": "0|1", "v2": "0|1"}
        self.mock_bg.variant_pool_phase_set = {"v1": "setA", "v2": "setA"}

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


if __name__ == "__main__":
    unittest.main(argv=["first-arg-is-ignored"], exit=False)