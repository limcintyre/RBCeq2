import unittest
from collections import defaultdict
from unittest.mock import MagicMock, patch

# --- Import the actual components from your project ---
from rbceq2.core_logic.alleles import Allele, BloodGroup, Pair
from rbceq2.core_logic.constants import AlleleState
from rbceq2.core_logic.utils import Zygosity
# Import the functions to be tested
from rbceq2.filters.knops import (
    ensure_co_existing_HET_SNP_used,
    filter_co_existing_in_other_allele,
    filter_co_existing_pairs,
    filter_co_existing_subsets,
    filter_co_existing_with_normal,
    filter_coexisting_pairs_on_antithetical_zygosity,
    parse_bio_info2,
    remove_unphased_co,
)
# This is a dependency of the functions under test
from rbceq2.filters.shared_filter_functionality import flatten_alleles


# --- Helper Mock Class for Testing ---
class MockAllele(Allele):
    def __init__(self, **kwargs):
        defaults = {
            "phenotype": ".", "genotype_alt": ".", "phenotype_alt": ".",
            "defining_variants": frozenset(), "null": False, "weight_geno": 1000,
            "reference": False, "sub_type": "default",
        }
        defaults.update(kwargs)
        super().__init__(**defaults)

    def __eq__(self, other):
        return isinstance(other, Allele) and self.genotype == other.genotype

    def __hash__(self):
        return hash(self.genotype)

class MockPair:
    def __init__(self, allele1, allele2):
        self.alleles = (allele1, allele2)
        self.allele1 = allele1
        self.allele2 = allele2
    @property
    def comparable(self): return [frozenset(a.genotype.split("+")) for a in self.alleles]
    @property
    def contains_reference(self): return any(a.reference for a in self.alleles)
    @property
    def same_subtype(self): return self.allele1.sub_type == self.allele2.sub_type
    def __iter__(self): return iter(self.alleles)
    def __eq__(self, other): return isinstance(other, MockPair) and self.alleles == other.alleles
    def __hash__(self): return hash(self.alleles)
    def __repr__(self): return f"Pair(Genotype: {self.allele1.genotype}/{self.allele2.genotype})"


# --- Final, Consolidated, and Complete Test Class ---
class TestKnopsCoExistingFilters(unittest.TestCase):
    def setUp(self):
        """Set up a fresh mock BloodGroup object for each test."""
        self.mock_bg = MagicMock(spec=BloodGroup)
        self.mock_bg.alleles = {
            AlleleState.NORMAL: [], AlleleState.CO: [], AlleleState.RAW: [], AlleleState.FILT: []
        }
        self.mock_bg.filtered_out = defaultdict(list)
        self.mock_bg.variant_pool = {}
        self.mock_bg.variant_phase_set = {}
        self.mock_bg.type = "KN"

        def mock_remove(items_to_remove, filter_name, state=AlleleState.CO):
            self.mock_bg.filtered_out[filter_name].extend(items_to_remove)
            if self.mock_bg.alleles.get(state):
                self.mock_bg.alleles[state] = [item for item in self.mock_bg.alleles[state] if item not in items_to_remove]
        
        self.mock_bg.remove_pairs.side_effect = mock_remove
        self.mock_bg.remove_alleles.side_effect = mock_remove

    # --- Tests for filter_co_existing_with_normal ---
    def test_filter_with_normal_removes_pair_not_in_normal_list(self):
        pair_to_filter = MockPair(MockAllele(genotype="A"), MockAllele(genotype="B"))
        pair_to_keep = MockPair(MockAllele(genotype="C"), MockAllele(genotype="D"))
        self.mock_bg.alleles[AlleleState.NORMAL] = [pair_to_keep]
        self.mock_bg.alleles[AlleleState.CO] = [pair_to_filter, pair_to_keep]
        self.mock_bg.filtered_out[AlleleState.NORMAL] = [pair_to_filter]
        filter_co_existing_with_normal({1: self.mock_bg})
        self.mock_bg.remove_pairs.assert_called_once_with([pair_to_filter], "filter_co_existing_with_normal", AlleleState.CO)

    # --- Tests for parse_bio_info2 ---
    def test_parse_bio_info2(self):
        a1, a2, a3 = MockAllele(genotype="A+B"), MockAllele(genotype="C"), MockAllele(genotype="D")
        result = parse_bio_info2([MockPair(a1, a2), MockPair(a2, a3)])
        self.assertEqual(result, [[frozenset({"A", "B"}), frozenset({"C"})], [frozenset({"C"}), frozenset({"D"})]])

    # --- Tests for filter_co_existing_subsets ---
    def test_filter_subsets_removes_subset_pair(self):
        # CORRECTED TEST DATA
        aA, aB, aC = MockAllele(genotype="A"), MockAllele(genotype="B"), MockAllele(genotype="C")
        aAB = MockAllele(genotype="A+B")
        
        # This pair's flattened alleles are {'A', 'B'}
        pair_subset = MockPair(aA, aB)
        # This pair's flattened alleles are {'A', 'B', 'C'}, a proper superset
        pair_superset = MockPair(aAB, aC)
        
        self.mock_bg.alleles[AlleleState.CO] = [pair_subset, pair_superset]
        self.mock_bg.alleles[AlleleState.RAW] = [aA, aB, aC, aAB]
        self.mock_bg.variant_pool = {}
        
        filter_co_existing_subsets({1: self.mock_bg})
        
        self.mock_bg.remove_pairs.assert_called_once_with([pair_subset], "filter_co_existing_subsets", AlleleState.CO)

    def test_filter_subsets_does_not_remove_permutations(self):
        aA, aB = MockAllele(genotype="A"), MockAllele(genotype="B")
        pair1, pair2 = MockPair(aA, aB), MockPair(aB, aA)
        self.mock_bg.alleles[AlleleState.CO] = [pair1, pair2]
        self.mock_bg.alleles[AlleleState.RAW] = [aA, aB]
        filter_co_existing_subsets({1: self.mock_bg})
        self.mock_bg.remove_pairs.assert_not_called()

    # --- Tests for filter_coexisting_pairs_on_antithetical_zygosity ---
    def test_co_antithetical_removes_mismatched_subtypes(self):
        self.mock_bg.type = "SystemX"
        antitheticals = {"SystemX": ["S1", "S2"]}
        a_s1, a_s1_alt = MockAllele(genotype="A", sub_type="S1"), MockAllele(genotype="B", sub_type="S1")
        a_s2 = MockAllele(genotype="C", sub_type="S2")
        self.mock_bg.alleles[AlleleState.NORMAL] = [MockPair(a_s1, a_s2)]
        pair_to_remove = MockPair(a_s1, a_s1_alt)
        self.mock_bg.alleles[AlleleState.CO] = [pair_to_remove]
        filter_coexisting_pairs_on_antithetical_zygosity({1: self.mock_bg}, antitheticals)
        self.mock_bg.remove_pairs.assert_called_once_with([pair_to_remove], "filter_co_pairs_on_antithetical_zygosity", AlleleState.CO)

    # --- Tests for ensure_co_existing_HET_SNP_used ---
    @patch("rbceq2.filters.knops.check_var", return_value=2)
    def test_ensure_co_het_removes_pair_with_multiple_hits(self, mock_check_var):
        pair_ambiguous = MockPair(MockAllele(genotype="A"), MockAllele(genotype="B"))
        self.mock_bg.type = "KN"
        self.mock_bg.variant_pool = {"het_var": Zygosity.HET}
        self.mock_bg.alleles[AlleleState.CO] = [pair_ambiguous]
        ensure_co_existing_HET_SNP_used({1: self.mock_bg})
        mock_check_var.assert_called_with(self.mock_bg, pair_ambiguous, AlleleState.CO, "het_var")
        self.mock_bg.remove_pairs.assert_called_once_with([pair_ambiguous], "ensure_co_existing_HET_SNP_used", AlleleState.CO)

    # --- Tests for remove_unphased_co ---
    @patch("rbceq2.filters.knops.flatten_alleles")
    @patch("rbceq2.filters.knops.identify_unphased")
    def test_remove_unphased_co_removes_unphased_pair(self, mock_identify, mock_flatten):
        unphased_pair = MockPair(MockAllele(genotype="A"), MockAllele(genotype="B"))
        mock_identify.return_value = [unphased_pair.allele1]
        self.mock_bg.alleles[AlleleState.CO] = [unphased_pair]
        self.mock_bg.type = "KN"
        remove_unphased_co({1: self.mock_bg}, phased=True)
        self.mock_bg.remove_pairs.assert_called_once_with([unphased_pair], "remove_unphased_co", AlleleState.CO)

    # --- Tests for filter_co_existing_pairs ---
    def test_filter_co_pairs_removes_mushed_with_slash_in_subtype(self):
        mushed_pair = MockPair(MockAllele(genotype="A"), MockAllele(genotype="B+C", sub_type="KN*01/KN*02", genotype_alt="mushed"))
        self.mock_bg.alleles[AlleleState.CO] = [mushed_pair]
        filter_co_existing_pairs({1: self.mock_bg})
        self.mock_bg.remove_pairs.assert_called_once_with([mushed_pair], "filter_co_existing_pairs", AlleleState.CO)

    # --- Tests for filter_co_existing_in_other_allele ---
    def test_filter_co_in_other_removes_mushed_matching_real_allele_variants(self):
        mushed_allele = MockAllele(genotype="A+B", defining_variants={"v1", "v2"}, genotype_alt="mushed")
        real_allele = MockAllele(genotype="C_REAL", defining_variants={"v1", "v2"})
        pair_to_remove = MockPair(mushed_allele, MockAllele(genotype="D", defining_variants={"v3"}))
        self.mock_bg.alleles[AlleleState.CO] = [pair_to_remove]
        self.mock_bg.alleles[AlleleState.FILT] = [real_allele]
        filter_co_existing_in_other_allele({1: self.mock_bg})
        self.mock_bg.remove_pairs.assert_called_once_with([pair_to_remove], "filter_co_existing_in_other_allele", AlleleState.CO)


if __name__ == "__main__":
    unittest.main(argv=["first-arg-is-ignored"], exit=False)
    
# import unittest
# from collections import defaultdict
# from unittest.mock import MagicMock, patch

# # --- Import the actual components from your project ---
# from rbceq2.core_logic.alleles import Allele, BloodGroup
# from rbceq2.core_logic.constants import AlleleState
# from rbceq2.core_logic.utils import Zygosity

# # Import the functions to be tested
# from rbceq2.filters.knops import (
#     ensure_co_existing_HET_SNP_used,
#     filter_co_existing_in_other_allele,
#     filter_co_existing_pairs,
#     filter_co_existing_subsets,
#     filter_co_existing_with_normal,
#     filter_coexisting_pairs_on_antithetical_zygosity,
#     parse_bio_info2,
#     remove_unphased_co,
# )
# # This is a dependency of the functions under test


# # --- Helper Mock Class for Testing ---
# class MockAllele(Allele):
#     def __init__(self, **kwargs):
#         defaults = {
#             "phenotype": ".",
#             "genotype_alt": ".",
#             "phenotype_alt": ".",
#             "defining_variants": frozenset(),
#             "null": False,
#             "weight_geno": 1000,
#             "reference": False,
#             "sub_type": "default",
#         }
#         defaults.update(kwargs)
#         super().__init__(**defaults)

#     def __eq__(self, other):
#         return isinstance(other, Allele) and self.genotype == other.genotype

#     def __hash__(self):
#         return hash(self.genotype)


# class MockPair:
#     def __init__(self, allele1, allele2):
#         self.alleles = (allele1, allele2)
#         self.allele1 = allele1
#         self.allele2 = allele2

#     @property
#     def comparable(self):
#         return [frozenset(a.genotype.split("+")) for a in self.alleles]

#     @property
#     def contains_reference(self):
#         return any(a.reference for a in self.alleles)

#     @property
#     def same_subtype(self):
#         return self.allele1.sub_type == self.allele2.sub_type

#     def __iter__(self):
#         return iter(self.alleles)

#     def __eq__(self, other):
#         return isinstance(other, MockPair) and self.alleles == other.alleles

#     def __hash__(self):
#         return hash(self.alleles)

#     def __repr__(self):
#         return f"Pair(Genotype: {self.allele1.genotype}/{self.allele2.genotype})"


# # --- Final, Consolidated, and Complete Test Class ---
# class TestKnopsCoExistingFilters(unittest.TestCase):
#     def setUp(self):
#         """Set up a fresh mock BloodGroup object for each test."""
#         self.mock_bg = MagicMock(spec=BloodGroup)
#         self.mock_bg.alleles = {
#             AlleleState.NORMAL: [],
#             AlleleState.CO: [],
#             AlleleState.RAW: [],
#             AlleleState.FILT: [],
#         }
#         self.mock_bg.filtered_out = defaultdict(list)
#         self.mock_bg.variant_pool = {}
#         self.mock_bg.variant_phase_set = {}
#         self.mock_bg.type = "KN"

#         def mock_remove(items_to_remove, filter_name, state=AlleleState.CO):
#             self.mock_bg.filtered_out[filter_name].extend(items_to_remove)
#             if self.mock_bg.alleles.get(state):
#                 self.mock_bg.alleles[state] = [
#                     item
#                     for item in self.mock_bg.alleles[state]
#                     if item not in items_to_remove
#                 ]

#         self.mock_bg.remove_pairs.side_effect = mock_remove
#         self.mock_bg.remove_alleles.side_effect = mock_remove

#     # --- Tests for filter_co_existing_with_normal ---
#     def test_filter_with_normal_removes_pair_not_in_normal_list(self):
#         pair_to_filter = MockPair(MockAllele(genotype="A"), MockAllele(genotype="B"))
#         pair_to_keep = MockPair(MockAllele(genotype="C"), MockAllele(genotype="D"))
#         self.mock_bg.alleles[AlleleState.NORMAL] = [pair_to_keep]
#         self.mock_bg.alleles[AlleleState.CO] = [pair_to_filter, pair_to_keep]
#         self.mock_bg.filtered_out[AlleleState.NORMAL] = [pair_to_filter]
#         filter_co_existing_with_normal({1: self.mock_bg})
#         self.mock_bg.remove_pairs.assert_called_once_with(
#             [pair_to_filter], "filter_co_existing_with_normal", AlleleState.CO
#         )

#     def test_filter_with_normal_skips_mushed_alleles(self):
#         mushed_allele = MockAllele(genotype="A", genotype_alt="mushed")
#         pair_to_skip = MockPair(mushed_allele, MockAllele(genotype="B"))
#         self.mock_bg.alleles[AlleleState.NORMAL] = []
#         self.mock_bg.alleles[AlleleState.CO] = [pair_to_skip]
#         filter_co_existing_with_normal({1: self.mock_bg})
#         self.mock_bg.remove_pairs.assert_not_called()

#     # --- Tests for parse_bio_info2 ---
#     def test_parse_bio_info2(self):
#         a1, a2, a3 = (
#             MockAllele(genotype="A+B"),
#             MockAllele(genotype="C"),
#             MockAllele(genotype="D"),
#         )
#         result = parse_bio_info2([MockPair(a1, a2), MockPair(a2, a3)])
#         self.assertEqual(
#             result,
#             [
#                 [frozenset({"A", "B"}), frozenset({"C"})],
#                 [frozenset({"C"}), frozenset({"D"})],
#             ],
#         )

#     # --- Tests for filter_co_existing_subsets ---
#     def test_filter_subsets_removes_subset_pair(self):
#         aA, aB, aAB = (
#             MockAllele(genotype="A", defining_variants={"vA"}),
#             MockAllele(genotype="B", defining_variants={"vB"}),
#             MockAllele(genotype="A+B", defining_variants={"vA", "vB"}),
#         )
#         pair_subset, pair_superset = MockPair(aA, aB), MockPair(aAB, aB)
#         self.mock_bg.alleles[AlleleState.CO] = [pair_subset, pair_superset]
#         self.mock_bg.alleles[AlleleState.RAW] = [aA, aB, aAB]
#         self.mock_bg.variant_pool = {"vA": Zygosity.HOM, "vB": Zygosity.HOM}
#         filter_co_existing_subsets({1: self.mock_bg})
#         self.mock_bg.remove_pairs.assert_called_once_with(
#             [pair_subset], "filter_co_existing_subsets", AlleleState.CO
#         )

#     def test_filter_subsets_does_not_remove_permutations(self):
#         aA, aB = (
#             MockAllele(genotype="A", defining_variants={"vA"}),
#             MockAllele(genotype="B", defining_variants={"vB"}),
#         )
#         pair1, pair2 = MockPair(aA, aB), MockPair(aB, aA)
#         self.mock_bg.alleles[AlleleState.CO] = [pair1, pair2]
#         self.mock_bg.alleles[AlleleState.RAW] = [aA, aB]
#         self.mock_bg.variant_pool = {"vA": Zygosity.HOM, "vB": Zygosity.HOM}
#         filter_co_existing_subsets({1: self.mock_bg})
#         self.mock_bg.remove_pairs.assert_not_called()

#     # --- Tests for filter_coexisting_pairs_on_antithetical_zygosity ---
#     def test_co_antithetical_removes_mismatched_subtypes(self):
#         self.mock_bg.type = "SystemX"
#         antitheticals = {"SystemX": ["S1", "S2"]}
#         a_s1, a_s1_alt = (
#             MockAllele(genotype="A", sub_type="S1"),
#             MockAllele(genotype="B", sub_type="S1"),
#         )
#         a_s2 = MockAllele(genotype="C", sub_type="S2")
#         self.mock_bg.alleles[AlleleState.NORMAL] = [MockPair(a_s1, a_s2)]
#         pair_to_remove = MockPair(a_s1, a_s1_alt)
#         self.mock_bg.alleles[AlleleState.CO] = [pair_to_remove]
#         filter_coexisting_pairs_on_antithetical_zygosity(
#             {1: self.mock_bg}, antitheticals
#         )
#         self.mock_bg.remove_pairs.assert_called_once_with(
#             [pair_to_remove], "filter_co_pairs_on_antithetical_zygosity", AlleleState.CO
#         )

#     # --- Tests for ensure_co_existing_HET_SNP_used ---
#     @patch("rbceq2.filters.knops.check_var")
#     def test_ensure_co_het_removes_pair_with_multiple_hits(self, mock_check_var):
#         mock_check_var.return_value = 2
#         pair_ambiguous = MockPair(
#             MockAllele(genotype="A", defining_variants={"v1"}),
#             MockAllele(genotype="B", defining_variants={"v2"}),
#         )
#         self.mock_bg.type = "KN"
#         self.mock_bg.variant_pool = {"het_var": Zygosity.HET}
#         self.mock_bg.alleles[AlleleState.CO] = [pair_ambiguous]
#         ensure_co_existing_HET_SNP_used({1: self.mock_bg})
#         mock_check_var.assert_called_with(
#             self.mock_bg, pair_ambiguous, AlleleState.CO, "het_var"
#         )
#         self.mock_bg.remove_pairs.assert_called_once_with(
#             [pair_ambiguous], "ensure_co_existing_HET_SNP_used", AlleleState.CO
#         )

#     # --- Tests for remove_unphased_co ---
#     @patch("rbceq2.filters.knops.flatten_alleles")
#     @patch("rbceq2.filters.knops.identify_unphased")
#     def test_remove_unphased_co_removes_unphased_pair(
#         self, mock_identify, mock_flatten
#     ):
#         unphased_pair = MockPair(
#             MockAllele(genotype="A", defining_variants={"vA"}),
#             MockAllele(genotype="B", defining_variants={"vB"}),
#         )
#         mock_identify.return_value = [unphased_pair.allele1]
#         self.mock_bg.alleles[AlleleState.CO] = [unphased_pair]
#         self.mock_bg.type = "KN"
#         remove_unphased_co({1: self.mock_bg}, phased=True)
#         self.mock_bg.remove_pairs.assert_called_once_with(
#             [unphased_pair], "remove_unphased_co", AlleleState.CO
#         )

#     # --- Tests for filter_co_existing_pairs ---
#     def test_filter_co_pairs_removes_mushed_with_slash_in_subtype(self):
#         mushed_pair = MockPair(
#             MockAllele(genotype="A"),
#             MockAllele(genotype="B+C", sub_type="KN*01/KN*02", genotype_alt="mushed"),
#         )
#         self.mock_bg.alleles[AlleleState.CO] = [mushed_pair]
#         filter_co_existing_pairs({1: self.mock_bg})
#         self.mock_bg.remove_pairs.assert_called_once_with(
#             [mushed_pair], "filter_co_existing_pairs", AlleleState.CO
#         )

#     # --- Tests for filter_co_existing_in_other_allele ---
#     def test_filter_co_in_other_removes_mushed_matching_real_allele_variants(self):
#         mushed_allele = MockAllele(
#             genotype="A+B", defining_variants={"v1", "v2"}, genotype_alt="mushed"
#         )
#         real_allele = MockAllele(genotype="C_REAL", defining_variants={"v1", "v2"})
#         pair_to_remove = MockPair(
#             mushed_allele, MockAllele(genotype="D", defining_variants={"v3"})
#         )
#         self.mock_bg.alleles[AlleleState.CO] = [pair_to_remove]
#         self.mock_bg.alleles[AlleleState.FILT] = [real_allele]
#         filter_co_existing_in_other_allele({1: self.mock_bg})
#         self.mock_bg.remove_pairs.assert_called_once_with(
#             [pair_to_remove], "filter_co_existing_in_other_allele", AlleleState.CO
#         )


# if __name__ == "__main__":
#     unittest.main(argv=["first-arg-is-ignored"], exit=False)

#TODO check all this is covered b4 del
# import unittest
# from unittest.mock import MagicMock
# from rbceq2.core_logic.alleles import Allele, BloodGroup, Pair

# from rbceq2.core_logic.constants import AlleleState
# from rbceq2.core_logic.utils import Zygosity
# from rbceq2.filters.knops import (
#     ensure_co_existing_HET_SNP_used,
#     filter_co_existing_in_other_allele,
#     filter_co_existing_pairs,
#     filter_co_existing_subsets,
#     filter_co_existing_with_normal,
#     filter_coexisting_pairs_on_antithetical_zygosity,
#     parse_bio_info2,
#     remove_unphased_co,
# )

# from collections import defaultdict


# # --- Mock Data Classes (retained for test isolation and simplicity) ---
# # Using these mock classes allows us to test the filter logic without
# # depending on the internal implementation of the real Allele/Pair classes.
# class MockAllele:
#     def __init__(
#         self,
#         genotype,
#         defining_variants,
#         sub_type=".",
#         genotype_alt=".",
#         reference=False,
#     ):
#         self.genotype = genotype
#         self.defining_variants = frozenset(defining_variants)
#         self.sub_type = sub_type
#         self.genotype_alt = genotype_alt
#         self.reference = reference

#     def __repr__(self):
#         return f"Allele(genotype='{self.genotype}')"

#     def __eq__(self, other):
#         return isinstance(other, MockAllele) and self.genotype == other.genotype

#     def __hash__(self):
#         return hash(self.genotype)


# class MockPair:
#     def __init__(self, allele1, allele2):
#         self.alleles = (allele1, allele2)
#         self.allele1 = allele1
#         self.allele2 = allele2

#     @property
#     def comparable(self):
#         return [frozenset(a.genotype.split("+")) for a in self.alleles]

#     @property
#     def contains_reference(self):
#         return any(a.reference for a in self.alleles)

#     @property
#     def same_subtype(self):
#         return self.allele1.sub_type == self.allele2.sub_type

#     def __iter__(self):
#         return iter(self.alleles)

#     def __eq__(self, other):
#         return isinstance(other, MockPair) and self.alleles == other.alleles

#     def __hash__(self):
#         return hash(self.alleles)

#     def __repr__(self):
#         return f"Pair(Genotype: {self.allele1.genotype}/{self.allele2.genotype})"


# # --- Final, Corrected Test Class ---
# class TestCoExistingFilters(unittest.TestCase):
#     def setUp(self):
#         """Set up a fresh mock BloodGroup object for each test."""
#         self.mock_bg = MagicMock(spec=BloodGroup)
#         self.mock_bg.alleles = {
#             AlleleState.NORMAL: [],
#             AlleleState.CO: [],
#             AlleleState.RAW: [],
#             AlleleState.FILT: [],
#         }
#         self.mock_bg.filtered_out = defaultdict(list)
#         self.mock_bg.variant_pool = {}
#         self.mock_bg.variant_phase_set = {}
#         self.mock_bg.type = "KN"

#     def test_filter_co_existing_with_normal(self):
#         """Test filtering a CO pair if its equivalent NORMAL pair was already filtered."""
#         allele1 = MockAllele("KN*01.06", {"varA"})
#         allele2 = MockAllele("KN*01.07", {"varB"})
#         pair_to_filter = MockPair(allele1, allele2)

#         # Setup: The pair exists in both states but is marked as filtered in NORMAL
#         self.mock_bg.alleles[AlleleState.NORMAL] = []  # It was removed from here
#         self.mock_bg.alleles[AlleleState.CO] = [pair_to_filter]
#         self.mock_bg.filtered_out[AlleleState.NORMAL] = [pair_to_filter]

#         # Action: Call the decorated function with a dictionary
#         filter_co_existing_with_normal({1: self.mock_bg})

#         # Assert: remove_pairs should be called on the CO state
#         self.mock_bg.remove_pairs.assert_called_once_with(
#             [pair_to_filter], "filter_co_existing_with_normal", AlleleState.CO
#         )

#     def test_parse_bio_info2(self):
#         """Test parsing of pairs into a list of frozensets (not decorated)."""
#         allele1 = MockAllele("A+B", {})
#         allele2 = MockAllele("C", {})
#         allele3 = MockAllele("D", {})
#         pair1 = MockPair(allele1, allele2)
#         pair2 = MockPair(allele2, allele3)

#         # Action: Call directly as it's not decorated
#         result = parse_bio_info2([pair1, pair2])

#         expected = [
#             [frozenset({"A", "B"}), frozenset({"C"})],
#             [frozenset({"C"}), frozenset({"D"})],
#         ]
#         self.assertEqual(result, expected)

#     def test_filter_co_existing_subsets(self):
#         """Test filtering of co-existing pairs that are subsets of others."""
#         allele_A = MockAllele("A", {"varA"})
#         allele_B = MockAllele("B", {"varB"})
#         allele_AB = MockAllele("A+B", {"varA", "varB"})
#         pair_subset = MockPair(allele_A, allele_B)
#         pair_superset = MockPair(allele_AB, allele_B)

#         self.mock_bg.alleles[AlleleState.CO] = [pair_subset, pair_superset]
#         self.mock_bg.alleles[AlleleState.RAW] = [allele_A, allele_B, allele_AB]
#         self.mock_bg.variant_pool = {"varA": Zygosity.HOM, "varB": Zygosity.HOM}

#         filter_co_existing_subsets({1: self.mock_bg})

#         self.mock_bg.remove_pairs.assert_called_once_with(
#             [pair_subset], "filter_co_existing_subsets", AlleleState.CO
#         )

#     def test_filter_coexisting_pairs_on_antithetical_zygosity(self):
#         """Test filtering based on antithetical zygosity rules."""
#         antitheticals = {"KN": ["KN*01", "KN*02"]}
#         allele_01 = MockAllele("A", {}, sub_type="KN*01")
#         allele_02 = MockAllele("B", {}, sub_type="KN*02")
#         allele_01_mut = MockAllele("C", {}, sub_type="KN*01")

#         valid_pair = MockPair(allele_01, allele_02)
#         invalid_pair = MockPair(allele_01, allele_01_mut)

#         # Setup: NORMAL alleles must have >1 subtype to trigger the filter logic.
#         # CO alleles contain the pair that should be filtered.
#         self.mock_bg.alleles[AlleleState.NORMAL] = [valid_pair]
#         self.mock_bg.alleles[AlleleState.CO] = [valid_pair, invalid_pair]

#         filter_coexisting_pairs_on_antithetical_zygosity(
#             {1: self.mock_bg}, antitheticals
#         )

#         self.mock_bg.remove_pairs.assert_called_once_with(
#             [invalid_pair], "filter_co_pairs_on_antithetical_zygosity", AlleleState.CO
#         )

#     def test_ensure_co_existing_HET_SNP_used(self):
#         """Test removal of pairs where a HET SNP could ambiguously form alleles."""
#         allele_A = MockAllele("A", {"var1"})
#         allele_B = MockAllele("B", {"var2"})
#         pair_ambiguous = MockPair(allele_A, allele_B)

#         self.mock_bg.variant_pool = {"het_var": Zygosity.HET}
#         self.mock_bg.alleles[AlleleState.CO] = [pair_ambiguous]

#         # Patch check_var where it's defined and make it return > 1
#         with unittest.mock.patch(
#             "rbceq2.filters.shared_filter_functionality.check_var", return_value=2
#         ) as mock_check_var:
#             ensure_co_existing_HET_SNP_used({1: self.mock_bg})

#         mock_check_var.assert_called_with(
#             self.mock_bg, pair_ambiguous, AlleleState.CO, "het_var"
#         )
#         self.mock_bg.remove_pairs.assert_called_once_with(
#             [pair_ambiguous], "ensure_co_existing_HET_SNP_used", AlleleState.CO
#         )

#     def test_remove_unphased_co(self):
#         """Test removal of co-existing pairs that are not correctly phased."""
#         allele1 = MockAllele("A", {"varA"})
#         allele2 = MockAllele("B", {"varB"})
#         unphased_pair = MockPair(allele1, allele2)

#         self.mock_bg.alleles[AlleleState.CO] = [unphased_pair]
#         self.mock_bg.type = "KN"  # Required for the function to proceed

#         # Patch identify_unphased where it's defined.
#         with unittest.mock.patch(
#             "rbceq2.filters.shared_filter_functionality.identify_unphased",
#             return_value=[allele1],
#         ) as mock_identify:
#             remove_unphased_co({1: self.mock_bg}, phased=True)

#         # Assert
#         self.mock_bg.remove_pairs.assert_called_once_with(
#             [unphased_pair], "remove_unphased_co", AlleleState.CO
#         )

#     def test_filter_co_existing_pairs(self):
#         """Test removal of pairs with invalid 'mushed' alleles containing a '/' in subtype."""
#         allele_good = MockAllele("A", {"varA"})
#         allele_mushed = MockAllele(
#             "B+C", {"varB", "varC"}, sub_type="KN*01/KN*02", genotype_alt="mushed"
#         )
#         pair_to_remove = MockPair(allele_good, allele_mushed)

#         self.mock_bg.alleles[AlleleState.CO] = [pair_to_remove]

#         filter_co_existing_pairs({1: self.mock_bg})

#         self.mock_bg.remove_pairs.assert_called_once_with(
#             [pair_to_remove], "filter_co_existing_pairs", AlleleState.CO
#         )

#     def test_filter_co_existing_in_other_allele(self):
#         """Test removal of a 'mushed' pair that matches a different existing allele's variants."""
#         mushed_variants = {"varA", "varB"}
#         allele_mushed = MockAllele("A+B", mushed_variants, genotype_alt="mushed")
#         allele_real = MockAllele(
#             "C_real", mushed_variants
#         )  # Same variants, different genotype
#         allele_other = MockAllele("D", {"varD"})
#         pair_to_remove = MockPair(allele_mushed, allele_other)

#         self.mock_bg.alleles[AlleleState.CO] = [pair_to_remove]
#         self.mock_bg.alleles[AlleleState.FILT] = [allele_real]

#         filter_co_existing_in_other_allele({1: self.mock_bg})

#         self.mock_bg.remove_pairs.assert_called_once_with(
#             [pair_to_remove], "filter_co_existing_in_other_allele", AlleleState.CO
#         )


# class TestEnsureCoExistingHetSnpUsed(unittest.TestCase):
#     def test_no_co_alleles(self):
#         """Test when bg.alleles[AlleleState.CO] is None."""
#         bg = BloodGroup(
#             type="KN",
#             alleles={AlleleState.CO: None},
#             sample="sample1",
#             variant_pool={"variant1": Zygosity.HET},
#             filtered_out=defaultdict(list),
#         )
#         result_bg = list(ensure_co_existing_HET_SNP_used({1: bg}).values())[0]
#         self.assertIs(result_bg, bg)
#         self.assertIsNone(result_bg.alleles[AlleleState.CO])

#     def test_no_het_variants(self):
#         """Test when there are no heterozygous variants in variant_pool."""
#         allele1 = Allele(
#             genotype="Allele1",
#             phenotype="Phenotype1",
#             genotype_alt=".",
#             phenotype_alt=".",
#             defining_variants=frozenset({"var1"}),
#             null=False,
#             weight_geno=1,
#             reference=False,
#             sub_type="Subtype1",
#         )
#         allele2 = Allele(
#             genotype="Allele2",
#             phenotype="Phenotype2",
#             genotype_alt=".",
#             phenotype_alt=".",
#             defining_variants=frozenset({"var2"}),
#             null=False,
#             weight_geno=1,
#             reference=False,
#             sub_type="Subtype1",
#         )
#         pair = Pair(allele1=allele1, allele2=allele2)
#         bg = BloodGroup(
#             type="KN",
#             alleles={AlleleState.CO: [pair]},
#             sample="sample2",
#             variant_pool={"variant3": Zygosity.HOM},
#             filtered_out=defaultdict(list),
#         )
#         result_bg = list(ensure_co_existing_HET_SNP_used({1: bg}).values())[0]
#         self.assertIs(result_bg, bg)
#         self.assertIn(pair, result_bg.alleles[AlleleState.CO])
#         self.assertNotIn(pair, result_bg.filtered_out["ensure_HET_SNP_used_CO"])

#     def test_het_variants_no_matching_alleles(self):
#         """Test when there are heterozygous variants but no matching alleles are found."""
#         allele1 = Allele(
#             genotype="Allele1",
#             phenotype="Phenotype1",
#             genotype_alt=".",
#             phenotype_alt=".",
#             defining_variants=frozenset({"var1"}),
#             null=False,
#             weight_geno=1,
#             reference=False,
#             sub_type="Subtype1",
#         )
#         allele2 = Allele(
#             genotype="Allele2",
#             phenotype="Phenotype2",
#             genotype_alt=".",
#             phenotype_alt=".",
#             defining_variants=frozenset({"var2"}),
#             null=False,
#             weight_geno=1,
#             reference=False,
#             sub_type="Subtype1",
#         )
#         pair = Pair(allele1=allele1, allele2=allele2)
#         allele3 = Allele(
#             genotype="Allele3",
#             phenotype="Phenotype3",
#             genotype_alt=".",
#             phenotype_alt=".",
#             defining_variants=frozenset({"var3"}),
#             null=False,
#             weight_geno=1,
#             reference=False,
#             sub_type="Subtype1",
#         )
#         pair_other = Pair(allele1=allele3, allele2=allele3)
#         bg = BloodGroup(
#             type="KN",
#             alleles={AlleleState.CO: [pair, pair_other]},
#             sample="sample3",
#             variant_pool={"variant1": Zygosity.HET},
#             filtered_out=defaultdict(list),
#         )
#         result_bg = list(ensure_co_existing_HET_SNP_used({1: bg}).values())[0]
#         self.assertIs(result_bg, bg)
#         self.assertIn(pair, result_bg.alleles[AlleleState.CO])
#         self.assertIn(pair_other, result_bg.alleles[AlleleState.CO])
#         self.assertNotIn(pair, result_bg.filtered_out["ensure_HET_SNP_used_CO"])

#     def test_hits_not_greater_than_one(self):
#         """Test when hits are not greater than 1, pairs are not removed."""
#         allele1 = Allele(
#             genotype="Allele1",
#             phenotype="Phenotype1",
#             genotype_alt=".",
#             phenotype_alt=".",
#             defining_variants=frozenset({"var1"}),
#             null=False,
#             weight_geno=1,
#             reference=False,
#             sub_type="Subtype1",
#         )
#         allele2 = Allele(
#             genotype="Allele2",
#             phenotype="Phenotype2",
#             genotype_alt=".",
#             phenotype_alt=".",
#             defining_variants=frozenset({"var2"}),
#             null=False,
#             weight_geno=1,
#             reference=False,
#             sub_type="Subtype1",
#         )
#         pair = Pair(allele1=allele1, allele2=allele2)
#         allele3 = Allele(
#             genotype="Allele3",
#             phenotype="Phenotype3",
#             genotype_alt=".",
#             phenotype_alt=".",
#             defining_variants=frozenset({"var1", "variant1"}),
#             null=False,
#             weight_geno=1,
#             reference=False,
#             sub_type="Subtype1",
#         )
#         pair_other = Pair(allele1=allele3, allele2=allele3)
#         bg = BloodGroup(
#             type="KN",
#             alleles={AlleleState.CO: [pair, pair_other]},
#             sample="sample4",
#             variant_pool={"variant1": Zygosity.HET},
#             filtered_out=defaultdict(list),
#         )
#         result_bg = list(ensure_co_existing_HET_SNP_used({1: bg}).values())[0]
#         self.assertIs(result_bg, bg)
#         self.assertIn(pair, result_bg.alleles[AlleleState.CO])
#         self.assertNotIn(pair, result_bg.filtered_out["ensure_HET_SNP_used_CO"])

#     def test_hits_greater_than_one_pair_removed(self):
#         """Test when hits > 1, pairs are removed."""
#         allele1 = Allele(
#             genotype="Allele1",
#             phenotype="Phenotype1",
#             genotype_alt=".",
#             phenotype_alt=".",
#             defining_variants=frozenset({"var1"}),
#             null=False,
#             weight_geno=1,
#             reference=False,
#             sub_type="Subtype1",
#         )
#         allele2 = Allele(
#             genotype="Allele2",
#             phenotype="Phenotype2",
#             genotype_alt=".",
#             phenotype_alt=".",
#             defining_variants=frozenset({"var2"}),
#             null=False,
#             weight_geno=1,
#             reference=False,
#             sub_type="Subtype1",
#         )
#         pair = Pair(allele1=allele1, allele2=allele2)
#         allele3 = Allele(
#             genotype="Allele3",
#             phenotype="Phenotype3",
#             genotype_alt=".",
#             phenotype_alt=".",
#             defining_variants=frozenset({"var1", "variant1"}),
#             null=False,
#             weight_geno=1,
#             reference=False,
#             sub_type="Subtype1",
#         )
#         allele4 = Allele(
#             genotype="Allele4",
#             phenotype="Phenotype4",
#             genotype_alt=".",
#             phenotype_alt=".",
#             defining_variants=frozenset({"var2", "variant1"}),
#             null=False,
#             weight_geno=1,
#             reference=False,
#             sub_type="Subtype1",
#         )
#         pair_other1 = Pair(allele1=allele3, allele2=allele3)
#         pair_other2 = Pair(allele1=allele4, allele2=allele4)
#         bg = BloodGroup(
#             type="KN",
#             alleles={AlleleState.CO: [pair, pair_other1, pair_other2]},
#             sample="sample5",
#             variant_pool={"variant1": Zygosity.HET},
#             filtered_out=defaultdict(list),
#         )
#         result_bg = list(ensure_co_existing_HET_SNP_used({1: bg}).values())[0]
#         self.assertIs(result_bg, bg)
#         self.assertNotIn(pair, result_bg.alleles[AlleleState.CO])
#         self.assertIn(pair, result_bg.filtered_out["ensure_co_existing_HET_SNP_used"])
#         self.assertIn(pair_other1, result_bg.alleles[AlleleState.CO])
#         self.assertIn(pair_other2, result_bg.alleles[AlleleState.CO])


# class TestFilterCoExistingInOtherAllele(unittest.TestCase):
#     def test_pairs_removed(self):
#         allele1 = Allele(
#             genotype="FUT2*01.03.01",
#             phenotype="Type1",
#             genotype_alt="mushed",
#             phenotype_alt=".",
#             defining_variants=frozenset({"var1"}),
#             null=False,
#             weight_geno=1,
#             reference=False,
#             sub_type="subtype1",
#         )
#         allele2 = Allele(
#             genotype="FUT2*01.06",
#             phenotype="Type2",
#             genotype_alt=".",
#             phenotype_alt=".",
#             defining_variants=frozenset({"var1"}),
#             null=False,
#             weight_geno=2,
#             reference=False,
#             sub_type="subtype2",
#         )
#         pair = Pair(allele1=allele1, allele2=allele2)
#         raw_allele = Allele(
#             genotype="FUT2*01.03.03",
#             phenotype="Type3",
#             genotype_alt=".",
#             phenotype_alt=".",
#             defining_variants=frozenset({"var1"}),
#             null=False,
#             weight_geno=3,
#             reference=False,
#             sub_type="subtype3",
#         )
#         bg = BloodGroup(
#             type="KN",
#             alleles={AlleleState.CO: [pair], AlleleState.FILT: [raw_allele]},
#             sample="sample1",
#             variant_pool={},
#         )
#         bg_filtered = list(filter_co_existing_in_other_allele({1: bg}).values())[0]
#         self.assertFalse(bg_filtered.alleles[AlleleState.CO])

#     def test_no_removal_required(self):
#         allele1 = Allele(
#             genotype="FUT2*01.03.01",
#             phenotype="Type1",
#             genotype_alt=".",
#             phenotype_alt=".",
#             defining_variants=frozenset({"var2"}),
#             null=False,
#             weight_geno=1,
#             reference=False,
#             sub_type="subtype1",
#         )
#         allele2 = Allele(
#             genotype="FUT2*01.06",
#             phenotype="Type2",
#             genotype_alt=".",
#             phenotype_alt=".",
#             defining_variants=frozenset({"var3"}),
#             null=False,
#             weight_geno=2,
#             reference=False,
#             sub_type="subtype2",
#         )
#         pair = Pair(allele1=allele1, allele2=allele2)
#         raw_allele = Allele(
#             genotype="FUT2*01.03.03",
#             phenotype="Type3",
#             genotype_alt=".",
#             phenotype_alt=".",
#             defining_variants=frozenset({"var1"}),
#             null=False,
#             weight_geno=3,
#             reference=False,
#             sub_type="subtype3",
#         )
#         bg = BloodGroup(
#             type="KN",
#             alleles={AlleleState.CO: [pair], "raw": [raw_allele]},
#             sample="sample1",
#             variant_pool={},
#         )

#         bg_filtered = list(filter_co_existing_in_other_allele({1: bg}).values())[0]
#         self.assertTrue(bg_filtered.alleles[AlleleState.CO])

#     def test_empty_co_existing_list(self):
#         bg = BloodGroup(
#             type="KN",
#             alleles={AlleleState.CO: [], "raw": []},
#             sample="sample1",
#             variant_pool={},
#         )

#         bg_filtered = list(filter_co_existing_in_other_allele({1: bg}).values())[0]
#         self.assertEqual(bg_filtered.alleles[AlleleState.CO], [])


# class TestFilterCoExistingPairs(unittest.TestCase):
#     def test_normal_case(self):
#         allele1 = Allele(
#             genotype="A1",
#             phenotype="M",
#             genotype_alt="mushed",
#             phenotype_alt=".",
#             defining_variants=frozenset({"var1"}),
#             null=False,
#             weight_geno=1,
#             reference=False,
#             sub_type="A/B",
#         )
#         allele2 = Allele(
#             genotype="A2",
#             phenotype="N",
#             genotype_alt=".",
#             phenotype_alt=".",
#             defining_variants=frozenset({"var2"}),
#             null=False,
#             weight_geno=2,
#             reference=False,
#             sub_type="C/D",
#         )
#         pair = Pair(allele1=allele1, allele2=allele2)
#         bg = BloodGroup(
#             type="KN",
#             alleles={AlleleState.CO: [pair]},
#             sample="sample1",
#             variant_pool={},
#         )

#         bg_filtered = list(filter_co_existing_pairs({1: bg}).values())[0]
#         self.assertFalse(bg_filtered.alleles[AlleleState.CO])

#     def test_no_removal_required(self):
#         allele1 = Allele(
#             genotype="A1",
#             phenotype="M",
#             genotype_alt=".",
#             phenotype_alt=".",
#             defining_variants=frozenset({"var1"}),
#             null=False,
#             weight_geno=1,
#             reference=False,
#             sub_type="A/B",
#         )
#         allele2 = Allele(
#             genotype="A2",
#             phenotype="N",
#             genotype_alt=".",
#             phenotype_alt=".",
#             defining_variants=frozenset({"var2"}),
#             null=False,
#             weight_geno=2,
#             reference=False,
#             sub_type="C/D",
#         )
#         pair = Pair(allele1=allele1, allele2=allele2)
#         bg = BloodGroup(
#             type="KN",
#             alleles={AlleleState.CO: [pair]},
#             sample="sample1",
#             variant_pool={},
#         )

#         bg_filtered = list(filter_co_existing_pairs({1: bg}).values())[0]
#         self.assertTrue(bg_filtered.alleles[AlleleState.CO])

#     def test_empty_co_existing_list(self):
#         bg = BloodGroup(
#             type="KN",
#             alleles={AlleleState.CO: []},
#             sample="sample1",
#             variant_pool={},
#         )

#         bg_filtered = list(filter_co_existing_pairs({1: bg}).values())[0]
#         self.assertEqual(bg_filtered.alleles[AlleleState.CO], [])


# class TestFilterCoExistingSubsets(unittest.TestCase):
#     def setUp(self):
#         self.allele1 = Allele(
#             genotype="KN*01.06",
#             phenotype="Type2",
#             genotype_alt=".",
#             phenotype_alt=".",
#             defining_variants=frozenset(
#                 {"207782856_A_G, 207782916_A_T, 207782889_A_G, 207782931_A_G"}
#             ),
#             null=False,
#             weight_geno=1,
#             reference=False,
#             sub_type="subtype1",
#         )
#         self.allele2 = Allele(
#             genotype="KN*01.07",
#             phenotype="Type2",
#             genotype_alt=".",
#             phenotype_alt=".",
#             defining_variants=frozenset({"207782916_A_T, 207782889_A_G"}),
#             null=False,
#             weight_geno=1,
#             reference=False,
#             sub_type="subtype1",
#         )
#         self.allele3 = Allele(
#             genotype="KN*01.10",
#             phenotype="Type2",
#             genotype_alt=".",
#             phenotype_alt=".",
#             defining_variants=frozenset({" 207782916_A_T, 207782931_A_G"}),
#             null=False,
#             weight_geno=1,
#             reference=False,
#             sub_type="subtype1",
#         )
#         self.allele4 = Allele(
#             genotype="KN*01.07+KN*01.10",
#             phenotype="Type2",
#             genotype_alt="mushed",
#             phenotype_alt=".",
#             defining_variants=frozenset(
#                 {"207782916_A_T, 207782889_A_G, 207782931_A_G"}
#             ),
#             null=False,
#             weight_geno=1,
#             reference=False,
#             sub_type="subtype1",
#         )
#         self.pair1 = Pair(allele1=self.allele1, allele2=self.allele4)  # Not a subset
#         self.pair2 = Pair(allele1=self.allele4, allele2=self.allele4)  # Subset

#     def test_pairs_removed(self):
#         bg = BloodGroup(
#             type="KN",
#             alleles={AlleleState.CO: [self.pair1, self.pair2], AlleleState.RAW: []},
#             sample="sample230",
#             variant_pool={
#                 "1:207782856_A_G": "Heterozygous",
#                 "1:207782889_A_G": "Homozygous",
#                 "1:207782916_A_T": "Homozygous",
#                 "1:207782931_A_G": "Homozygous",
#             },
#             filtered_out=defaultdict(list),
#         )
#         self.assertTrue(self.pair2 in bg.alleles[AlleleState.CO])
#         filtered_bg = list(filter_co_existing_subsets({1: bg}).values())[0]

#         self.assertTrue(self.pair2 not in filtered_bg.alleles[AlleleState.CO])
#         self.assertTrue(
#             self.pair2 in filtered_bg.filtered_out["filter_co_existing_subsets"]
#         )

#     def test_no_pairs_removed(self):
#         bg = BloodGroup(
#             type="KN",
#             alleles={AlleleState.CO: [self.pair1, self.pair2], AlleleState.RAW: []},
#             sample="sample230",
#             variant_pool={},
#             filtered_out=defaultdict(list),
#         )
#         # Reverse the condition so no subsets exist
#         filtered_bg = list(filter_co_existing_subsets({1: bg}).values())[0]
#         self.assertTrue(self.pair1 in filtered_bg.alleles[AlleleState.CO])
#         self.assertTrue(
#             self.pair1 not in filtered_bg.filtered_out["filter_co_existing_subsets"]
#         )  ###TODO wtf - is 'in' not working for Pair objects?

#     def test_empty_co_existing_list(self):
#         bg = BloodGroup(
#             type="KN",
#             alleles={AlleleState.CO: []},
#             sample="sample230",
#             variant_pool={},
#             filtered_out=defaultdict(list),
#         )
#         filtered_bg = list(filter_co_existing_subsets({1: bg}).values())[0]
#         self.assertEqual(filtered_bg.alleles[AlleleState.CO], [])

#     def test_permutation_pairs(self):
#         # Create Allele objects
#         alleleA = Allele(
#             genotype="AlleleA",
#             phenotype="TypeA",
#             genotype_alt=".",
#             phenotype_alt=".",
#             defining_variants=frozenset({"varA"}),
#             null=False,
#             weight_geno=1,
#             reference=False,
#             sub_type="subtype1",
#         )
#         alleleB = Allele(
#             genotype="AlleleB",
#             phenotype="TypeB",
#             genotype_alt=".",
#             phenotype_alt=".",
#             defining_variants=frozenset({"varB"}),
#             null=False,
#             weight_geno=1,
#             reference=False,
#             sub_type="subtype1",
#         )

#         # Create pairs that are permutations of each other
#         pair1 = Pair(allele1=alleleA, allele2=alleleB)
#         pair2 = Pair(allele1=alleleB, allele2=alleleA)

#         # Prepare the BloodGroup object
#         bg = BloodGroup(
#             type="KN",
#             alleles={
#                 AlleleState.CO: [pair1, pair2],
#                 AlleleState.RAW: [alleleA, alleleB],
#             },
#             sample="sample231",
#             variant_pool={"varA": Zygosity.HOM, "varB": Zygosity.HOM},
#             filtered_out=defaultdict(list),
#         )

#         # Now run the function
#         filtered_bg = list(filter_co_existing_subsets({1: bg}).values())[0]

#         # Both pairs should remain because they are permutations
#         self.assertIn(pair1, filtered_bg.alleles[AlleleState.CO])
#         self.assertIn(pair2, filtered_bg.alleles[AlleleState.CO])

#     def test_alleles_with_one_het_variant_not_same_subtype(self):
#         # Create Allele objects
#         allele1 = Allele(
#             genotype="Allele1",
#             phenotype="Type1",
#             genotype_alt=".",
#             phenotype_alt=".",
#             defining_variants=frozenset({"var1"}),
#             null=False,
#             weight_geno=1,
#             reference=False,
#             sub_type="subtype1",
#         )
#         allele2 = Allele(
#             genotype="Allele2",
#             phenotype="Type2",
#             genotype_alt=".",
#             phenotype_alt=".",
#             defining_variants=frozenset({"var2", "var3"}),
#             null=False,
#             weight_geno=1,
#             reference=False,
#             sub_type="subtype2",  # Different subtype
#         )

#         pair = Pair(allele1=allele1, allele2=allele2)

#         # Prepare the BloodGroup object
#         bg = BloodGroup(
#             type="KN",
#             alleles={AlleleState.CO: [pair], AlleleState.RAW: [allele1, allele2]},
#             sample="sample232",
#             variant_pool={
#                 "var1": Zygosity.HET,
#                 "var2": Zygosity.HOM,
#                 "var3": Zygosity.HOM,
#             },
#             filtered_out=defaultdict(list),
#         )

#         # Run the function
#         filtered_bg = list(filter_co_existing_subsets({1: bg}).values())[0]

#         # The pair should not be removed due to the 'continue' in the branch
#         self.assertIn(pair, filtered_bg.alleles[AlleleState.CO])
#         self.assertNotIn(pair, filtered_bg.filtered_out["filter_co_existing_subsets"])

#     def test_contains_reference_pair_removed(self):
#         # Create Reference Allele
#         ref_allele = Allele(
#             genotype="RefAllele",
#             phenotype="RefType",
#             genotype_alt=".",
#             phenotype_alt=".",
#             defining_variants=frozenset({"ref_var"}),
#             null=False,
#             weight_geno=1,
#             reference=True,  # Mark as reference
#             sub_type="subtype1",
#         )

#         # Create Non-Reference Alleles
#         allele1 = Allele(
#             genotype="Allele1",
#             phenotype="Type1",
#             genotype_alt=".",
#             phenotype_alt=".",
#             defining_variants=frozenset({"var1"}),
#             null=False,
#             weight_geno=1,
#             reference=False,
#             sub_type="subtype1",
#         )
#         allele2 = Allele(
#             genotype="Allele2",
#             phenotype="Type2",
#             genotype_alt=".",
#             phenotype_alt=".",
#             defining_variants=frozenset({"var2"}),
#             null=False,
#             weight_geno=1,
#             reference=False,
#             sub_type="subtype1",
#         )

#         # Create Combined Allele
#         combined_allele = Allele(
#             genotype="RefAllele+Allele1",
#             phenotype="CombinedType",
#             genotype_alt=".",
#             phenotype_alt=".",
#             defining_variants=frozenset({"ref_var", "var1"}),
#             null=False,
#             weight_geno=1,
#             reference=False,
#             sub_type="subtype1",
#         )

#         # Create Pairs
#         pair_with_ref = Pair(allele1=ref_allele, allele2=allele1)
#         other_pair = Pair(allele1=combined_allele, allele2=allele2)

#         # Prepare the BloodGroup object
#         bg = BloodGroup(
#             type="KN",
#             alleles={
#                 AlleleState.CO: [pair_with_ref, other_pair],
#                 AlleleState.RAW: [ref_allele, allele1, allele2, combined_allele],
#             },
#             sample="sample233",
#             variant_pool={
#                 "ref_var": Zygosity.HOM,
#                 "var1": Zygosity.HOM,
#                 "var2": Zygosity.HOM,
#             },
#             filtered_out=defaultdict(list),
#         )

#         # Run the function
#         filtered_bg = list(filter_co_existing_subsets({1: bg}).values())[0]

#         # Assert that pair_with_ref is removed
#         self.assertNotIn(pair_with_ref, filtered_bg.alleles[AlleleState.CO])
#         self.assertIn(
#             pair_with_ref, filtered_bg.filtered_out["filter_co_existing_subsets"]
#         )

#         # Assert that other_pair remains
#         self.assertIn(other_pair, filtered_bg.alleles[AlleleState.CO])


# class TestFilterCoExistingWithNormal(unittest.TestCase):
#     def setUp(self):
#         self.allele1 = Allele(
#             genotype="KN1",
#             phenotype="Type1",
#             genotype_alt=".",
#             phenotype_alt=".",
#             defining_variants=frozenset({"var1"}),
#             null=False,
#             weight_geno=1,
#             reference=False,
#             sub_type="subtype1",
#         )

#         self.allele2 = Allele(
#             genotype="KN2",
#             phenotype="Type2",
#             genotype_alt=".",
#             phenotype_alt=".",
#             defining_variants=frozenset({"var2"}),
#             null=False,
#             weight_geno=2,
#             reference=False,
#             sub_type="subtype2",
#         )
#         # Normal pair
#         self.pair1 = Pair(allele1=self.allele1, allele2=self.allele2)
#         self.allele3 = Allele(
#             genotype="KN1+KN2",
#             phenotype="Type2",
#             genotype_alt="mushed",
#             phenotype_alt=".",
#             defining_variants=frozenset({"var1, var2"}),
#             null=False,
#             weight_geno=2,
#             reference=False,
#             sub_type="subtype2",
#         )
#         # Pair with coexisting alleles
#         self.pair2 = Pair(allele1=self.allele1, allele2=self.allele3)

#     def test_normal_case(self):
#         bg = BloodGroup(
#             type="KN",
#             alleles={AlleleState.NORMAL: [self.pair1], AlleleState.CO: [self.pair1]},
#             sample="sample230",
#             variant_pool={},
#             filtered_out=defaultdict(list),
#         )
#         self.assertTrue(self.pair1 in bg.alleles[AlleleState.NORMAL])
#         self.assertTrue(self.pair1 in bg.alleles[AlleleState.CO])
#         bg.remove_pairs([self.pair1], "normal_filter", AlleleState.NORMAL)
#         self.assertTrue(self.pair1 not in bg.alleles[AlleleState.NORMAL])
#         filtered_bg = list(filter_co_existing_with_normal({1: bg}).values())[0]
#         self.assertTrue(
#             self.pair1 in filtered_bg.filtered_out["filter_co_existing_with_normal"]
#         )
#         self.assertEqual(filtered_bg.alleles, {"co_existing": [], "pairs": []})

#     def test_no_removal_needed(self):
#         # Modify pair2 to meet exclusion criteria
#         bg = BloodGroup(
#             type="KN",
#             alleles={AlleleState.NORMAL: [self.pair1], AlleleState.CO: [self.pair2]},
#             sample="sample230",
#             variant_pool={},
#             filtered_out=defaultdict(list),
#         )
#         filtered_bg = list(filter_co_existing_with_normal({1: bg}).values())[0]
#         self.assertTrue(self.pair2 in filtered_bg.alleles[AlleleState.CO])
#         self.assertEqual(len(filtered_bg.alleles[AlleleState.CO]), 1)
#         self.assertEqual(len(filtered_bg.alleles[AlleleState.NORMAL]), 1)

#     def test_empty_co_existing_list(self):
#         bg = BloodGroup(
#             type="KN",
#             alleles={AlleleState.NORMAL: [self.pair1], AlleleState.CO: []},
#             sample="sample230",
#             variant_pool={},
#             filtered_out=defaultdict(list),
#         )
#         filtered_bg = list(filter_co_existing_with_normal({1: bg}).values())[0]
#         self.assertEqual(filtered_bg.alleles[AlleleState.CO], [])


# class TestFilterCoexistingPairsOnAntitheticalZygosity(unittest.TestCase):
#     def test_no_co_alleles(self):
#         """Test when bg.alleles[AlleleState.CO] is None."""
#         bg = BloodGroup(
#             type="KN",
#             alleles={AlleleState.CO: None},
#             sample="sample1",
#         )
#         antitheticals = {"ExampleType": ["antigen1", "antigen2"]}

#         # result_bg = filter_coexisting_pairs_on_antithetical_zygosity(bg, antitheticals)
#         result_bg = list(
#             filter_coexisting_pairs_on_antithetical_zygosity(
#                 {1: bg}, antitheticals
#             ).values()
#         )[0]

#         self.assertIs(result_bg, bg)
#         self.assertIsNone(result_bg.alleles[AlleleState.CO])

#     def test_bg_type_not_in_antitheticals(self):
#         """Test when bg.type is not in antitheticals."""
#         allele1 = Allele(
#             genotype="Allele1",
#             phenotype="Phenotype1",
#             genotype_alt=".",
#             phenotype_alt=".",
#             defining_variants=frozenset(),
#             null=False,
#             weight_geno=1,
#             reference=False,
#             sub_type="Subtype1",
#         )
#         allele2 = Allele(
#             genotype="Allele2",
#             phenotype="Phenotype2",
#             genotype_alt=".",
#             phenotype_alt=".",
#             defining_variants=frozenset(),
#             null=False,
#             weight_geno=1,
#             reference=False,
#             sub_type="Subtype2",
#         )
#         pair = Pair(allele1=allele1, allele2=allele2)

#         bg = BloodGroup(
#             type="OtherType",
#             alleles={AlleleState.CO: [pair], AlleleState.NORMAL: [allele1, allele2]},
#             sample="sample2",
#         )
#         antitheticals = {"ExampleType": ["antigen1", "antigen2"]}

#         # result_bg = filter_coexisting_pairs_on_antithetical_zygosity(bg, antitheticals)
#         result_bg = list(
#             filter_coexisting_pairs_on_antithetical_zygosity(
#                 {1: bg}, antitheticals
#             ).values()
#         )[0]
#         self.assertIs(result_bg, bg)
#         self.assertIn(pair, result_bg.alleles[AlleleState.CO])

#     def test_flattened_sub_types_length_one(self):
#         """Test when length of flattened_sub_types is 1."""
#         allele1 = Allele(
#             genotype="Allele1",
#             phenotype="Phenotype1",
#             genotype_alt=".",
#             phenotype_alt=".",
#             defining_variants=frozenset(),
#             null=False,
#             weight_geno=1,
#             reference=False,
#             sub_type="Subtype1",
#         )
#         allele2 = Allele(
#             genotype="Allele2",
#             phenotype="Phenotype2",
#             genotype_alt=".",
#             phenotype_alt=".",
#             defining_variants=frozenset(),
#             null=False,
#             weight_geno=1,
#             reference=False,
#             sub_type="Subtype1",
#         )
#         pair = Pair(allele1=allele1, allele2=allele2)

#         # For AlleleState.NORMAL, we use Pair objects
#         pair_normal = Pair(allele1=allele1, allele2=allele2)

#         bg = BloodGroup(
#             type="KN",
#             alleles={
#                 AlleleState.CO: [pair],
#                 AlleleState.NORMAL: [pair_normal],  # Now contains Pair objects
#             },
#             sample="sample3",
#         )
#         antitheticals = {"ExampleType": ["antigen1", "antigen2"]}

#         result_bg = list(
#             filter_coexisting_pairs_on_antithetical_zygosity(
#                 {1: bg}, antitheticals
#             ).values()
#         )[0]

#         self.assertIs(result_bg, bg)
#         self.assertIn(pair, result_bg.alleles[AlleleState.CO])

#     def test_pairs_with_flat_sub_types_equal(self):
#         """Test when flat_sub_types == flattened_sub_types."""
#         allele1 = Allele(
#             genotype="Allele1",
#             phenotype="Phenotype1",
#             genotype_alt=".",
#             phenotype_alt=".",
#             defining_variants=frozenset(),
#             null=False,
#             weight_geno=1,
#             reference=False,
#             sub_type="Subtype1",
#         )
#         allele2 = Allele(
#             genotype="Allele2",
#             phenotype="Phenotype2",
#             genotype_alt=".",
#             phenotype_alt=".",
#             defining_variants=frozenset(),
#             null=False,
#             weight_geno=1,
#             reference=False,
#             sub_type="Subtype2",
#         )
#         pair = Pair(allele1=allele1, allele2=allele2)

#         # For AlleleState.NORMAL, we use Pair objects
#         pair_normal = Pair(allele1=allele1, allele2=allele2)

#         bg = BloodGroup(
#             type="KN",
#             alleles={
#                 AlleleState.CO: [pair],
#                 AlleleState.NORMAL: [pair_normal],  # Now contains Pair objects
#             },
#             sample="sample4",
#         )
#         antitheticals = {"ExampleType": ["antigen1", "antigen2"]}

#         result_bg = list(
#             filter_coexisting_pairs_on_antithetical_zygosity(
#                 {1: bg}, antitheticals
#             ).values()
#         )[0]

#         self.assertIs(result_bg, bg)
#         self.assertIn(pair, result_bg.alleles[AlleleState.CO])

#     def test_pairs_with_flat_sub_types_not_equal(self):
#         """Test when flat_sub_types != flattened_sub_types."""
#         allele1 = Allele(
#             genotype="Allele1",
#             phenotype="Phenotype1",
#             genotype_alt=".",
#             phenotype_alt=".",
#             defining_variants=frozenset(),
#             null=False,
#             weight_geno=1,
#             reference=False,
#             sub_type="Subtype1",
#         )
#         allele2 = Allele(
#             genotype="Allele2",
#             phenotype="Phenotype2",
#             genotype_alt=".",
#             phenotype_alt=".",
#             defining_variants=frozenset(),
#             null=False,
#             weight_geno=1,
#             reference=False,
#             sub_type="Subtype1",
#         )
#         allele3 = Allele(
#             genotype="Allele3",
#             phenotype="Phenotype3",
#             genotype_alt=".",
#             phenotype_alt=".",
#             defining_variants=frozenset(),
#             null=False,
#             weight_geno=1,
#             reference=False,
#             sub_type="Subtype2",
#         )
#         # Create pairs for NORMAL state
#         pair_normal1 = Pair(allele1=allele1, allele2=allele2)
#         pair_normal2 = Pair(allele1=allele1, allele2=allele3)  # Different subtypes

#         pair_co = Pair(allele1=allele1, allele2=allele2)

#         bg = BloodGroup(
#             type="KN",
#             alleles={
#                 AlleleState.CO: [pair_co],
#                 AlleleState.NORMAL: [pair_normal1, pair_normal2],
#             },
#             sample="sample5",
#             filtered_out=defaultdict(list),
#         )
#         antitheticals = {"ExampleType": ["antigen1", "antigen2"]}

#         # Adjusted function call to match the decorator
#         result_bg = list(
#             filter_coexisting_pairs_on_antithetical_zygosity(
#                 {1: bg}, antitheticals
#             ).values()
#         )[0]

#         self.assertIs(result_bg, bg)
#         self.assertNotIn(pair_co, result_bg.alleles[AlleleState.CO])
#         self.assertIn(
#             pair_co, result_bg.filtered_out["filter_co_pairs_on_antithetical_zygosity"]
#         )


# class TestParseBioInfo2(unittest.TestCase):
#     def test_multiple_pairs(self):
#         allele1 = Allele(
#             genotype="A1+A2",
#             phenotype="M",
#             genotype_alt=".",
#             phenotype_alt=".",
#             defining_variants=frozenset({"var1"}),
#             null=False,
#             weight_geno=1,
#             reference=False,
#             sub_type="subtype1",
#         )
#         allele2 = Allele(
#             genotype="B1+B2",
#             phenotype="N",
#             genotype_alt=".",
#             phenotype_alt=".",
#             defining_variants=frozenset({"var2"}),
#             null=False,
#             weight_geno=2,
#             reference=False,
#             sub_type="subtype2",
#         )
#         pair1 = Pair(allele1=allele1, allele2=allele2)
#         pair2 = Pair(allele1=allele2, allele2=allele1)
#         result = parse_bio_info2([pair1, pair2])
#         expected = [
#             [frozenset({"A1", "A2"}), frozenset({"B1", "B2"})],
#             [frozenset({"B1", "B2"}), frozenset({"A1", "A2"})],
#         ]
#         self.assertEqual(result, expected)

#     def test_empty_list(self):
#         result = parse_bio_info2([])
#         expected = []
#         self.assertEqual(result, expected)


# if __name__ == "__main__":
#     unittest.main(argv=["first-arg-is-ignored"], exit=False)
