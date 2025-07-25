from __future__ import annotations

from rbceq2.core_logic.alleles import BloodGroup, Pair
from rbceq2.core_logic.constants import AlleleState
from rbceq2.core_logic.utils import Zygosity, apply_to_dict_values

from rbceq2.filters.shared_filter_functionality import (
    flatten_alleles,
    identify_unphased,
    check_var
)

from icecream import ic


@apply_to_dict_values
def filter_co_existing_with_normal(bg: BloodGroup) -> BloodGroup:
    """Filter co-existing allele pairs based on normal allele pairs.

    This function checks the co-existing allele pairs (AlleleState.CO) and removes
    those that are not present in the filtered normal allele pairs. Basically, if it
    was filtered out in any context it should be filtred out in all contexts.
    Allele pairs with any allele having 'mushed' as genotype_alt are ignored.

    Args:
        bg (BloodGroup): A BloodGroup object containing alleles, including normal and
            co-existing pairs, as well as filtered-out pairs.

    Returns:
        BloodGroup: The updated BloodGroup after filtering co-existing pairs.
    
    Example:

    KN*01.07/KN*01.07 was removed by ensure_HET_SNP_used, co that is carried over to the 
    co_existing set
    #Results:
    Genotypes count: 1
    Genotypes: KN*01.06/KN*01.07
    Phenotypes (numeric): 
    KN:1,-2,3,-4,5,6,7,8,-9,10,11,-12,13
    Phenotypes (alphanumeric):
    Kn(a+b-),McC(a+b+),Sl1-,Yk(a+),Vil+,Sl3+,KCAM-,KDAS+,DACY+,YCAD-,KNMB+

    #Data:
    Vars: 
    1:207782916_A_T : Homozygous
    1:207782769_ref : Homozygous
    1:207782856_A_G : Heterozygous
    1:207782931_A_G : Homozygous
    1:207782889_A_G : Homozygous

    Raw: 
    Allele 
    genotype: KN*01 
    defining_variants: 
            1:207782916_A_T
            1:207782769_ref 
    weight_geno: 1000 
    phenotype: KN:1,-2 or Kn(a+),Kn(b-) 
    reference: True 

    Allele 
    genotype: KN*01.06 
    defining_variants: 
            1:207782889_A_G
            1:207782856_A_G
            1:207782931_A_G 
    weight_geno: 1000 
    phenotype: KN:1,-2,-3,-4,6,7,-9,10 
    reference: False 

    Allele 
    genotype: KN*01.07 
    defining_variants: 
            1:207782889_A_G
            1:207782931_A_G 
    weight_geno: 1000 
    phenotype: KN:1,-2,-4,7,-9,10
    reference: False 

    Allele 
    genotype: KN*01.10 
    defining_variants: 
            1:207782931_A_G 
    weight_geno: 1000 
    phenotype: KN:1,-2,-9,10
    reference: False 

    #Filters applied:

    ensure_HET_SNP_used: 
    Pair(Genotype: KN*01.07/KN*01.07 
    Pair(Genotype: KN*01.07/KN*01.10

    filter_co_existing_with_normal: 
    Pair(Genotype: KN*01.07/KN*01.07
    """
    if bg.alleles[AlleleState.CO] is None:
        return bg
    to_remove = []
    filtered_pairs = [
        pair
        for pair in bg.alleles[AlleleState.NORMAL]
        if pair not in bg.filtered_out[AlleleState.NORMAL]
    ]
    for pair in bg.alleles[AlleleState.CO]:
        if any(a.genotype_alt == "mushed" for a in pair.alleles):
            continue
        if pair in filtered_pairs:
            continue
        to_remove.append(pair)

    if to_remove:
        bg.remove_pairs(to_remove, "filter_co_existing_with_normal", AlleleState.CO)

    return bg


def parse_bio_info2(pairs: list[Pair]) -> list[list[frozenset[str]]]:
    """Parse biological information string into a set of frozensets representing each
    side of '/'.

    Args:
        info (str): Biological information string.

    Returns:
        A set of frozensets, each representing unique substrings
    """

    return [pair.comparable for pair in pairs]

@apply_to_dict_values
def filter_co_existing_subsets(bg: BloodGroup) -> BloodGroup:
    """Filter co-existing allele pairs that are subsets of larger allele combinations.

    This filtering ensures that no combination smaller than either 2/2 or 1/3 (in this
    case) exists. For example, KN*01.-13+KN*01.06/KN*01.07 is a subset of
    KN*01.-13+KN*01.06+KN*01.12/KN*01.07 and therefore cannot exist.

    Args:
        bg (BloodGroup): A BloodGroup object containing allele pairs in various states.

    Returns:
        BloodGroup: The updated BloodGroup after filtering co-existing allele pairs.
        
    Example:
    ----------
    KN*01.06/KN*01.07+KN*01.10
    KN*01.07+KN*01.10/KN*01.07+KN*01.10 - not possible,
    KN*01.07+KN*01.10 + 207782856_A_G = KN*01.06

    KN*01.06 = 207782856_A_G, HET
              207782916_A_T,
              207782889_A_G,
              207782931_A_G
    KN*01.07 = 207782889_A_G,
              207782916_A_T
    KN*01.10 = 207782931_A_G,
              207782916_A_T

    variant_pool = {
        '1:207782856_A_G': 'Heterozygous',
        '1:207782889_A_G': 'Homozygous',
        '1:207782916_A_T': 'Homozygous',
        '1:207782931_A_G': 'Homozygous'
    }
    """

    def compare_to_all(pair_to_check: Pair, all_co_pairs: list[Pair]) -> bool:
        """
        Compare a pair against a list of other pairs to check if it is a proper subset.
        """
        # Flatten the allele genotypes from the pair we are checking into a single set.
        # e.g., Pair('A', 'B+C') -> {'A', 'B', 'C'}
        flat_alleles_to_check = {g for side in pair_to_check.comparable for g in side}

        for other_pair in all_co_pairs:
            if pair_to_check == other_pair:
                continue

            flat_other_alleles = {g for side in other_pair.comparable for g in side}

            # --- THE CORE FIX IS HERE ---
            # Use the proper subset operator '<' on the flattened sets of allele names.
            # This ensures the pair is strictly smaller than the one it's compared to,
            # which resolves the symmetrical removal bug.
            if flat_alleles_to_check < flat_other_alleles:
                # This secondary logic for HET variants is preserved from the original.
                different_alleles = [
                    allele
                    for allele in bg.alleles[AlleleState.RAW]
                    if allele.genotype in flat_other_alleles.difference(flat_alleles_to_check)
                ]
                alleles_with_1_HET_var = [
                    allele
                    for allele in different_alleles
                    if [
                        bg.variant_pool.get(variant) for variant in allele.defining_variants
                    ].count(Zygosity.HET)
                    == 1
                ]
                if alleles_with_1_HET_var and not pair_to_check.same_subtype:
                    continue
                
                # If the checks pass, we've found a superset, so this pair should be removed.
                return True
        
        return False

    if bg.alleles[AlleleState.CO] is None:
        return bg

    to_remove = []
    all_co_pairs = bg.alleles[AlleleState.CO]
    all_co_pairs_without_ref = [pair for pair in all_co_pairs if not pair.contains_reference]
    
    for pair in all_co_pairs:
        if pair.contains_reference:
            if compare_to_all(pair, all_co_pairs):
                to_remove.append(pair)
        elif compare_to_all(pair, all_co_pairs_without_ref):
            to_remove.append(pair)

    if to_remove:
        bg.remove_pairs(to_remove, "filter_co_existing_subsets", AlleleState.CO)

    return bg

#this was replaced by the above (AI generated) version after unit tests detected a logic flaw:
# You are absolutely correct. Thank you for providing the class definitions. This context is critical, and you've pinpointed the exact area of complexity that explains the behavior. My previous analysis was incomplete because it was missing this context.
# The issue is a fascinating and subtle conflict between two different ways of defining a "subset" relationship:
# The Scientific Definition (Allele.__contains__): The Allele class defines the in operator based on defining_variants. allele_A in allele_B is True only if allele_A's variants are a proper subset of allele_B's variants. This is the scientifically correct and intuitive way to determine if one allele is a simpler version of another.
# The String-Based Logic (in the filter): The filter_co_existing_subsets function, as written, completely ignores the Allele.__contains__ method. Instead, it uses its own internal helper function, is_subset, which operates on the Pair.comparable property. This property is derived from the allele's genotype string, splitting it by +.
# This is the root of the problem. The test fails because the function's internal logic for determining a "subset" is flawed and does not behave like a true mathematical subset operation.
# The Flaw in the Code's Logic
# Let's re-examine the test case with this new understanding:
# pair_subset: Contains alleles A and B. Its comparable property effectively represents the sets {'A'} and {'B'}.
# pair_superset: Contains alleles A+B and B. Its comparable property represents the sets {'A', 'B'} and {'B'}.
# When the function's internal is_subset logic runs:
# It compares pair_subset to pair_superset: It checks if the elements {'A'} and {'B'} from the subset exist within the flattened elements of the superset ({'A', 'B'}). They do. So, pair_subset is correctly marked for removal.
# It compares pair_superset to pair_subset: It checks if the elements {'A', 'B'} and {'B'} from the superset exist within the flattened elements of the subset ({'A', 'B'}). They do. The flawed logic concludes that the superset is a "subset" of the subset and incorrectly marks pair_superset for removal as well.
# The function should only remove a pair if it is a proper subset of another. The current implementation removes any pair that has overlapping allele components, leading to the incorrect symmetrical removal that your test caught.
# The Suggested Code Change
# The best way to fix the function is to make its internal logic behave like a proper subset check, respecting its design of working with the flattened genotype strings. We can achieve this by comparing the size of the flattened allele sets. A pair is only a subset if its total set of allele components is a proper subset of the other's.
# Here is the corrected version of the filter_co_existing_subsets function. This version fixes the flawed logic, and with this change, your original unit test will pass as expected.
# E2e tests seem fine with it but I'll keep this here for a few iterations just in case 

# @apply_to_dict_values  
# def filter_co_existing_subsets(bg: BloodGroup) -> BloodGroup:
#     """Filter co-existing allele pairs that are subsets of larger allele combinations.

#     This filtering ensures that no combination smaller than either 2/2 or 1/3 (in this
#     case) exists. For example, KN*01.-13+KN*01.06/KN*01.07 is a subset of
#     KN*01.-13+KN*01.06+KN*01.12/KN*01.07 and therefore cannot exist.

#     Args:
#         bg (BloodGroup): A BloodGroup object containing allele pairs in various states.

#     Returns:
#         BloodGroup: The updated BloodGroup after filtering co-existing allele pairs.
#      Example:
#     ----------
#     KN*01.06/KN*01.07+KN*01.10
#     KN*01.07+KN*01.10/KN*01.07+KN*01.10 - not possible,
#     KN*01.07+KN*01.10 + 207782856_A_G = KN*01.06

#     KN*01.06 = 207782856_A_G, HET
#               207782916_A_T,
#               207782889_A_G,
#               207782931_A_G
#     KN*01.07 = 207782889_A_G,
#               207782916_A_T
#     KN*01.10 = 207782931_A_G,
#               207782916_A_T

#     variant_pool = {
#         '1:207782856_A_G': 'Heterozygous',
#         '1:207782889_A_G': 'Homozygous',
#         '1:207782916_A_T': 'Homozygous',
#         '1:207782931_A_G': 'Homozygous'
#     }
#     """

#     def compare_to_all(pair_compare, comparables):
#         """Compare a pair against a list of pairs to check if it is a subset.

#         Args:
#             pair_compare (tuple): A pair representing the allele combination to compare.
#             comparables (list[tuple]): A list of allele pairs to compare against.

#         Returns:
#             bool: True if the pair is a subset of any other pair, False otherwise.
#         """

#         def number_of_alleles(pair_to_compare):
#             return sum(len(bit) for bit in pair_to_compare)

#         def is_subset(allele_to_compare):
#             return all(
#                 co_allele in other_pair[0] or co_allele in other_pair[1]
#                 for co_allele in allele_to_compare
#             )

#         for other_pair in comparables:
#             if pair_compare == other_pair:
#                 continue
#             no_alleles_in_mushed_pair = number_of_alleles(pair_compare)
#             no_alleles_in_other_mushed_pair = number_of_alleles(other_pair)
#             flat_alleles = flatten_alleles(pair_compare)
#             flat_other_alleles = flatten_alleles(other_pair)

#             if (
#                 no_alleles_in_mushed_pair == no_alleles_in_other_mushed_pair
#                 and flat_alleles == flat_other_alleles
#             ):
#                 continue
#                 # don't remove due to het permutations ie
#                 # KN*01.10', 'KN*01.-05+KN*01.07
#                 # KN*01.-05', 'KN*01.07+KN*01.10
#             if is_subset(pair_compare[0]) and is_subset(pair_compare[1]):
#                 different_alleles = [
#                     allele
#                     for allele in bg.alleles[AlleleState.RAW]
#                     if allele.genotype in flat_other_alleles.difference(flat_alleles)
#                 ]
#                 alleles_with_1_HET_var = [
#                     allele
#                     for allele in different_alleles
#                     if [
#                         bg.variant_pool[variant] for variant in allele.defining_variants
#                     ].count(Zygosity.HET)
#                     == 1
#                 ]
#                 if alleles_with_1_HET_var and not pair.same_subtype:
#                     continue
#                 return True
#         return False

#     if bg.alleles[AlleleState.CO] is None:
#         return bg
#     to_remove = []
#     all_comparable = parse_bio_info2(bg.alleles[AlleleState.CO])
#     all_comparable_without_ref = parse_bio_info2(
#         [pair for pair in bg.alleles[AlleleState.CO] if not pair.contains_reference]
#     )
#     for pair in bg.alleles[AlleleState.CO]:
#         if pair.contains_reference:
#             if compare_to_all(pair.comparable, all_comparable):
#                 to_remove.append(pair)
#         elif compare_to_all(pair.comparable, all_comparable_without_ref):
#             to_remove.append(pair)
#     if to_remove:
#         bg.remove_pairs(to_remove, "filter_co_existing_subsets", AlleleState.CO)

#     return bg


@apply_to_dict_values
def filter_coexisting_pairs_on_antithetical_zygosity(
    bg: BloodGroup, antitheticals: dict[str, list[str]]
) -> BloodGroup:
    """Process genetic data to identify alleles and genotypes.

    Args:
    ----------
    res (dict[str, list[Allele]]):
        A dictionary mapping genotypes to lists of Allele objects.
    variant_pool_numeric (dict[str, int]):
        A dictionary mapping variants to their counts.

    Returns:
    ----------
    list[str]:
        A list of called genotypes based on processed data.
    
    Example:
    
    Genotypes count: 4
    Genotypes: 
    KN*01.-05+KN*01.10/KN*02
    KN*01.-05/KN*02.10
    KN*01.10/KN*02  
    KN*01/KN*02.10
    # note, the last two options are the case where 1:207760773_C_T is with
    # 1:207782769_G_A, RBCeq2 considers this possible, albiet unlikely
    # use of phased data will resolve this

    Phenotypes (numeric): KN:1,2,3,4,5,-6,-7,8,9,10,11,-12,13
    Phenotypes (alphanumeric): 
    Kn(a+b+),McC(a+b-),Sl1+,Yk(a+),Vil-,Sl3+,KCAM+,KDAS+,DACY+,YCAD-,KNMB+

    #Data:
    Vars: 
    1:207782916_A_T : Homozygous
    1:207782769_ref : Heterozygous
    1:207760773_C_T : Heterozygous
    1:207782931_A_G : Heterozygous
    1:207782769_G_A : Heterozygous

    Raw: 
    Allele 
    genotype: KN*01 
    defining_variants: 
            1:207782769_ref
            1:207782916_A_T #HOM
    weight_geno: 1000 
    phenotype: KN:1,-2 or Kn(a+),Kn(b-) 
    reference: True 

    Allele 
    genotype: KN*01.-05 
    defining_variants: 
            1:207760773_C_T 
    weight_geno: 1000 
    phenotype: KN:1,-2,-5 or Kn(a+),Kn(b-),Yk(a-) 
    reference: False 

    Allele 
    genotype: KN*01.10 
    defining_variants: 
            1:207782931_A_G 
    weight_geno: 1000 
    phenotype: KN:1,-2,-9,10 or Kn(a+),Kn(b-),KCAM-,KDAS+ 
    reference: False 

    Allele 
    genotype: KN*02 
    defining_variants: 
            1:207782769_G_A 
    weight_geno: 1000 
    phenotype: KN:-1,2 or Kn(a-),Kn(b+) 
    reference: False 

    Allele 
    genotype: KN*02.10 
    defining_variants: 
            1:207782931_A_G
            1:207782769_G_A 
    weight_geno: 1000 
    phenotype: KN:-1,2,-9,10 or Kn(a-),Kn(b+),KCAM-,KDAS+ 
    reference: False 

    #Filters applied:

    filter_pairs_on_antithetical_zygosity: 
    Pair(Genotype: KN*01.-05/KN*01.10 
    Pair(Genotype: KN*01/KN*01.-05 
    Pair(Genotype: KN*01/KN*01.10

    filter_co_pairs_on_antithetical_zygosity: 
    Pair(Genotype: KN*01.-05/KN*01.10 
    Pair(Genotype: KN*01/KN*01.-05+KN*01.10 #the point of this filter
    Pair(Genotype: KN*01/KN*01.-05 
    Pair(Genotype: KN*01/KN*01.10

    """
    
    if bg.alleles[AlleleState.CO] is None:
        return bg
    if bg.type not in antitheticals:
        return bg
    
    to_remove = []
    flattened_sub_types = {
        allele.sub_type
        for allele in flatten_alleles(bg.alleles[AlleleState.NORMAL])
    }
    if len(flattened_sub_types) > 1:
        for pair in bg.alleles[AlleleState.CO]:
            flat_sub_types = {allele.sub_type for allele in pair}
            if flat_sub_types == flattened_sub_types:
                continue
            else:
                to_remove.append(pair)
    if to_remove:
        bg.remove_pairs(
        to_remove, "filter_co_pairs_on_antithetical_zygosity", AlleleState.CO
    )

    return bg



@apply_to_dict_values
def ensure_co_existing_HET_SNP_used(bg: BloodGroup) -> BloodGroup:
    """
    Ensures that heterozygous variants are utilized in allele pairs if they can form
    existing alleles.

    This function iterates over heterozygous variants in the variant pool of a
    BloodGroup object.
    For each heterozygous variant, it checks whether adding this variant to the defining
    variants of alleles in each pair results in alleles that already exist outside the
    pair. If multiple such matches are found, the pair is considered invalid and is
    removed from the allele pairs.

    Args:
        bg (BloodGroup): The BloodGroup object containing allele pairs and variant pool
        information.

    Returns:
        BloodGroup: The updated BloodGroup object with inconsistent allele pairs
        removed.

    Example:
        Suppose you have a BloodGroup object with allele pairs and a variant pool that
        includes heterozygous variants. This function will remove pairs where the
        heterozygous variants are not properly used in the allele definitions, but
        could form existing alleles when added.

    Genotypes count: 5
    Genotypes: 
    KN*01.-05+KN*01.10/KN*02
    KN*01.-05/KN*02##### not possible, 207782931_A_G must be on one side
    KN*01.-05/KN*02.10
    KN*01.10/KN*02
    KN*01/KN*02.10
    Phenotypes (numeric): KN:1,2,3,4,5,-6,-7,8,9,-10,11,-12,13
    KN:1,2,3,4,5,-6,-7,8,9,10,11,-12,13
    Phenotypes (alphanumeric): Kn(a+b+),McC(a+b-),Sl1+,Yk(a+),Vil-,Sl3+,KCAM+,KDAS+,DACY+,YCAD-,KNMB+
    Kn(a+b+),McC(a+b-),Sl1+,Yk(a+),Vil-,Sl3+,KCAM+,KDAS-,DACY+,YCAD-,KNMB+

    #Data:
    Vars: 
    1:207782916_A_T : Homozygous
    1:207782769_ref : Heterozygous
    1:207760773_C_T : Heterozygous
    1:207782931_A_G : Heterozygous
    1:207782769_G_A : Heterozygous
  
    Raw: 
    Allele 
    genotype: KN*01 
    defining_variants: 
            1:207782916_A_T
            1:207782769_ref 
    weight_geno: 1000 
    phenotype: KN:1,-2 or Kn(a+),Kn(b-) 
    reference: True 

    Allele 
    genotype: KN*01.-05 
    defining_variants: 
            1:207760773_C_T 
    weight_geno: 1000 
    phenotype: KN:1,-2,-5 or Kn(a+),Kn(b-),Yk(a-) 
    reference: False 

    Allele 
    genotype: KN*01.10 
    defining_variants: 
            1:207782931_A_G 
    weight_geno: 1000 
    phenotype: KN:1,-2,-9,10 or Kn(a+),Kn(b-),KCAM-,KDAS+ 
    reference: False 

    Allele 
    genotype: KN*02 
    defining_variants: 
            1:207782769_G_A 
    weight_geno: 1000 
    phenotype: KN:-1,2 or Kn(a-),Kn(b+) 
    reference: False 

    Allele 
    genotype: KN*02.10 
    defining_variants: 
            1:207782769_G_A
            1:207782931_A_G 
    weight_geno: 1000 
    phenotype: KN:-1,2,-9,10 or Kn(a-),Kn(b+),KCAM-,KDAS+ 
    reference: False 


    """
    if bg.type != "KN":
        return bg
    to_remove = []
    for variant, zygo in bg.variant_pool.items():
        if zygo == Zygosity.HET:
            for pair in bg.alleles[AlleleState.CO]:
                hits = check_var(bg, pair, AlleleState.CO, variant)
                if hits > 1:
                    to_remove.append(pair)

    if to_remove:
        bg.remove_pairs(to_remove, "ensure_co_existing_HET_SNP_used", AlleleState.CO)

    return bg


@apply_to_dict_values
def remove_unphased_co(bg: BloodGroup, phased: bool) -> BloodGroup:
    """Remove unphased alleles from the BloodGroup's bg.alleles[AlleleState.CO] 
    state if phased flag is set.

    This function iterates through the alleles in the FILT state and checks their
    phasing. Alleles with more than two distinct phases trigger a warning. If an
    allele has exactly two phases and no placeholder ('.') is present, it is marked
    for removal. Alleles with a single phase remain, as they are assumed to align
    with the reference.

    Args:
        bg (BloodGroup): A BloodGroup object containing allele states and phasing
            information.
        phased (bool): A flag indicating whether phasing should be enforced.

    Returns:
        BloodGroup: The updated BloodGroup with improperly phased alleles removed.
    
    Example:
            
    #Results:
    Genotypes count: 1
    Genotypes: KN*01.-05/KN*01.06
    Phenotypes (numeric): KN:1,-2,3,4,5,6,7,8,9,10,11,-12,13
    Phenotypes (alphanumeric): 
    Kn(a+b-),McC(a+b+),Sl1+,Yk(a+),Vil+,Sl3+,KCAM+,KDAS+,DACY+,YCAD-,KNMB+

    #Data:
    Vars: 
    1:207609424_ref : Homozygous
    1:207609571_A_T : Homozygous
    1:207587428_C_T : Heterozygous
    1:207609544_A_G : Heterozygous
    1:207609586_A_G : Heterozygous
    1:207609511_A_G : Heterozygous
    Vars_phase: 
    1:207609424_ref : 1/1
    1:207609571_A_T : 1/1
    1:207587428_C_T : 1|0
    1:207609544_A_G : 0|1
    1:207609586_A_G : 0|1
    1:207609511_A_G : 0|1
    Vars_phase_set: 
    1:207609424_ref : .
    1:207609571_A_T : .
    1:207587428_C_T : 202971590
    1:207609544_A_G : 202971590
    1:207609586_A_G : 202971590
    1:207609511_A_G : 202971590

    Raw: 
    Allele 
    genotype: KN*01 
    defining_variants: 
            1:207609571_A_T
            1:207609424_ref 
    weight_geno: 1000 
    phenotype: KN:1,-2 or Kn(a+),Kn(b-) 
    reference: True 

    Allele 
    genotype: KN*01.-05 
    defining_variants: 
            1:207587428_C_T 1|0
    weight_geno: 1000 
    phenotype: KN:1,-2,-5 or Kn(a+),Kn(b-),Yk(a-) 
    reference: False 

    Allele 
    genotype: KN*01.06 
    defining_variants: 
            1:207609586_A_G
            1:207609511_A_G
            1:207609544_A_G 
    weight_geno: 1000 
    phenotype: KN:1,-2,-3,-4,6,7,-9,10 or 
    Kn(a+),Kn(b-),Sl1-,McC(a-),McC(b+),Vil+,KCAM-,KDAS+ 
    reference: False 

    Allele 
    genotype: KN*01.07 
    defining_variants: 
            1:207609586_A_G
            1:207609544_A_G 
    weight_geno: 1000 
    phenotype: KN:1,-2,-4,7,-9,10 or
    Kn(a+),Kn(b-),Sl1-,Vil+,KCAM-,KDAS+ 
    reference: False 

    Allele 
    genotype: KN*01.10 
    defining_variants: 
            1:207609586_A_G 0|1
    weight_geno: 1000 
    phenotype: KN:1,-2,-9,10 or 
    Kn(a+),Kn(b-),KCAM-,KDAS+ 
    reference: False 

    #Filters applied:

    2025-07-24 11:42:24.409 | DEBUG    | 
    remove_unphased_co: 
    Pair(Genotype: KN*01/KN*01.-05+KN*01.06
    Pair(Genotype: KN*01/KN*01.-05+KN*01.07
    Pair(Genotype: KN*01/KN*01.-05+KN*01.10

    ie, KN*01.-05+KN*01.10 not possible:
        1:207609586_A_G = 0|1
        1:207587428_C_T 1|0
    """

    if not phased or bg.type != "KN":
        return bg
        
    co_alleles = bg.alleles[AlleleState.CO]
    if co_alleles is None:
        return bg
    to_remove = []
    alleles_to_remove = identify_unphased(bg, flatten_alleles(co_alleles))
    for pair in bg.alleles[AlleleState.CO]:
        if pair.allele1 in alleles_to_remove or pair.allele2 in alleles_to_remove:
            to_remove.append(pair)

    if to_remove:
        bg.remove_pairs(to_remove, "remove_unphased_co", AlleleState.CO)
    return bg


# doesn't come up as it was written to handle & in A4GALT
# keeping due to chance that future additions to KN will
# require same logic
@apply_to_dict_values
def filter_co_existing_pairs(bg: BloodGroup) -> BloodGroup:
    """Example:
    [Allele(genotype='A4GALT*01.02',
           genotype_alt='.',
           defining_variants=frozenset({'22:43089849_T_C'}),
           weight_geno=1000,
           reference=False,
           sub_type='A4GALT*01',
    Allele(genotype='A4GALT*02',
           genotype_alt='.',
           defining_variants=frozenset({'22:43113793_C_A'}),
           weight_geno=1000,
           reference=False,
           sub_type='A4GALT*02',
    Allele(genotype='A4GALT*02.02',
           genotype_alt='.',
           defining_variants=frozenset({'22:43089849_T_C',
                                        '22:43113793_C_A'}),
           weight_geno=1000,
           reference=False,
           sub_type='A4GALT*02',.

    variant_pool={'22:43089849_T_C': 'Heterozygous',
                   '22:43113793_C_A': 'Heterozygous'},

    filtered_out=defaultdict(<class 'list'>,  #modifying SNP has to be on one side
    {'filter_pairs_by_context': [[Allele(genotype='A4GALT*02',
                                         defining_variants=frozenset({'22:43113793_C_A'}),
                                  Allele(genotype='A4GALT*01',
                                         defining_variants=frozenset({'22:43113793_ref'}),

     'filter_pairs_on_antithetical_zygosity': #22:43113793 is HET so can't
     have A4GALT*01* on both sides
                                [[Allele(genotype='A4GALT*01.02',
                                        defining_variants=frozenset({'22:43089849_T_C'}),
                                Allele(genotype='A4GALT*01',
                                        defining_variants=frozenset({'22:43113793_ref'}),

    sample: 008Kenya A4GALT
    A4GALT*01.02/A4GALT*02
    A4GALT*01/A4GALT*02.02
     me:
     A4GALT*01.02/A4GALT*02
     A4GALT*01/A4GALT*02.02
     A4GALT*01/A4GALT*01.02&A4GALT*02 not possible:
        can't have A4GALT*01 and A4GALT*02 on same side,
        can't have ref/ref and
        have to have only A4GALT*01* on one side and only A4GALT*02* on the other side,
         when anithetical SNP HET.
        DIFFERENT logic to both existing filters
    """  # noqa: D401, D205
    if bg.alleles[AlleleState.CO] is not None:
        to_remove = []
        for pair in bg.alleles[AlleleState.CO]:
            for allele in pair:
                if allele.genotype_alt == "mushed" and "/" in allele.sub_type:
                    to_remove.append(pair)
                    break

        bg.remove_pairs(to_remove, "filter_co_existing_pairs", AlleleState.CO)

    return bg


# doesn't come up as its about & and FUT2*01.06
# keeping incase this logic is needed for KN one day
@apply_to_dict_values
def filter_co_existing_in_other_allele(bg: BloodGroup) -> BloodGroup:
    """Example:
    ----------
    1:

    226
    FUT2*01.03.01/FUT2*01.06
    FUT2*01/FUT2*01.03.03
    FUT2*01/FUT2*01.03.01&FUT2*01.06 - not possible
    FUT2*01.03.01 and FUT2*01.06 have 1 defining variant each,
    together they define FUT2*01.03.03
    """  # noqa: D401, D415, D400, D205
    if bg.alleles[AlleleState.CO] is not None:
        to_remove = []
        for pair in bg.alleles[AlleleState.CO]:
            for allele in pair:
                if allele.genotype_alt == "mushed":
                    for raw_allele in bg.alleles[AlleleState.FILT]:
                        if (
                            raw_allele.defining_variants == allele.defining_variants
                            and raw_allele.genotype != allele.genotype
                        ):
                            to_remove.append(pair)
        bg.remove_pairs(to_remove, "filter_co_existing_in_other_allele", AlleleState.CO)

    return bg
