from __future__ import annotations

from functools import partial
from collections.abc import Callable
from rbceq2.core_logic.alleles import BloodGroup, Pair
from rbceq2.core_logic.constants import (
    ABO_DELG_VARIANTS,
    SYNTHESISED_HOM_REF_GT,
    UNDETERMINED_SLOT,
    AlleleState,
)
from rbceq2.core_logic.utils import (
    Zygosity,
    apply_to_dict_values,
)
from rbceq2.filters.shared_filter_functionality import (
    flatten_alleles,
    all_hom,
    identify_unphased,
    proceed,
)
from rbceq2.core_logic.alleles import Allele
from icecream import ic

@apply_to_dict_values
def remove_unphased(bg: BloodGroup, phased: bool) -> BloodGroup:
    """Remove unphased alleles from the BloodGroup's FILT state if phased flag is set.

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

    #Example:
    Genotypes count: 1
    Genotypes: JK*01W.06/JK*02
    Phenotypes (numeric): JK:1w,2
    Phenotypes (alphanumeric): Jk(a+wb+)

    #Data:
    Vars:
    18:45739554_ref : Heterozygous
    18:45730450_G_A : Heterozygous
    18:45736573_A_G : Homozygous
    18:45739554_G_A : Heterozygous
    Vars_phase:
    18:45739554_ref : 1|0
    18:45730450_G_A : 1|0
    18:45736573_A_G : 1/1
    18:45739554_G_A : 0|1
    Vars_phase_set:
    18:45739554_ref : 20911244
    18:45730450_G_A : 20911244
    18:45736573_A_G : .
    18:45739554_G_A : 20911244


    #Filters applied:
    remove_unphased:
     Vars_phase:
    18:45739554_ref : 1|0
    18:45730450_G_A : 1|0
    18:45736573_A_G : 1/1
    18:45739554_G_A : 0|1
    Allele
    genotype: JK*02W.03
    defining_variants:
            18:45736573_A_G 1/1
            18:45730450_G_A 1|0
            18:45739554_G_A 0|1
    weight_geno: 1000
    phenotype: JK:-1,2w or Jk(a-),Jk(b+w)
    reference: False

    Allele
    genotype: JK*02W.04
    defining_variants:
            18:45730450_G_A 1|0
            18:45739554_G_A 0|1
    weight_geno: 1000
    phenotype: JK:-1,2w or Jk(a-),Jk(b+w)
    reference: False
    """

    if not phased:
        return bg

    to_remove = identify_unphased(bg, bg.alleles[AlleleState.FILT])
    if to_remove:
        bg.remove_alleles(to_remove, "remove_unphased", AlleleState.FILT)
    return bg


def _get_allele_phase_info(allele, phase_dict):
    """ """
    return [phase_dict[variant] for variant in allele.defining_variants]


@apply_to_dict_values
def filter_if_all_HET_vars_on_same_side_and_phased(
    bg: BloodGroup, phased: bool
) -> BloodGroup:
    """
    All HET vars on same side so can only have ref/HOMs on other side

    Example
    Sample: GM18501 BG Name: GYPB

    #Results:
        Genotypes count: 3
        Genotypes:
        GYPB*03/GYPB*03.06
        GYPB*03/GYPB*06.02
        GYPB*03.06/GYPB*06.02 # not possible
        Phenotypes (numeric): MNS:3,-4 | MNS:3,-4,6
        Phenotypes (alphanumeric): S+,s- | S+,s-,He+

    #Data:
    Vars:
    4:143999443_G_A : Homozygous
    4:143997537_C_T : Heterozygous
    4:144001261_T_C : Heterozygous
    4:144001254_T_A : Heterozygous
    4:144001262_A_C : Heterozygous
    4:144001249_C_A : Heterozygous
    4:144001250_T_C : Heterozygous
    Vars_phase:
    4:143999443_G_A : 1/1
    4:143997537_C_T : 0|1
    4:144001261_T_C : 0|1
    4:144001254_T_A : 0|1
    4:144001262_A_C : 0|1
    4:144001249_C_A : 0|1
    4:144001250_T_C : 0|1
    Vars_phase_set:
    4:143999443_G_A : .
    4:143997537_C_T : 129362934
    4:144001261_T_C : 129362934
    4:144001254_T_A : 129362934
    4:144001262_A_C : 129362934
    4:144001249_C_A : 129362934
    4:144001250_T_C : 129362934
    Raw:
    Allele
    genotype: GYPB*03
    defining_variants:
            4:143999443_G_A #hom
    weight_geno: 1000
    phenotype: MNS:3,-4 or S+,s-
    reference: False

    Allele
    genotype: GYPB*03.06
    defining_variants:
            4:143997537_C_T #this is on same side as
            4:143999443_G_A #hom
    weight_geno: 1000
    phenotype: MNS:3,-4 or S+,s-
    reference: False

    Allele
    genotype: GYPB*06.02
    defining_variants:
            4:144001250_T_C #all these, so cant be opposite
            4:144001249_C_A
            4:143999443_G_A #hom
            4:144001261_T_C
            4:144001262_A_C
            4:144001254_T_A
    weight_geno: 1000
    phenotype: MNS:3,-4,6 or S+,s-,He+
    reference: False

    """

    if not phased:
        return bg

    for allele_state in [AlleleState.NORMAL, AlleleState.CO]:
        if not proceed(bg, allele_state):
            continue
        to_remove = []
        for pair in bg.alleles[allele_state]:
            for variant in pair.allele1.defining_variants:
                if bg.variant_pool.get(variant) != Zygosity.HET:
                    continue
                phase = bg.variant_pool_phase[variant]
                for variant2 in pair.allele2.defining_variants:
                    if bg.variant_pool.get(variant2) != Zygosity.HET:
                        continue
                    phase2 = bg.variant_pool_phase[variant2]
                    if phase == phase2 and "|" in phase:
                        to_remove.append(pair)
        if to_remove:
            bg.remove_pairs(
                to_remove,
                "filter_if_all_HET_vars_on_same_side_and_phased",
                allele_state,
            )

    return bg


@apply_to_dict_values
def filter_on_in_relationship_if_HET_vars_on_dif_side_and_phased(
    bg: BloodGroup, phased: bool
) -> BloodGroup:
    """
    If an allele is HOM and it's 'in' every other properly phased allele
    AND the there's at least 1 of those on each side, it can't exist

    Sample: HG03774 BG Name: LU

    #Results:
    Genotypes count: 3
    Genotypes:
    LU*02.-13/LU*02.19
    LU*02/LU*02.-13 #not possible, HET SNPs on opposite sides so LU*02 'in'
    LU*02/LU*02.19  #not possible, HET SNPs on opposite sides
    Phenotypes (numeric): LU:-1,2,13 | LU:-1,2,13,18,19 | LU:-1,2,18,19
    Phenotypes (alphanumeric): Lu(a-b+),Au(a+b+) | Lu(a-b+),Au(a+b+),Lu13+ | Lu(a-b+),Lu13+

    #Data:

    Vars:
    19:44812188_ref : Homozygous
    19:44819059_C_T : Heterozygous
    19:44819705_A_T : Heterozygous
    19:44819634_C_T : Heterozygous
    19:44819487_A_G : Heterozygous
    Vars_phase:
    19:44812188_ref : 1/1
    19:44819059_C_T : 1|0
    19:44819705_A_T : 1|0
    19:44819634_C_T : 1|0
    19:44819487_A_G : 0|1
    Vars_phase_set:
    19:44812188_ref : .
    19:44819059_C_T : 43975436
    19:44819705_A_T : 43975436
    19:44819634_C_T : 43975436
    19:44819487_A_G : 43975436

    Raw:
    Allele
    genotype: LU*02
    defining_variants:
            19:44812188_ref
    weight_geno: 1000
    phenotype: LU:-1,2 or Lu(a-),Lu(b+)
    reference: True

    Allele
    genotype: LU*02.-13
    defining_variants:
            19:44812188_ref
            19:44819634_C_T #1|0
            19:44819705_A_T #1|0
            19:44819059_C_T #1|0
    weight_geno: 1000
    phenotype: LU:-1,2,-13 or Lu(a-),Lu(b+),Lu13-
    reference: False

    Allele
    genotype: LU*02.19
    defining_variants:
            19:44812188_ref
            19:44819487_A_G #0|1
    weight_geno: 1000
    phenotype: LU:-1,2,-18,19 or Lu(a-),Lu(b+),Au(a-),Au(b+)
    reference: False
    """

    # If an allele is HOM and it's 'in' every other properly phased allele
    # AND the there's at least 1 of those on each side, it can't exist

    if not phased:
        return bg
    for allele_state in [AlleleState.NORMAL, AlleleState.CO]:
        if not proceed(bg, allele_state):
            continue
        to_remove = []
        pairs_with_HET = []
        for pair in bg.alleles[allele_state]:
            if all_hom(bg.variant_pool, pair.allele1) or all_hom(
                bg.variant_pool, pair.allele2
            ):
                continue
            if not allele_phased(pair.allele1, bg.variant_pool_phase_set):
                continue  # TODO - next refactor this type of functionality
            # should move into a new PhasedAllele class
            if not allele_phased(pair.allele2, bg.variant_pool_phase_set):
                continue
            phase1 = find_phase(bg.variant_pool_phase, pair.allele1)
            phase2 = find_phase(bg.variant_pool_phase, pair.allele2)

            if phase1 == {None} or phase2 == {None}:
                continue
            if phase1 == {"unknown"} or phase2 == {"unknown"}:
                continue
            if len(phase1) == 1 and len(phase2) == 1:
                pairs_with_HET.append(pair)

        if pairs_with_HET:
            for pair_with_HET in pairs_with_HET:
                flattened_alleles = flatten_alleles([pair_with_HET])
                for pair in bg.alleles[allele_state]:
                    for allele in pair:
                        if all_hom(bg.variant_pool, allele) and all(
                            allele in flat_allele for flat_allele in flattened_alleles
                        ):
                            to_remove.append(pair)

        if to_remove:
            bg.remove_pairs(
                to_remove,
                "filter_on_in_relationship_if_HET_vars_on_dif_side_and_phased",
                allele_state,
            )

    return bg


@apply_to_dict_values
def filter_on_in_relationship_when_HOM_cant_be_on_one_side(
    bg: BloodGroup, phased: bool
) -> BloodGroup:
    """
    Example
    025-10-21 13:35:31.854 Sample: NA18913.vcf BG Name: ABO

    #Results:
    Genotypes count: 4
    Genotypes:
    ABO*O.01.83/ABO*O.01.83 #no - handled in
        filter_on_in_relationship_if_all_HOM_and_phased
    ABO*O.01.24/ABO*O.01.44
    ABO*O.01.24/ABO*O.01.83
    ABO*O.01.44/ABO*O.01.83 #no - because; ABO*O.01.83 and ABO*O.01.44 can be on the
    same side (neither is in the other) but ABO*O.01.83 and ABO*O.01.24 can't be on
    the same side and this ABO*O.01.44/ABO*O.01.83 means that they are on the same side
    because ABO*O.01.44 is 0|1 and ABO*O.01.24 is 1|0

    Phenotypes (numeric):
    Phenotypes (alphanumeric): O

    #Data:
    Vars:
    9:133257521_ref : Homozygous
    9:133257486_T_C : Homozygous
    9:133256074_G_A : Heterozygous
    9:133261367_C_A : Homozygous
    9:133259833_G_A : Homozygous
    9:133256028_C_T : Heterozygous
    9:133259834_C_T : Homozygous
    9:133256205_G_C : Heterozygous
    9:133255935_G_T : Heterozygous
    9:133255801_C_T : Heterozygous
    9:133255928_C_G : Heterozygous
    9:133255902_C_T : Heterozygous
    9:133256085_A_T : Heterozygous
    9:133255960_G_A : Heterozygous
    Vars_phase:
    9:133257521_ref : 1/1
    9:133257486_T_C : 1/1
    9:133256074_G_A : 1|0
    9:133261367_C_A : 1/1
    9:133259833_G_A : 1/1
    9:133256028_C_T : 1|0
    9:133259834_C_T : 1/1
    9:133256205_G_C : 1|0
    9:133255935_G_T : 1|0
    9:133255801_C_T : 1|0
    9:133255928_C_G : 1|0
    9:133255902_C_T : 0|1
    9:133256085_A_T : 0|1
    9:133255960_G_A : 0|1
    Vars_phase_set:  (all 133104364)

    Allele
    genotype: ABO*O.01.83 (in ABO*O.01.24)
    defining_variants:
            9:133257486_T_C 1/1
            9:133257521_ref 1/1
            9:133261367_C_A 1/1
            9:133259834_C_T 1/1
    weight_geno: 1000
    phenotype: . or O
    reference: False

    Allele
    genotype: ABO*O.01.24
    defining_variants:
            9:133255928_C_G 1|0
            9:133257521_ref 1/1
            9:133255935_G_T 1|0
            9:133256074_G_A 1|0
            9:133257486_T_C 1/1
            9:133255801_C_T 1|0
            9:133261367_C_A 1/1
            9:133259833_G_A 1/1
            9:133256028_C_T 1|0
            9:133259834_C_T 1/1
            9:133256205_G_C 1|0
    weight_geno: 1000
    phenotype: . or O
    reference: False

    Allele
    genotype: ABO*O.01.44
    defining_variants:
            9:133257521_ref 1/1
            9:133257486_T_C 1/1
            9:133255960_G_A 0|1
            9:133256085_A_T 0|1
            9:133255902_C_T 0|1
    weight_geno: 1000
    phenotype: . or O
    reference: False
    """

    if not phased:
        return bg
    for allele_state in [AlleleState.NORMAL, AlleleState.CO]:
        if not proceed(bg, allele_state):
            continue
        to_remove = []
        fully_phased_pairs = []
        for pair in bg.alleles[allele_state]:
            if not allele_phased(pair.allele1, bg.variant_pool_phase_set):
                continue  # TODO - next refactor this type of functionality
            # should move into a new PhasedAllele class
            if not allele_phased(pair.allele2, bg.variant_pool_phase_set):
                continue
            phase1 = find_phase(bg.variant_pool_phase, pair.allele1)
            phase2 = find_phase(bg.variant_pool_phase, pair.allele2)
            if phase1 == {None} or phase2 == {None}:
                continue
            if phase1 == {"unknown"} or phase2 == {"unknown"}:
                continue
            fully_phased_pairs.append(pair)
        if fully_phased_pairs:
            flattened_alleles = flatten_alleles(fully_phased_pairs)
            for pair in fully_phased_pairs:
                if all_hom(bg.variant_pool, pair.allele1) or all_hom(
                    bg.variant_pool, pair.allele2
                ):
                    if all_hom(bg.variant_pool, pair.allele1):
                        homs_partner_allele = pair.allele2
                        hom_allele = pair.allele1
                    else:
                        homs_partner_allele = pair.allele1
                        hom_allele = pair.allele2
                    phase_of_homs_partner = find_phase(
                        bg.variant_pool_phase, homs_partner_allele
                    )
                    for flat_allele in flattened_alleles:
                        if flat_allele in pair.alleles:
                            continue
                        phase_of_flat_allele = find_phase(
                            bg.variant_pool_phase, flat_allele
                        )
                        if (
                            phase_of_flat_allele != phase_of_homs_partner
                            and hom_allele in flat_allele
                        ):
                            to_remove.append(pair)
        if to_remove:
            bg.remove_pairs(
                to_remove,
                "filter_on_in_relationship_when_HOM_cant_be_on_one_side",
                allele_state,
            )

    return bg


def find_phase(variant_pool_phase: dict[str, str], allele: Allele) -> set[str | None]:
    """Extract phasing information for an allele's defining variants.

    Collects the phase assignments for all defining variants of an allele,
    excluding heterozygous unphased variants (0/1 or 1/0). This helps determine
    which haplotype (phase group) the allele belongs to.

    Args:
        variant_pool_phase: Dictionary mapping variant identifiers to their
            phase/genotype strings (e.g., "1|0", "0|1", "0/1").
        allele: Allele object containing the defining variants to check.

    Returns:
        Set of phase strings found for the allele's variants, excluding
        heterozygous unphased calls. May contain None for variants not in
        the phase pool.

    Example:
        >>> phase_pool = {"var1": "1|0", "var2": "0|1", "var3": "0/1"}
        >>> allele = Allele(defining_variants=["var1", "var2", "var3"])
        >>> find_phase(phase_pool, allele)
        {'1|0', '0|1'}
    """
    return set(
        [
            variant_pool_phase.get(variant)
            for variant in allele.defining_variants
            if variant_pool_phase.get(variant) not in ["0/1", "1/0"]
        ]
    )


def allele_phased(allele: Allele, phase_dict: dict[str, str]) -> bool:
    """Check if an allele's variants belong to a consistent phase set.

    Determines whether all defining variants of an allele are consistently phased
    by verifying they belong to the same phase set. Alleles with all homozygous
    variants (indicated by '.') are considered phased by definition.

    Args:
        allele: Allele object containing defining variants to check for phasing.
        phase_dict: Dictionary mapping variant identifiers to their phase set
            identifiers (numeric strings) or '.' for homozygous variants.

    Returns:
        True if all variants belong to the same phase set or are all homozygous,
        False if variants belong to different phase sets (indicating inconsistent
        phasing that prevents reliable allele determination).

    Example:
        >>> phase_dict = {"var1": "12345", "var2": "12345", "var3": "."}
        >>> allele = Allele(defining_variants=["var1", "var2"])
        >>> allele_phased(allele, phase_dict)
        True
    """
    phase_sets = set(
        [phase_dict.get(variant, "None") for variant in allele.defining_variants]
    )
    if phase_sets == set("."):
        return True  # all hom
    return len([phase_set for phase_set in phase_sets if phase_set.isdigit()]) == 1


@apply_to_dict_values
def filter_on_in_relationship_if_all_HOM_and_phased(
    bg: BloodGroup, phased: bool
) -> BloodGroup:
    """Filter out allele pairs inconsistent with variant phasing when all variants are phased.

    Removes allele pairs that violate phase relationships when all defining variants
    are either homozygous or phased heterozygous. This prevents invalid combinations
    like pairing alleles that should be on opposite haplotypes.

    Args:
        bg: BloodGroup object containing allele pairs and variant information.
        phased: Boolean indicating whether the sample variants are phased.

    Returns:
        The modified BloodGroup object with phase-inconsistent pairs removed,
        or the original object unchanged if phasing information is unavailable
        or variants are not fully phased.

    Example:
    2025-10-21 13:35:31.854 | DEBUG    | Sample: NA18913.vcf BG Name: ABO
    Genotypes: ABO*O.01.24/ABO*O.01.44
    #Results:
    Genotypes count: 4
    Genotypes:
    ABO*O.01.83/ABO*O.01.83 #no
    ABO*O.01.24/ABO*O.01.44
    ABO*O.01.24/ABO*O.01.83
    ABO*O.01.44/ABO*O.01.83 #no - but dif filter needed
    Phenotypes (numeric):
    Phenotypes (alphanumeric): O

    #Data:
    Vars:
    9:133257521_ref : Homozygous
    9:133257486_T_C : Homozygous
    9:133256074_G_A : Heterozygous
    9:133261367_C_A : Homozygous
    9:133259833_G_A : Homozygous
    9:133256028_C_T : Heterozygous
    9:133259834_C_T : Homozygous
    9:133256205_G_C : Heterozygous
    9:133255935_G_T : Heterozygous
    9:133255801_C_T : Heterozygous
    9:133255928_C_G : Heterozygous
    9:133255902_C_T : Heterozygous
    9:133256085_A_T : Heterozygous
    9:133255960_G_A : Heterozygous
    Vars_phase:
    9:133257521_ref : 1/1
    9:133257486_T_C : 1/1
    9:133256074_G_A : 1|0
    9:133261367_C_A : 1/1
    9:133259833_G_A : 1/1
    9:133256028_C_T : 1|0
    9:133259834_C_T : 1/1
    9:133256205_G_C : 1|0
    9:133255935_G_T : 1|0
    9:133255801_C_T : 1|0
    9:133255928_C_G : 1|0
    9:133255902_C_T : 0|1
    9:133256085_A_T : 0|1
    9:133255960_G_A : 0|1
    Vars_phase_set:  (all 133104364)
    """

    if not phased:
        return bg
    for allele_state in [AlleleState.NORMAL, AlleleState.CO]:
        if not proceed(bg, allele_state):
            continue
        if len(bg.alleles[allele_state]) < 2:
            continue
        to_remove = []
        flattened_alleles = flatten_alleles(bg.alleles[allele_state])
        for pair in bg.alleles[allele_state]:
            if (
                pair.allele1 == pair.allele2
                and all_hom(bg.variant_pool, pair.allele1)
                and any(
                    pair.allele1 in flat_allele for flat_allele in flattened_alleles
                )
            ):
                to_remove.append(pair)

        if to_remove:
            bg.remove_pairs(
                to_remove,
                "filter_on_in_relationship_if_all_HOM_and_phased",
                allele_state,
            )

    return bg


@apply_to_dict_values
def filter_pairs_by_phase(
    bg: BloodGroup, phased: bool, reference_alleles
) -> BloodGroup:
    """
    Filters out allele pairs where both alleles are in the same phase.

    This function is intended to remove allele pairs from a BloodGroup object when both
    alleles in a pair are in the same phase, indicating they are on the same
    chromosome and cannot be inherited together. The function operates under the
    following logic:

    - If `phased` is False, the function returns the BloodGroup object unchanged.
    - For each allele pair in `bg.alleles[AlleleState.NORMAL]`:
        - If the pair contains a reference allele, it is retained.
        - Extract the phase sets (`p1` and `p2`) for each allele in the pair.
        - If both alleles are homozygous (phase sets are {"."}), the pair is retained.
        - If the phase sets are identical, the pair is removed.
        - If the non-homozygous phase sets differ, the pair is retained.
        - If the non-homozygous phase sets are identical, the pair is removed.
        - If none of the above conditions are met, a ValueError is raised.

    - If all pairs are removed and there were pairs with phase information, new pairs
    are created by pairing each allele with the reference allele for the blood group
    type.

    Args:
        bg (BloodGroup): The BloodGroup object containing allele pairs.
        phased (bool): A flag indicating whether phase information is available.
        reference_alleles (dict): A dictionary mapping blood group types to reference
        alleles.

    Returns:
        BloodGroup: The updated BloodGroup object with inconsistent allele pairs removed.


    Example:
    ----------
    Suppose you have allele pairs where both alleles are on the same phase strand.
    This function will remove such pairs, ensuring that only valid allele combinations
    are retained. If all pairs are removed and phase information is present, it will
    create new pairs with the reference allele to represent possible allele
    combinations.


    Meant to remove pairs where both alleles are on the same strand ie
    to_remove: [[Allele(genotype='FUT2*01N.16',
                        defining_variants=frozenset({'19:48703728_G_A'}),
                        weight_geno=8,
                        weight_pheno=5,
                        reference=False,
                        sub_type='FUT2*01',
                        phase=0|1),
                 Allele(genotype='FUT2*01N.02',
                        defining_variants=frozenset({'19:48703417_G_A'}),
                        weight_geno=1,
                        weight_pheno=5,
                        reference=False,
                        sub_type='FUT2*01',
                        phase=0|1)]]

    dont remove if ref in pair
    if there is only 1 pair and they are phased then change to 2 pairs (or &) with ref

    The equality test needs both sides to be real phase, which is what carries_phase
    decides. Two values being equal only means the alleles are on one chromosome if the
    values locate a chromosome in the first place, and three that do not compare equal
    to each other all the time: 'unknown' for a variant the caller did not phase, an
    unphased genotype like '0/1' for a heterozygote in a partly phased file, and a
    homozygous genotype, which is on both chromosomes rather than one. The phase set
    check above does not catch them, because it counts phase sets rather than reading
    them, and a pair with nothing phased has one shared phase set of 'unknown'.

    Sample: HG00099 BG Name: FUT3, where both defining variants are heterozygous and
    neither was phased:

        Vars_phase:
        19:5844526_ref : unknown
        19:5844638_ref : unknown

    FUT3*01.04 needs 19:5844526_ref and FUT3*01N.03.01 needs 19:5844638_ref, so both
    alleles' phase read 'unknown', matched, and the pair was removed as though the two
    had been seen on the same chromosome. Removing it also took FUT3*01.04 out of the
    pool that cant_pair_with_ref_cuz_SNPs_must_be_on_other_side reads later, so the
    reference pair FUT3*01.01/FUT3*01N.03.01 - which leaves the heterozygous
    19:5844526_ref with nowhere to sit - was reported instead of being excluded. The
    call changed rather than narrowing: unphased says
    FUT3*01.04/FUT3*01N.03.01,FUT3*01.01/FUT3*01N.03.02 and phased said
    FUT3*01.01/FUT3*01N.03.01,FUT3*01.01/FUT3*01N.03.02.

    A pair holding the same allele in both slots is skipped for a different reason. It
    is not two alleles competing for one chromosome, it is one allele written the way
    pairs are written everywhere else, so the same-strand question does not apply and
    the equality test can only ever say yes. The all hom escape below covers the
    diploid form of this and nothing covered the single copy form.

    Sample: HG01873 BG Name: XK, a single copy region:

        Vars:
        X:37727623_C_G : c.509-13C>G : Hemizygous
        Vars_phase:
        X:37727623_C_G : c.509-13C>G : 1

    XK*N.33 sits on the one chromosome the locus has, which is what the phase '1' says,
    and the pair is XK*N.33/XK*N.33 because that is how a single copy region is carried
    through the filters - the second slot only becomes HAPLOID_SECOND_SLOT at reporting.
    Both sides read '1', matched, and the pair went. It was the last one, so the branch
    below paired the reference with allele1 and with allele2 - the same allele twice -
    which put XK*01/XK*N.33 in twice and reinstated the pair excluded_due_to_rank_ref
    had already excluded by name. Kx- became Kx+.
    """

    if not phased:
        return bg
    to_remove = []
    for pair in bg.alleles[AlleleState.NORMAL]:
        if pair.contains_reference:
            continue
        if pair.allele1 == pair.allele2:
            continue  # one allele in two slots, not two alleles on one chromosome
        p1_phases = set(_get_allele_phase_info(pair.allele1, bg.variant_pool_phase))
        p1_zygo = set(_get_allele_phase_info(pair.allele1, bg.variant_pool))
        p1_phase_sets = set(
            _get_allele_phase_info(pair.allele1, bg.variant_pool_phase_set)
        )
        p2_phases = set(_get_allele_phase_info(pair.allele2, bg.variant_pool_phase))
        p2_zygo = set(_get_allele_phase_info(pair.allele2, bg.variant_pool))
        p2_phase_sets = set(
            _get_allele_phase_info(pair.allele2, bg.variant_pool_phase_set)
        )

        phase_set = p1_phase_sets.union(p2_phase_sets)
        if len(phase_set) != 1:
            continue  # can't use phasing info

        if not all(carries_phase(phase) for phase in p1_phases | p2_phases):
            continue  # no value here says which chromosome anything is on

        if p1_zygo == {Zygosity.HOM} and p2_zygo == {Zygosity.HOM}:  # all hom
            continue
        elif p1_phases == p2_phases:
            to_remove.append(pair)
    if len(bg.alleles[AlleleState.NORMAL]) == len(to_remove):
        for pair in to_remove:
            bg.alleles[AlleleState.NORMAL].append(
                Pair(reference_alleles[bg.type], pair.allele1)
            )
            bg.alleles[AlleleState.NORMAL].append(
                Pair(reference_alleles[bg.type], pair.allele2)
            )
    if to_remove:
        bg.remove_pairs(to_remove, "filter_pairs_by_phase")

    return bg


@apply_to_dict_values
def impossible_alleles_phased(bg: BloodGroup, phased: bool) -> BloodGroup:
    """
    Filters alleles in a BloodGroup object based on phasing consistency and subsumption.

    When `phased` is True, this function attempts to simplify the set of possible
    alleles (`bg.alleles[AlleleState.]`) by looking for alleles that are impossibe
    due to being 'in' (defining variants are subset of) another allele:
        1. A1 is A2, all vars HOM (should be removed above - not phased dependant)
        2. A1 is A2, all vars HET and in same phase
        2. A1 is A2, vars are mix of HET and HOM, all HET are in same phase


    The function modifies `bg.alleles[AlleleState.]` in-place by removing the
    alleles deemed "impossible" under these phased conditions. Details of removed
    allele pairs due to subsumption are stored in
    `bg.filtered_out["allele_subsumed_by_other_phased"]`.

    Args:
        bg (BloodGroup): A BloodGroup object containing allele states, variant pool
            (with observed zygosities), and phasing information within Allele objects.
        phased (bool): A flag indicating whether phasing rules should be applied.
            If False, the function returns the BloodGroup object unmodified.

    Returns:
        BloodGroup: The modified BloodGroup object. The `alleles[AlleleState.]`
            list may be reduced, and `filtered_out` may be updated.


    Example (same phase mainly HET):
    GYPB*03/GYPB*03N.03 - in GYPB*03N.04
    GYPB*03/GYPB*03N.04 - is the only posibility as all vars in phase
    GYPB*03/GYPB*06.02 - in GYPB*03N.04
    Phenotypes (numeric): MNS:3,-4,5 | MNS:3,-4,6
    Phenotypes (alphanumeric): S+,s-,He+ | S+,s-,U+

    #Data:
    Vars: {
    '4:143999443_G_A': 'Homozygous',
    '4:144001261_T_C': 'Heterozygous',
    '4:144001254_T_A': 'Heterozygous',
    '4:144001250_T_C': 'Heterozygous',
    '4:144001249_C_A': 'Heterozygous',
    '4:144001262_A_C': 'Heterozygous',
    '4:143997535_C_A': 'Heterozygous'}

    Filtered:
    Allele
    genotype: GYPB*03
    defining_variants:
            4:143999443_G_A #hom
    weight_geno: 1000
    phenotype: MNS:3,-4 or S+,s-
    reference: False
    phases: ('.',)

    Allele
    genotype: GYPB*06.02
    defining_variants:
            4:144001249_C_A
            4:143999443_G_A #hom
            4:144001254_T_A
            4:144001262_A_C
            4:144001261_T_C
            4:144001250_T_C
    weight_geno: 1000
    phenotype: MNS:3,-4,6 or S+,s-,He+
    reference: False
    phases: ('143997535', '143997535', '143997535',
      '143997535', '143997535', '.')

    Allele
    genotype: GYPB*03N.03
    defining_variants:
            4:143997535_C_A
            4:143999443_G_A #hom
    weight_geno: 1000
    phenotype: MNS:-3,-4,5w or S-,s-,U+w
    reference: False
    phases: ('.', '143997535')

    Allele
    genotype: GYPB*03N.04
    defining_variants:
            4:143997535_C_A
            4:144001249_C_A
            4:143999443_G_A #hom
            4:144001254_T_A
            4:144001262_A_C
            4:144001261_T_C
            4:144001250_T_C
    weight_geno: 1000
    phenotype: MNS:-3,-4,5w or S-,s-,U+w
    reference: False
    phases: ('143997535', '143997535', '143997535',
      '143997535', '143997535', '143997535', '.')


    """

    if not phased:
        return bg

    for allele_state in [AlleleState.NORMAL, AlleleState.CO]:
        if not proceed(bg, allele_state):
            continue
        if len(bg.alleles[allele_state]) in [1, 0]:
            return bg
        # process alleles
        alleles = list(flatten_alleles(bg.alleles[allele_state]))

        # A variant at a single-copy locus is on the one chromosome that
        # locus has. That is not phase and it needs no phase set - the checks below
        # ask which of two chromosomes a variant sits on, and there is only one.
        # Recognised here rather than given a synthesised phase set upstream: a
        # fabricated set is a second value in the pool, which breaks the 'all variants
        # share one phase set' precondition modify_variant_phase_pool_if_large_indel
        # needs, silently stopping that repair.
        def all_hemizygous(allele) -> bool:
            zygos = [
                bg.variant_pool.get(variant) for variant in allele.defining_variants
            ]
            return bool(zygos) and all(z == Zygosity.HEM for z in zygos)

        alleles_with_variants_in_same_phase_set = [
            allele
            for allele in alleles
            if all_hemizygous(allele)
            or check_phase(bg.variant_pool_phase_set, allele, ".")
        ]

        alleles_with_variants_in_same_phase = [
            allele
            for allele in alleles_with_variants_in_same_phase_set
            if all_hemizygous(allele)
            or check_phase(bg.variant_pool_phase, allele, "1/1")
        ]

        # split by phase
        l1, l2, l3 = [], [], []  # 1|0, 0|1, or 1 (hemi)
        for allele in sorted(
            alleles_with_variants_in_same_phase,
            key=lambda allele: len(allele.defining_variants),
            reverse=True,
        ):
            phases = set([])
            for variant in allele.defining_variants:
                phase = bg.variant_pool_phase[variant]
                if phase != "1/1":
                    phases.add(bg.variant_pool_phase[variant])

            assert len(phases) == 1
            phase = phases.pop()
            if phase == "1|0":
                l1.append(allele)
            elif phase == "0|1":
                l2.append(allele)
            elif phase == "1":
                l3.append(allele)
            else:
                assert phase in "unknown" or "/" in phase
        # figure out what to remove

        alleles_to_remove = iterate_over_list(l1)
        alleles_to_remove += iterate_over_list(l2)
        alleles_to_remove += iterate_over_list(l3)
        to_remove = []
        for pair in bg.alleles[allele_state]:
            if pair.allele1 in alleles_to_remove or pair.allele2 in alleles_to_remove:
                to_remove.append(pair)
        assert len(to_remove) != len(bg.alleles[allele_state])

        if to_remove:
            bg.remove_pairs(to_remove, "filter_impossible_alleles_phased", allele_state)
        assert bg.alleles[allele_state]

    return bg


@apply_to_dict_values
def narrow_second_slot_candidates_by_phase(
    bg: BloodGroup, phased: bool
) -> BloodGroup:
    """Drop a candidate another candidate already accounts for.

    cant_name_second_slot_cuz_ref_impossible names one chromosome and refuses the other,
    reporting one genotype per candidate allele. Nothing then narrows that list, because
    it removed the pairs to build it and what it leaves behind is not a pair - so the
    phase based 'in' filters, which run later and work on pairs, never see it.

    Where the candidates stand in a subset relation and phase puts every heterozygous
    variant on one side, the subset under-describes the chromosome. Choosing it would
    leave the superset's extra variants on that same chromosome with nothing to carry
    them, which is not a reading of the data, it is a gap in one.

    HG01527 RHCE in the long read set: five heterozygous variants, every one '0|1' in
    phase set 25233074, and six candidates of which RHCE*01.20.04.02 holds all five and
    the other five are strict subsets of it. Reported as six genotypes whose phenotypes
    disagree, so the phenotype columns come out empty - the cell has no answer rather
    than an imprecise one.

    Deliberately narrow. It acts only when **every** heterozygous variant across every
    candidate shares one side, which is the case where the superset is unambiguously the
    whole story. Two sides means the candidates describe different chromosomes and the
    subset may be the right one; unphased means nothing is known and ambiguity is the
    correct output. Homozygous variants are not consulted at all - they are on both
    chromosomes and locate nothing.

    Args:
        bg (BloodGroup): A BloodGroup whose second slot was refused by name.
        phased (bool): The --phased flag. Without it there are no sides to read.

    Returns:
        BloodGroup: With the superseded candidates dropped from single_slot_genotypes and
        recorded in filtered_out under this filter's name.
    """
    if not phased or len(bg.single_slot_genotypes) < 2:
        return bg

    by_genotype = {
        allele.genotype: allele for allele in bg.alleles[AlleleState.RAW]
    }
    candidates: dict[str, Allele] = {}
    for rendered in bg.single_slot_genotypes:
        name = rendered.split("/")[0]
        allele = by_genotype.get(name)
        if allele is None:
            # A candidate whose allele is not to hand cannot be reasoned about, and
            # dropping the others on a partial view would be worse than not acting.
            return bg
        candidates[name] = allele

    sides = {
        bg.variant_pool_phase.get(variant)
        for allele in candidates.values()
        for variant in allele.defining_variants
        if bg.variant_pool.get(variant) == Zygosity.HET
    }
    if len(sides) != 1:
        return bg
    side = sides.pop()
    if side is None or not carries_phase(side):
        return bg

    superseded = [
        allele
        for name, allele in candidates.items()
        if any(
            allele.defining_variants < other.defining_variants
            for other_name, other in candidates.items()
            if other_name != name
        )
    ]
    if not superseded or len(superseded) == len(candidates):
        return bg

    dropped = {allele.genotype for allele in superseded}
    bg.filtered_out["narrow_second_slot_candidates_by_phase"].extend(superseded)
    bg.single_slot_genotypes = [
        rendered
        for rendered in bg.single_slot_genotypes
        if rendered.split("/")[0] not in dropped
    ]

    return bg


def carries_phase(value: str) -> bool:
    """Whether a phase pool value says anything about which chromosome a variant is on.

    Three ways a value says nothing, and only the first was recognised before:

    '.' marks a homozygote in the phase set pool. 'unknown' marks a variant the caller
    did not phase - data_procesing already treats the two alike (see the ['unknown', '.']
    test there) and this brings the phase filters into line.

    An unphased genotype - '0/1', '1/1', './.' - has no bar, so it does not say which
    chromosome anything is on. Every heterozygote in a partly phased file looks like
    '0/1', and they all look like each other, which is how two different alleles came to
    be read as identically phased.

    A homozygous genotype says nothing either, and it can be written with a bar: '1|1'
    and '0|0' are on both chromosomes. Only '1/1' was filtered out by name, so '1|1'
    survived as if it located something.

    SYNTHESISED_HOM_REF_GT is rejected by name for a different reason. The final branch
    has to accept a bare token, because a caller's phase set id is one - '25233074' has
    no separator to test. The sentinel is a bare token too, so it fell through there and
    was read as a phase set. Named rather than shape checked: nothing guarantees a phase
    set id is numeric, so testing for digits would risk refusing a real one.
    """
    if value in {".", "unknown", "", SYNTHESISED_HOM_REF_GT}:
        return False
    if "/" in value:
        return False
    if "|" in value:
        left, _, right = value.partition("|")
        return left != right
    return True


def check_phase(variant_pool: dict[str, str], current_allele: Allele, hom: str) -> bool:
    """Check if an allele's heterozygous variants belong to a single phase set.

    Verifies that all non-homozygous defining variants of an allele belong to the
    same phase set, which is necessary for determining if variants are on the same
    haplotype. Homozygous variants are excluded from this check as they don't
    provide phasing information.

    Args:
        variant_pool: Dictionary mapping variant identifiers to their phase set IDs.
        current_allele: Allele object containing defining variants to check.
        hom: String identifier representing homozygous variants (e.g., "." or "1/1")
            to exclude from phase checking.

    Returns:
        True if all heterozygous variants belong to exactly one phase set (indicating
        consistent phasing), False if variants belong to multiple phase sets or no
        heterozygous variants are found.
    """

    phase_sets = [
        phase
        for variant, phase in variant_pool.items()
        if variant in current_allele.defining_variants and phase != hom
    ]
    if not all(carries_phase(phase) for phase in phase_sets):
        return False

    return len(set(phase_sets)) == 1


def iterate_over_list(allele_list: list[Allele]) -> list[Allele]:
    """Identify alleles that are subsets of other alleles in the list.

    Finds all alleles whose defining variants are completely contained within
    another allele's defining variants. This helps identify redundant or less
    specific allele calls that should be removed in favor of more complete
    allele definitions.

    Args:
        allele_list: List of Allele objects to compare against each other.

    Returns:
        List of alleles that are subsets of at least one other allele in the
        input list. An allele may appear multiple times if it's a subset of
        multiple other alleles.

    Note:
        Uses the 'in' operator which should be defined in the Allele class
        to check subset relationships between alleles.
    """
    alleles_to_remove: list[Allele] = []
    for allele in allele_list:
        for allele2 in allele_list:
            if allele in allele2:
                alleles_to_remove.append(allele)
    return alleles_to_remove


@apply_to_dict_values
def rm_ref_if_2x_HET_phased(bg: BloodGroup, phased: bool) -> BloodGroup:
    """Remove reference allele pairs when phased heterozygous variants support non-reference pairs.

    Eliminates reference allele pairs that are artifacts of NoHomMultiVariantStrategy
    when all variants are heterozygous. If proper phasing information confirms that
    two non-reference alleles are on opposite haplotypes, the reference-containing
    pairs can be safely removed as they represent incomplete interpretations.

    Args:
        bg: BloodGroup object containing allele pairs and phasing information.
        phased: Boolean indicating whether the sample variants are phased.

    Returns:
        The modified BloodGroup object with reference-containing pairs removed if
        a phased non-reference pair exists, or the original object unchanged if
        not phased or no valid phased pairs exist.

    Example:
    2025-09-08 10:22:21.891 | DEBUG    | Sample: HG02737.vcf BG Name: ABCC4

    #Results:
    Genotypes count: 3
    Genotypes:
    ABCC4*01.02W/ABCC4*01.03W
    ABCC4*01/ABCC4*01.02W #remove
    ABCC4*01/ABCC4*01.03W #remove
    Phenotypes (numeric): PEL:1
    PEL:1w
    Phenotypes (alphanumeric): PEL+
    PEL+w

    #Data:
    Vars:
    13:95206781_C_A : Heterozygous
    13:95163161_C_T : Heterozygous
    Vars_phase:
    13:95206781_C_A : 1|0
    13:95163161_C_T : 0|1
    Vars_phase_set:
    13:95206781_C_A : 94972116
    13:95163161_C_T : 94972116
    Raw:
    Allele
    genotype: ABCC4*01.02W
    defining_variants:
            13:95206781_C_A
    weight_geno: 1000
    phenotype: PEL:1w or PEL+w
    reference: False

    Allele
    genotype: ABCC4*01.03W
    defining_variants:
            13:95163161_C_T
    weight_geno: 1000
    phenotype: PEL:1w or PEL+w
    reference: False
    """

    if not phased:
        return bg
    to_remove = []
    phased_ref_free_pair_exists = False
    same_phase_set = partial(check_phase, bg.variant_pool_phase_set)
    same_phase = partial(check_phase, bg.variant_pool_phase)
    for pair in bg.alleles[AlleleState.NORMAL]:
        if pair.allele1.reference or pair.allele2.reference:
            to_remove.append(pair)
            continue
        if possible_to_use_phase(same_phase_set, same_phase, pair):
            phase1 = allele_phase(bg.variant_pool_phase, pair.allele1)
            phase2 = allele_phase(bg.variant_pool_phase, pair.allele2)
            assert phase1 != phase2
            phased_ref_free_pair_exists = True
    if to_remove and phased_ref_free_pair_exists:
        bg.remove_pairs(to_remove, "rm_ref_if_2x_HET_phased")

    return bg


def allele_phase(variant_pool, allele):
    """Extract the set of phase assignments for an allele's defining variants.

    Retrieves all phase information (e.g., "1|0", "0|1") for variants that define
    the given allele. This allows determination of which haplotype(s) the allele's
    variants belong to.

    Args:
        variant_pool: Dictionary mapping variant identifiers to their phase strings.
        allele: Allele object containing the defining variants to check.

    Returns:
        Set of phase strings for the allele's defining variants. An allele with
        variants on a single haplotype should return a single phase value, while
        variants on different haplotypes or unphased variants may return multiple
        values.

    Example:
        >>> phase_pool = {"var1": "1|0", "var2": "1|0", "var3": "0|1"}
        >>> allele = Allele(defining_variants=["var1", "var2"])
        >>> allele_phase(phase_pool, allele)
        {'1|0'}
    """
    return set(
        [
            phase
            for variant, phase in variant_pool.items()
            if variant in allele.defining_variants
        ]
    )


def possible_to_use_phase(same_phase_set: Callable, same_phase: Callable, pair: Pair):
    """Check if a pair of alleles can be phased using available phasing information.

    Determines whether both alleles in a pair have consistent phase set assignments
    and phase information, excluding homozygous variants. Both alleles must have
    their heterozygous variants in a single phase set and a single phase to be
    considered phaseable.

    Args:
        same_phase_set: Callable that checks if an allele's variants belong to one
            phase set, excluding a specified homozygous indicator.
        same_phase: Callable that checks if an allele's variants have consistent
            phase values, excluding a specified homozygous genotype.
        pair: Pair object containing two alleles to check for phasing consistency.

    Returns:
        True if both alleles have consistent phase set and phase assignments for
        their heterozygous variants, False otherwise.

    Note:
        Uses "." to represent homozygous variants in phase sets and "1/1" to
        represent homozygous genotypes in phase values.
    """
    return (same_phase_set(pair.allele1, ".") and same_phase(pair.allele2, "1/1")) and (
        same_phase_set(pair.allele2, ".") and same_phase(pair.allele2, "1/1")
    )


@apply_to_dict_values
def low_weight_hom(bg: BloodGroup, phased: bool) -> BloodGroup:
    """Remove allele pairs when phasing identifies a higher-weight pair on opposite
    haplotypes.

    Handles cases where a homozygous variant exists but isn't in the top 2 ranked
    chunk,
    requiring SomeHomMultiVariantStrategy to defer filtering. When phasing information
    shows the two highest-weight alleles are on opposite haplotypes, other pairs can
    be eliminated. Selects the pair with the lowest combined weight that has alleles
    on different phases.

    Case where there's a hom but it isn't in the top 2 ranked chunk,
    so SomeHomMultiVariantStrategy has to let it pass to here.

    Args:
        bg: BloodGroup object containing allele pairs and phasing information.
        phased: Boolean indicating whether the sample variants are phased.

    Returns:
        The modified BloodGroup object retaining only the lowest-weight phased pair
        if multiple phased pairs with different weights exist, or the original object
        unchanged if not phased, only one phased pair exists, or all pairs have equal
        weight.

    Example:
        in this example the 2 highest weights are in in opposite phase.
    no example where they're in phase, yet
    2025-09-17 08:08:58.947 | DEBUG | Sample: HG03437.vcf BG Name: FUT2

    #Results:
    Genotypes count: 4
    Genotypes:
    FUT2*01N.02/FUT2*01N.04 #only this is possible as in opposite phase
                            #and are the 2 biggest weights
    FUT2*01N.02/FUT2*01N.14
    FUT2*01N.02/FUT2*01N.16

    Phenotypes (numeric):
    Phenotypes (alphanumeric): Se-

    #Data:
    Vars:
    19:48703417_G_A : Heterozygous
    19:48703560_C_T : Heterozygous
    19:48703939_C_T : Heterozygous
    19:48703728_G_A : Homozygous
    Vars_phase:
    19:48703417_G_A : 0|1
    19:48703560_C_T : 1|0
    19:48703939_C_T : 1|0
    19:48703728_G_A : 1/1
    Vars_phase_set:
    19:48703417_G_A : 48646783
    19:48703560_C_T : 48646783
    19:48703939_C_T : 48646783
    19:48703728_G_A : .
    Raw:
    Allele
    genotype: FUT2*01N.02
    defining_variants:
            19:48703417_G_A
    weight_geno: 1
    phenotype: . or Se-
    reference: False

    Allele
    genotype: FUT2*01N.04
    defining_variants:
            19:48703560_C_T
    weight_geno: 2
    phenotype: . or Se-
    reference: False

    Allele
    genotype: FUT2*01N.14
    defining_variants:
            19:48703939_C_T
    weight_geno: 8
    phenotype: . or Se-
    reference: False

    Allele
    genotype: FUT2*01N.16
    defining_variants:
            19:48703728_G_A
    weight_geno: 8
    phenotype: . or Se-
    reference: False
    """

    if not phased:
        return bg
    same_phase_set = partial(check_phase, bg.variant_pool_phase_set)
    same_phase = partial(check_phase, bg.variant_pool_phase)
    # store as list of tuples: (weight, pair)
    pairs: list[tuple[float, Pair]] = []
    for pair in bg.alleles[AlleleState.NORMAL]:
        if possible_to_use_phase(same_phase_set, same_phase, pair):
            phase1 = allele_phase(bg.variant_pool_phase, pair.allele1)
            phase2 = allele_phase(bg.variant_pool_phase, pair.allele2)
            assert phase1 != phase2
            pairs.append((pair.allele1.weight_geno + pair.allele2.weight_geno, pair))
    if not pairs:
        return bg
    if len(pairs) == 1:
        return bg
    weights = set([pair_tup[0] for pair_tup in pairs])
    if len(weights) == 1:
        return bg
    # select the lowest-weighted pair
    best_weight, best_pair = min(pairs, key=lambda x: x[0])

    to_remove = [pair for pair in bg.alleles[AlleleState.NORMAL] if pair != best_pair]
    if to_remove:
        bg.remove_pairs(to_remove, "low_weight_hom")

    return bg


def locus_has_a_copy_number(variant: str, variant_pool: dict[str, str]) -> bool:
    """Whether the sample has a copy number at this variant's position.

    A token is absent from the pool for two opposite reasons, and only one of them is
    evidence. The alternate at the locus being homozygous leaves no reference copy for a
    '_ref' token to sit on, which is a contradiction. The caller never looking leaves the
    same gap and contradicts nothing.

    Telling them apart is what the rest of the locus is for: a sibling token carrying a
    copy number means the position was measured. Zygosity.NO_DATA is the one state that
    does not count, which is the same rule variant_pool_numeric applies when it omits it -
    'not measured' has no copy number. Zygosity.NO_COPIES does count: zero copies is a
    measurement of absence rather than an absence of measurement, and a locus with no
    copies has no reference copy either.

    Args:
        variant (str): A defining variant of the reference allele, 'chrom:pos_REF_ALT'
        or 'chrom:pos_ref'.
        variant_pool (dict[str, str]): The blood group's pool, variant to zygosity.

    Returns:
        bool: True if some token at that position carries a copy number.

    Example:
        HG00109 HPA3 on an array. The probe did not call 17:42453065, so HPA3*02 goes to
        no_call_at_defining_variant and 17:42453065_ref is not in the pool either:

        17:42453065_A_C : (None) : No_data

        The only token at the locus is NO_DATA, so this returns False and HPA3*01/HPA3*01
        stands. Reading that absence as impossibility would report Undetermined for a
        sample nobody measured - 547 of them on this input.

        NA19679 RHCE, short read. RHCE*01 needs 1:25408711_ref, which is absent:

        1:25408711_G_A : c.307C>T : Homozygous

        The alternate is there with a copy number, so this returns True.
    """
    chrom, rest = variant.split(":", 1)
    prefix = f"{chrom}:{rest.split('_', 1)[0]}_"

    return any(
        token.startswith(prefix) and zygosity != Zygosity.NO_DATA
        for token, zygosity in variant_pool.items()
    )


@apply_to_dict_values
def no_defining_variant(bg: BloodGroup) -> BloodGroup:
    """Remove allele pairs where reference allele variants are contradicted by the pool.

    Eliminates pairs containing reference alleles that have defining variants not
    present in the sample's variant pool, where the pool has something to say about
    the locus. This occurs when a reference variant is
    impossible because the alternate allele is homozygous. Skips alleles defined
    only by absence markers (variants ending in '.') and specific known insertions.

    Runs in both arms. Nothing in here reads phase - the test is against bg.variant_pool
    and nothing else - so gating it on --phased meant the same sample got the pair removed
    with the flag and reported without it, a difference decided by a flag that is about
    something else. That is the same reason cant_name_second_slot_cuz_ref_impossible runs
    in both arms. Measured over nine datasets: 65 pairs in 59 cells, all of them narrowings
    and none of them on the array.

    Absence alone is not enough, and locus_has_a_copy_number is the gate. This filter walks
    both alleles of a pair, so it reaches reference/reference pairs too, and on an array a
    locus with no call makes every reference pair look impossible - 547 cells on the array
    input, against 281 pairs where a real call contradicts the reference. The two
    populations separate perfectly on whether the locus carries a copy number.

    The ABO c.261delG insertion is one of those, and for a different reason from the rest:
    the database treats the deletion as the reference sequence, so ABO*A1.01 is a reference
    allele defined by an *alternate*. Its absence from the pool is then the ordinary state
    of a group O sample rather than a contradiction, and removing the pair over it would
    make ABO uncallable for anyone who is not group A or B. Both builds are exempted from
    ABO_DELG_VARIANTS - the GRCh37 form used to be written here without its chromosome
    prefix, so it matched no token and the exemption only ever worked on GRCh38.

    An empty variant pool is not that case and is skipped. Absence of evidence is not
    evidence that the reference variant is impossible - the pool is empty here because
    every allele the blood group had was removed by a filter, most often because the
    caller doubted the variants they were built from. Falling back to the reference allele
    when nothing is left is the convention this tool follows, and it is what a lab
    scientist does by hand; removing the pair instead turns 'we found nothing' into 'we
    cannot say', which is a different and stronger claim. Measured on HG02308 KN and
    HG03673 RHD, both of which reported no call where the reference allele was the answer.

    Note this filter only ever sees a blood group that had alleles - one with no variants
    at all never enters the pipeline and gets its reference genotype from add_refs - so an
    empty pool here always means something was taken away.

    Args:
        bg: BloodGroup object containing allele pairs and variant pool information.

    Returns:
        The modified BloodGroup object with impossible reference allele pairs removed,
        or the original object unchanged if no invalid pairs were found.

    Example
    need to rm ref as 1:25390874_ref not possible
    as 1:25390874_C_G: Homozygous

    bg.variant_pool_phase: {'1:25390874_C_G': '1/1',
                            '1:25408711_G_A': '0|1',
                            '1:25408711_ref': '1|0',
                            '1:25408815_T_C': '0|1',
                            '1:25408817_T_C': '0|1',
                            '1:25408840_G_T': '0|1',
                            '1:25408868_G_A': '0|1',
                            '1:25420739_G_C': '1|0',
                            '1:25420739_ref': '0|1'}
    bg.variant_pool: {'1:25390874_C_G': 'Homozygous',
                      '1:25408711_G_A': 'Heterozygous',
                      '1:25408711_ref': 'Heterozygous',
                      '1:25408815_T_C': 'Heterozygous',
                      '1:25408817_T_C': 'Heterozygous',
                      '1:25408840_G_T': 'Heterozygous',
                      '1:25408868_G_A': 'Heterozygous',
                      '1:25420739_G_C': 'Heterozygous',
                      '1:25420739_ref': 'Heterozygous'}
    ic| allele: Allele
                genotype: RHCE*04
                defining_variants:
                        1:25390874_C_G
                        1:25408815_T_C
                        1:25408868_G_A
                        1:25408711_G_A
                        1:25420739_ref
                        1:25408840_G_T
                        1:25408817_T_C
                weight_geno: 1000
                phenotype: RH:2,3,-4,-5,22 or C+,E+,c-,e-,CE+
                reference: False
    ic| allele: Allele
                genotype: RHCE*01
                defining_variants:
                        1:25420739_G_C
                        1:25408711_ref
                        1:25390874_ref #not possible as 25390874_C_G HOM
                weight_geno: 1000
                phenotype: RH:-2,-3,4,5,6 or C-,E-,c+,e+,f+
                reference: True
    """

    if not bg.variant_pool:
        return bg
    to_remove = []

    for pair in bg.alleles[AlleleState.NORMAL]:
        for allele in pair.alleles:
            if not allele.reference:
                continue
            if all(variant.endswith(".") for variant in allele.defining_variants):
                continue
            contradicted = [
                variant
                for variant in allele.defining_variants
                if variant not in ABO_DELG_VARIANTS
                and variant not in bg.variant_pool
                and locus_has_a_copy_number(variant, bg.variant_pool)
            ]
            if contradicted:
                to_remove.append(pair)
                break
    if to_remove:
        bg.remove_pairs(to_remove, "no_defining_variant")

    return bg


@apply_to_dict_values
def ref_not_phased(bg: BloodGroup, phased: bool) -> BloodGroup:
    """Remove allele pairs containing reference alleles previously identified as unphased.

    Eliminates pairs where the reference allele was already filtered out for being
    unphased but was subsequently re-added during allele pair generation. This
    prevents inconsistent phasing information from producing invalid genotype calls.

    Args:
        bg: BloodGroup object containing allele pairs and filtering history.
        phased: Boolean indicating whether the sample variants are phased.

    Returns:
        The modified BloodGroup object with pairs containing previously filtered
        unphased reference alleles removed, or the original object unchanged if
        not phased or no such pairs exist.

    Example:
    need to rm ref due to not being phased
    ie
    2025-11-12 09:09:21.312 | WARNING  |
    rbceq2.core_logic.alleles:remove_pairs:350 - all pairs removed!:
      HG00128.vcf RHCE filter_if_all_HET_vars_on_same_side_and_phased
    #Results:
    Genotypes count: 1
    Genotypes:
    Phenotypes (numeric):
    Phenotypes (alphanumeric):

    #Data:
    Vars:
    1:25390874_ref : Heterozygous
    1:25420739_G_C : Heterozygous
    1:25408711_ref : Heterozygous
    1:25420739_ref : Heterozygous
    1:25408711_G_A : Heterozygous
    1:25390874_C_G : Heterozygous
    Vars_phase:
    1:25390874_ref : 0|1
    1:25420739_G_C : 1|0
    1:25408711_ref : 1|0
    1:25420739_ref : 0|1
    1:25408711_G_A : 0|1
    1:25390874_C_G : 1|0
    Vars_phase_set:
    1:25390874_ref : 25214110
    1:25420739_G_C : 25214110
    1:25408711_ref : 25214110
    1:25420739_ref : 25214110
    1:25408711_G_A : 25214110
    1:25390874_C_G : 25214110

    Raw:
    Allele
    genotype: RHCE*01
    defining_variants:
            1:25390874_ref 0|1
            1:25408711_ref 1|0
            1:25420739_G_C 1|0
    weight_geno: 1000
    phenotype: RH:-2,-3,4,5,6 or C-,E-,c+,e+,f+
    reference: True

    Allele
    genotype: RHCE*03
    defining_variants:
            1:25408711_ref 1|0
            1:25420739_G_C 1|0
            1:25390874_C_G 1|0
    weight_geno: 1000
    phenotype: RH:-2,3,4,-5,27 or C-,E+,c+,e-,cE+
    reference: False
    """

    if not phased:
        return bg
    to_remove = []

    for pair in bg.alleles[AlleleState.NORMAL]:
        for allele in pair.alleles:
            if not allele.reference:
                continue
            if allele in bg.filtered_out["remove_unphased"]:
                to_remove.append(pair)
                # refs get added back in even if they've been previoulsy
                # identified as unphased, its a bit clunky ...
    if to_remove:
        bg.remove_pairs(to_remove, "ref_not_phased")

    return bg


@apply_to_dict_values
def cant_be_hom_ref_due_to_HET_SNP(bg: BloodGroup, phased: bool) -> BloodGroup:
    """Remove homozygous reference pairs when any defining variant is heterozygous.

    Eliminates allele pairs where both alleles are reference but at least one
    defining variant is heterozygous in the sample. A homozygous reference genotype
    is incompatible with any heterozygous variants, as heterozygosity indicates
    the presence of at least one alternate allele.

    Args:
        bg: BloodGroup object containing allele pairs and variant zygosity information.
        phased: Boolean indicating whether the sample variants are phased.

    Returns:
        The modified BloodGroup object with impossible homozygous reference pairs
        removed, or the original object unchanged if not phased or no invalid
        pairs exist.

    Example:
    2025-11-12 13:47:12.043 | DEBUG    | Sample: HG00365.vcf BG Name: RHCE

    #Results:
    Genotypes count: 1
    Genotypes:
    Phenotypes (numeric):
    Phenotypes (alphanumeric):

    #Data:
    Vars:
    1:25408711_ref : Heterozygous
    1:25390874_ref : Homozygous
    1:25420739_G_C : Heterozygous
    1:25420739_ref : Heterozygous
    1:25408711_G_A : Heterozygous
    Vars_phase:
    1:25408711_ref : 0|1
    1:25390874_ref : 1/1
    1:25420739_G_C : 0|1
    1:25420739_ref : 1|0
    1:25408711_G_A : 1|0
    Vars_phase_set:
    1:25408711_ref : 25211850
    1:25390874_ref : .
    1:25420739_G_C : 25211850
    1:25420739_ref : 25211850
    1:25408711_G_A : 25211850

    Raw:
    Allele
    genotype: RHCE*01
    defining_variants:
            1:25408711_ref : 0|1
            1:25420739_G_C : 0|1
            1:25390874_ref : 1/1
    weight_geno: 1000
    phenotype: RH:-2,-3,4,5,6 or C-,E-,c+,e+,f+
    reference: True
    """

    if not phased:
        return bg
    to_remove = []

    for pair in bg.alleles[AlleleState.NORMAL]:
        if pair.all_reference and any(
            bg.variant_pool.get(variant) == Zygosity.HET
            for variant in pair.allele1.defining_variants
        ):
            to_remove.append(pair)
    if to_remove:
        # reverts_to_reference is False because the pair being removed *is* the
        # reference pair, and it is removed precisely because the sample cannot be
        # homozygous reference. The default warning would promise a revert to the one
        # answer this filter has just ruled out. Same reasoning as
        # cant_name_second_slot_cuz_ref_impossible, which the parameter's own docstring
        # cites.
        bg.remove_pairs(
            to_remove, "cant_be_hom_ref_due_to_HET_SNP", reverts_to_reference=False
        )

    return bg


@apply_to_dict_values
def cant_name_second_slot_cuz_hom_ref_impossible(
    bg: BloodGroup, phased: bool
) -> BloodGroup:
    """Name the reference slot when the only pair was hom reference and is impossible.

    The mirror of cant_name_second_slot_cuz_ref_impossible, and the other half of the
    same rule: one slot identified and the other matching no allele is written
    'X/Undetermined'. That filter covers the shape where the *reference* is the slot
    that cannot be named. This one covers the shape where the reference is the slot the
    data settles and the partner is what cannot be named.

    The shape: cant_be_hom_ref_due_to_HET_SNP has just removed the only pair the blood
    group had, because it was reference/reference and a defining variant is
    heterozygous. That removal is right - the sample is not homozygous reference. But a
    heterozygote is one chromosome carrying the reference base and one carrying the
    alternate, so removing the pair throws away a slot the data has settled. The blood
    group ends up 'Undetermined/Undetermined' when one of the two is known.

    Phased only, and the phase is what makes it safe rather than a guess. Every
    defining variant of the reference that says which chromosome it is on has to say
    the *same* chromosome. Two heterozygous defining variants on opposite sides means
    neither chromosome carries the whole reference allele, and then declining both
    slots is the honest answer - so that case is left alone. A partly phased file where
    the heterozygote is written without a bar says nothing about sides either, and is
    also left alone.

    Nothing is excluded here, so nothing is recorded in filtered_out: the pair was
    already removed and recorded by cant_be_hom_ref_due_to_HET_SNP, whose name this
    reads to find it. This filter only names what that one left unnamed.

    Args:
        bg (BloodGroup): The BloodGroup object, after cant_be_hom_ref_due_to_HET_SNP.
        phased (bool): Whether the sample's variants are phased.

    Returns:
        BloodGroup: The BloodGroup with single_slot_genotypes set to one
        'reference/Undetermined' string, or unchanged if any gate is not met.

    Example:
        NA18571 RHCE, phased. FILTER drops RHCE*02 for two LowQual defining variants
        and remove_unphased drops RHCE*01.01 and RHCE*01.36, leaving only the
        reference:

        bg.variant_pool_phase: {'1:25390874_ref': '1/1',
                                '1:25408711_G_A': '0|1',
                                '1:25408711_ref': '1|0',
                                '1:25420739_G_C': '1|0',
                                '1:25420739_ref': '0|1'}

        RHCE*01 is 1:25390874_ref, 1:25408711_ref and 1:25420739_G_C. The two that
        carry phase are both '1|0', so the left chromosome is RHCE*01 outright and the
        right carries 25408711_G_A and 25420739_ref - the RHCE*02 signature minus the
        variants FILTER discarded, which is why nothing can name it.

        single_slot_genotypes -> ['RHCE*01/Undetermined']

        Was 'Undetermined/Undetermined'.
    """
    if not phased:
        return bg
    # .get, not [], because the co-existing stages have not created the key yet. Kept
    # so that moving this filter later cannot silently override a Knops result.
    if bg.alleles.get(AlleleState.CO) is not None:
        return bg
    # Only where nothing survived. If any pair is left the blood group has an answer
    # and this must not touch it.
    if bg.alleles[AlleleState.NORMAL] or bg.single_slot_genotypes:
        return bg

    removed = bg.filtered_out.get("cant_be_hom_ref_due_to_HET_SNP", [])
    references = [
        pair.allele1
        for pair in removed
        if isinstance(pair, Pair) and pair.all_reference
    ]
    # One reference allele, or there is nothing unambiguous to name.
    if len({reference.genotype for reference in references}) != 1:
        return bg
    reference = references[0]

    sides = {
        bg.variant_pool_phase.get(variant)
        for variant in reference.defining_variants
    }
    sides = {side for side in sides if side is not None and carries_phase(side)}
    if len(sides) != 1:
        return bg

    bg.single_slot_genotypes = [f"{reference.genotype}/{UNDETERMINED_SLOT}"]

    return bg
