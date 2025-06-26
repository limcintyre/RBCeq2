from __future__ import annotations

import operator
from collections import defaultdict
from functools import partial
from rbceq2.core_logic.alleles import Allele, BloodGroup, Pair
from rbceq2.core_logic.constants import LOW_WEIGHT, AlleleState
from rbceq2.core_logic.utils import (
    Zygosity,
    apply_to_dict_values,
    check_available_variants,
    one_HET_variant,
)
from icecream import ic


def flatten_alleles(pairs: list[Pair]) -> set[Allele]:
    """Flatten the pairs into a set of alleles.

    Args:
        pairs (list[Pair]): A list of Pair objects, where each Pair is an
            iterable containing Allele objects.

    Returns:
        set[Allele]: A set containing all Allele objects from the given pairs.
    """
    return {allele for pair in pairs for allele in pair}


def split_pair_by_ref(pair: Pair) -> tuple[Allele, Allele]:
    """Split a pair of alleles into reference and non-reference.

    This function assumes that exactly one of the alleles in the pair is marked as a
    reference.

    Args:
        pair (Pair): The pair of alleles to split.

    Returns:
        tuple[Allele, Allele]: A tuple where the first element is the reference allele
                                and the second is the non-reference allele.

    Raises:
        ValueError: If both or neither alleles are marked as reference.
    """

    if pair.allele1.reference and not pair.allele2.reference:
        ref, allele = pair.allele1, pair.allele2
    elif not pair.allele1.reference and pair.allele2.reference:
        ref, allele = pair.allele2, pair.allele1
    else:
        raise ValueError("Both ref")

    return ref, allele


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
    """
    if bg.alleles[AlleleState.CO] is not None:
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

    This filtering ensures that no combination smaller than either 2/2 or 1/3 (in this
    case) exists. For example, KN*01.-13+KN*01.06/KN*01.07 is a subset of
    KN*01.-13+KN*01.06+KN*01.12/KN*01.07 and therefore cannot exist.

    Args:
        bg (BloodGroup): A BloodGroup object containing allele pairs in various states.

    Returns:
        BloodGroup: The updated BloodGroup after filtering co-existing allele pairs.
    """

    def compare_to_all(pair_compare, comparables):
        """Compare a pair against a list of pairs to check if it is a subset.

        Args:
            pair_compare (tuple): A pair representing the allele combination to compare.
            comparables (list[tuple]): A list of allele pairs to compare against.

        Returns:
            bool: True if the pair is a subset of any other pair, False otherwise.
        """

        def number_of_alleles(pair_to_compare):
            return sum(len(bit) for bit in pair_to_compare)

        def is_subset(allele_to_compare):
            return all(
                co_allele in other_pair[0] or co_allele in other_pair[1]
                for co_allele in allele_to_compare
            )

        for other_pair in comparables:
            if pair_compare == other_pair:
                continue
            no_alleles_in_mushed_pair = number_of_alleles(pair_compare)
            no_alleles_in_other_mushed_pair = number_of_alleles(other_pair)
            flat_alleles = flatten_alleles(pair_compare)
            flat_other_alleles = flatten_alleles(other_pair)

            if (
                no_alleles_in_mushed_pair == no_alleles_in_other_mushed_pair
                and flat_alleles == flat_other_alleles
            ):
                continue
                # don't remove due to het permutations ie
                # KN*01.10', 'KN*01.-05+KN*01.07
                # KN*01.-05', 'KN*01.07+KN*01.10

            if is_subset(pair_compare[0]) and is_subset(pair_compare[1]):
                different_alleles = [
                    allele
                    for allele in bg.alleles[AlleleState.RAW]
                    if allele.genotype in flat_other_alleles.difference(flat_alleles)
                ]
                alleles_with_1_HET_var = [
                    allele
                    for allele in different_alleles
                    if [
                        bg.variant_pool[variant] for variant in allele.defining_variants
                    ].count(Zygosity.HET)
                    == 1
                ]
                if alleles_with_1_HET_var and not pair.same_subtype:
                    # this is potentially still undercooked -
                    # what about a 1 HET that could be on the other side without
                    # defining anything? is that covered elsewhere?
                    # no examples are coming up
                    continue

                return True
        return False

    if bg.alleles[AlleleState.CO] is not None:
        to_remove = []
        all_comparable = parse_bio_info2(bg.alleles[AlleleState.CO])
        all_comparable_without_ref = parse_bio_info2(
            [pair for pair in bg.alleles[AlleleState.CO] if not pair.contains_reference]
        )
        for pair in bg.alleles[AlleleState.CO]:
            if pair.contains_reference:
                if compare_to_all(pair.comparable, all_comparable):
                    to_remove.append(pair)
            elif compare_to_all(pair.comparable, all_comparable_without_ref):
                to_remove.append(pair)
        bg.remove_pairs(to_remove, "filter_co_existing_subsets", AlleleState.CO)

    return bg


@apply_to_dict_values
def filter_pairs_on_antithetical_zygosity(
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
    ----------
    013Kenya
    hashed ones not possible because loci that defines the antithetical antigens is HET,
    meaning there will def be 1 FY*01* and 1 FY*02*
    (and in this case FY*02/FY*01 not possible as a modifying SNP exists)
    'FY*01N.01/FY*02'
    #'FY*01N.01/FY*01' not possible
    'FY*02N.01/FY*01'
    #'FY*02/FY*01'
    """

    to_remove = []
    var_pool_no_chrom = {k.split(":")[1]: v for k, v in bg.variant_pool.items()}
    if bg.type in antitheticals:
        # ic(1111,antitheticals[bg.type]) RHCe issue
        anti1, anti2 = antitheticals[bg.type]
        if var_pool_no_chrom.get(anti1) == Zygosity.HOM:
            return bg
        if var_pool_no_chrom.get(anti2) == Zygosity.HOM:
            return bg
        flattened_sub_types = {
            allele.sub_type
            for allele in flatten_alleles(bg.alleles[AlleleState.NORMAL])
        }
        if len(flattened_sub_types) > 1:  # doesn't work for MNS or others with multiple
            # subtypes that aren't associted with an antitetical
            for pair in bg.alleles[AlleleState.NORMAL]:
                flat_sub_types = {allele.sub_type for allele in pair}
                if flat_sub_types == flattened_sub_types:
                    continue
                else:
                    to_remove.append(pair)
    # if to_remove:
    #     ic(1111, var_pool_no_chrom,to_remove, antitheticals['GYPB'], bg.variant_pool, anti1, anti2, bg.variant_pool.get(anti1), Zygosity.HOM)
    bg.remove_pairs(to_remove, "filter_pairs_on_antithetical_zygosity")

    return bg


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

    """

    to_remove = []
    if bg.alleles[AlleleState.CO] is not None:
        if bg.type in antitheticals:
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

        bg.remove_pairs(
            to_remove, "filter_co_pairs_on_antithetical_zygosity", AlleleState.CO
        )

    return bg


@apply_to_dict_values
def filter_pairs_on_antithetical_modyfying_SNP(
    bg: BloodGroup, antitheticals: dict[str, list[str]]
) -> BloodGroup:
    """Process genetic data to identify alleles and genotypes.

    Parameters:
    ----------
    res (dict[str, list[Allele]]):
        A dictionary mapping genotypes to lists of Allele objects.
    variant_pool_numeric (dict[str, int]):
        A dictionary mapping variants to their counts.

    Returns
    -------
    list[str]:
        A list of called genotypes based on processed data.

    Example:
    ----------

    'LU*01.19/LU*02' not possible because modifying SNP (45322744_A_G) is hom
    Allele(genotype='LU*02',
            phenotype='LU:2',
            genotype_alt='LU*B',
            phenotype_alt='Lu(a-b+)',
            sample='',
            defining_variants=frozenset({'19:45315445_ref'}),
            weight_geno=LOW_WEIGHT,
            weight_pheno=2,
            reference=True,
            Sub_type='LU*02'),
    Allele(genotype='LU*02.19',
            phenotype='LU:-18,19',
            genotype_alt='.',
            phenotype_alt='Au(a-b+)',
            sample='',
            defining_variants=frozenset({'19:45315445_ref',
                                        '19:45322744_A_G'}), hom
            weight_geno=LOW_WEIGHT,
            weight_pheno=1,
            reference=False,
            Sub_type='LU*02')]},
    sample='128',
    variant_pool={'19:45315445_G_A': 'Heterozygous',
                '19:45315445_ref': 'Heterozygous',
                '19:45322744_A_G': 'Homozygous'},
    genotypes=['LU*01.19/LU*02', 'LU*01.19/LU*02.19'],
    """

    to_remove = []
    if bg.type in antitheticals:
        modifying_SNP = None
        flattened_alleles = flatten_alleles(bg.alleles[AlleleState.NORMAL])
        d = defaultdict(set)
        for allele in flattened_alleles:
            if allele.number_of_defining_variants > 1:
                for variant in allele.defining_variants:
                    if variant.split(":")[1] not in antitheticals[bg.type]:
                        d[allele.sub_type].add(variant)

        if len(d) > 1:
            assert len(d) == 2
            putative_mod_SNPs = set.union(*d.values())
            if len(putative_mod_SNPs) == 1:
                modifying_SNP = putative_mod_SNPs.pop()
        if modifying_SNP is not None and bg.variant_pool[modifying_SNP] == Zygosity.HOM:
            for pair in bg.alleles[AlleleState.NORMAL]:
                for allele in pair:
                    if allele.number_of_defining_variants == 1:
                        variant = list(allele.defining_variants)[0]
                        if variant.split(":")[1] in antitheticals[bg.type]:
                            to_remove.append(pair)

    bg.remove_pairs(to_remove, "filter_pairs_on_antithetical_modyfying_SNP")
    return bg


@apply_to_dict_values
def cant_pair_with_ref_cuz_SNPs_must_be_on_other_side(bg: BloodGroup) -> BloodGroup:
    """Filter out allele pairs where reference alleles cannot pair with non-
    reference alleles due to SNP strand requirements.

    This function examines allele pairs from the BloodGroup's NORMAL allele
    state. For each pair that contains a reference allele but not all alleles are
    reference, it splits the pair to analyze the non-reference allele. It then
    determines which SNPs must be on the other strand based on the variant pool's
    zygosity annotations. If a non-reference allele's defining variants are all
    present in the set of SNPs that must be on the other strand (or additional
    criteria are met), the pair is removed.

    Args:
        bg (BloodGroup): A BloodGroup object containing allele states, variant pool,
            and other related data.

    Returns:
        BloodGroup: The updated BloodGroup object after filtering out invalid allele
            pairs.

    Example:
        Given allele pairs such as:
          - JK*01W.03 and JK*01W.04 cannot pair with a reference allele due to SNPs
            like '18:43310313_G_A' and '18:43311054_G_A' requiring specific strand
            orientation.
        Such pairs will be removed from the BloodGroup.

    pair: [Allele(genotype='JK*01W.04',
                  phenotype='.',
                  genotype_alt='.',
                  phenotype_alt='Jk(a+ᵂ)',
                  defining_variants=frozenset({'18:43311054_G_A'}),
                  weight_geno=1000,
                  weight_pheno=3,
                  reference=False,
                  sub_type='JK*01',
                  phases=None,
                  number_of_defining_variants=1),
           Allele(genotype='JK*01', #means HOM for 43319519_ref
                  defining_variants=frozenset({'18:43319519_ref'}),
                  reference=True,

    allele3: Allele(genotype='JK*01W.03',
                    defining_variants=frozenset({'18:43310313_G_A'}),

    bg.variant_pool: {'18:43310313_G_A': 'Heterozygous',
                      '18:43311054_G_A': 'Heterozygous',
                      '18:43311131_G_A': 'Heterozygous',
                      '18:43316538_A_G': 'Heterozygous'}
    SNPS_that_need_to_be_on_other_strand: ['18:43310313_G_A']
    flattened_alleles: {Allele(genotype='JK*01W.03',
                               defining_variants=frozenset({'18:43310313_G_A'}),

                        Allele(genotype='JK*01W.04',
                               defining_variants=frozenset({'18:43311054_G_A'}),

                        Allele(genotype='JK*01N.20',
                               defining_variants=frozenset({'18:43310313_G_A',
                                                            '18:43311054_G_A',
                                                            '18:43311131_G_A',
                                                            '18:43316538_A_G'}),

                        Allele(genotype='JK*01W.11',
                               defining_variants=frozenset({'18:43311054_G_A',
                                                      '18:43310313_G_A'}),
    """
    flattened_alleles = {
        allele
        for allele in flatten_alleles(bg.alleles[AlleleState.NORMAL])
        if not allele.reference
    }
    for allele_type in [AlleleState.NORMAL]:
        if bg.alleles[allele_type] is None:
            continue
        to_remove = []

        # step1 - for non ref allele find which SNPs have to be on other side
        # (because if they were on same side the allele in question woud be something
        # else)
        for pair in bg.alleles[allele_type]:
            if pair.contains_reference and not pair.all_reference:
                ref, allele = split_pair_by_ref(pair)
                for allele2 in flattened_alleles:
                    if allele in allele2:
                        SNPS_that_need_to_be_on_other_strand = [
                            SNP
                            for SNP, zygosity in bg.variant_pool.items()
                            if zygosity == Zygosity.HOM
                        ]

                        het_SNPS = []
                        for SNP in allele2.defining_variants.difference(
                            allele.defining_variants
                        ):
                            if bg.variant_pool[SNP] == Zygosity.HET:
                                het_SNPS.append(SNP)
                        if len(het_SNPS) == 1:
                            SNPS_that_need_to_be_on_other_strand += het_SNPS
                            for SNP in allele.defining_variants:
                                for SNP2 in bg.variant_pool:
                                    if SNP == SNP2:
                                        continue
                                    if SNP.split("_")[0] in SNP2:
                                        SNPS_that_need_to_be_on_other_strand.append(
                                            SNP2
                                        )
                            # step2 - if the SNPs that must be on the ref side + and
                            # HOMs define anything then ref not possible
                            for allele3 in flattened_alleles:
                                if all(
                                    SNP in SNPS_that_need_to_be_on_other_strand
                                    for SNP in allele3.defining_variants
                                ):
                                    to_remove.append(pair)
        bg.remove_pairs(
            to_remove, "cant_pair_with_ref_cuz_SNPs_must_be_on_other_side", allele_type
        )

    return bg


@apply_to_dict_values
def ABO_cant_pair_with_ref_cuz_261delG_HET(bg: BloodGroup) -> BloodGroup:
    """Ensure ABO allele pairs obey SNP strand constraints for reference alleles.

    When a reference allele is present (defined by the Lane variant 261del,
    i.e. 9:136132908_T_TC in GRCh37), the paired non-reference allele must be an
    O allele. Some pairings (e.g. ABO*A1.01 with ABO*AEL or ABO*AW variants) are
    not allowed if 136132908_T_TC is present in the defining variants. However,
    pairings with ABO*O alleles are allowed.

    Examples:
        ABO*A1.01/ABO*AEL.02 - Not allowed as 136132908_T_TC is in defining vars.
        ABO*A1.01/ABO*AEL.07 - Not allowed as 136132908_T_TC is in defining vars.
        ABO*A1.01/ABO*AW.25 - Not allowed as 136132908_T_TC is in defining vars.
        ABO*A1.01/ABO*AW.31.01 - Not allowed as 136132908_T_TC is in defining vars.
        ABO*A1.01/ABO*O.01.05 - Allowed.
        ABO*A1.01/ABO*O.01.22 - Allowed.
        ABO*A1.01/ABO*O.01.45 - Allowed.
        ABO*A1.01/ABO*O.01.71 - Allowed.

    Args:
        bg (BloodGroup): A BloodGroup object containing allele pairs, a variant
            pool, and other related genetic data.

    Returns:
        BloodGroup: The updated BloodGroup after filtering out allele pairs where a
            reference allele is improperly paired.
    """

    for allele_type in [AlleleState.NORMAL]:
        if bg.alleles[allele_type] is None:
            continue
        to_remove = []
        for pair in bg.alleles[allele_type]:
            if pair.contains_reference and not pair.all_reference:
                ref, allele = split_pair_by_ref(pair)
                tmp_pool2 = bg.variant_pool_numeric
                for variant_on_other_strand in ref.defining_variants:
                    if variant_on_other_strand in tmp_pool2:
                        tmp_pool2[variant_on_other_strand] -= 1
                check_vars_other_strand = partial(
                    check_available_variants, 0, tmp_pool2, operator.gt
                )
                if all(check_vars_other_strand(allele)):
                    # ie, can they exist given other strand
                    continue
                to_remove.append(pair)

        bg.remove_pairs(
            to_remove, "ABO_cant_pair_with_ref_cuz_261delG_HET", allele_type
        )

    return bg


@apply_to_dict_values
def cant_pair_with_ref_cuz_trumped(bg: BloodGroup) -> BloodGroup:
    """Filter out allele pairs where a reference allele is trumped by a superior allele.

    This function checks allele pairs in the NORMAL state that contain a reference
    allele. If a non-reference allele in the pair is trumped (i.e. has a higher
    weight_geno compared to another allele with one HET variant in the same subtype),
    the pair is removed.


    Args:
        bg (BloodGroup): A BloodGroup object containing allele pairs, variant pool,
            and related genetic data.

    Returns:
        BloodGroup: The updated BloodGroup after filtering out trumped allele pairs.


    Examples:
    ----------
    summary:
        BG Name: FUT3
        Pairs:
            - FUT3*01N.01.02/FUT3*01N.01.12
            - FUT3*01.16/FUT3*01N.01.12
            - FUT3*01/FUT3*01N.01.12
            - FUT3*01/FUT3*01N.01.02
        Filtered Out:
            - FUT3*01/FUT3*01.16 is removed because FUT3*01.16 is trumped by a better allele,
              e.g. FUT3*01N.01.02 or FUT3*01N.01.12.

    BG Name: FUT3
    Pairs:
    Pair(Genotype: FUT3*01N.01.02/FUT3*01N.01.12 Phenotype: ./.)
    Pair(Genotype: FUT3*01.16/FUT3*01N.01.12 Phenotype: ./.)
    Pair(Genotype: FUT3*01/FUT3*01N.01.12 Phenotype: ./.)
    Pair(Genotype: FUT3*01/FUT3*01N.01.02 Phenotype: ./.)

    Filtered Out: defaultdict(<class 'list'>, {'cant_pair_with_ref_cuz_trumped':
    [Pair(Genotype: FUT3*01/FUT3*01.16 Phenotype: ./.)]})

    FUT3*01.16 not with ref due to FUT3*01N.01.02 and FUT3*01N.01.12

    Current truth: FUT3 126

    FUT3*01.16/FUT3*01N.01.12
    FUT3*01/FUT3*01N.01.02
    FUT3*01/FUT3*01N.01.12
    FUT3*01N.01.02/FUT3*01N.01.12
    [Allele(genotype='FUT3*01.15',
            defining_variants=frozenset({'19:5844184_C_T',
                        '19:5844367_C_T'}), HOM
            weight_geno=1000,
            reference=False,
    Allele(genotype='FUT3*01.16',
            defining_variants=frozenset({'19:5844043_C_T',
                                        '19:5844184_C_T',
                                        '19:5844367_C_T'}), HOM
            weight_geno=1000,
            reference=False,
    Allele(genotype='FUT3*01N.01.02',
            defining_variants=frozenset({'19:5844184_C_T',
                                        '19:5844367_C_T', HOM
                                        '19:5844838_C_T'}), HOM
            weight_geno=1,
            reference=False,
    Allele(genotype='FUT3*01N.01.12',
            defining_variants=frozenset({'19:5843883_C_G',
                                        '19:5844367_C_T', HOM
                                        '19:5844838_C_T'}), HOM
            weight_geno=1,
            reference=False,

    res[goi].variant_pool: {'19:5843883_C_G': 'Heterozygous',
                            '19:5844043_C_T': 'Heterozygous',
                            '19:5844184_C_T': 'Heterozygous',
                            '19:5844367_C_T': 'Homozygous',
                            '19:5844838_C_T': 'Homozygous'}
    """

    for allele_type in [AlleleState.NORMAL]:
        # TODO remove evolutionary baggage from co_exists everywhere - after making
        # unit
        # tests
        if bg.alleles[allele_type] is None:
            continue
        to_remove = []
        flattened_alleles = flatten_alleles(bg.alleles[allele_type])

        if not any(allele.reference for allele in flattened_alleles):
            continue
        alleles_without_ref = [
            allele for allele in flattened_alleles if not allele.reference
        ]
        if allele_type == AlleleState.NORMAL:
            alleles_with_1_HET_var = [
                allele
                for allele in alleles_without_ref
                if [
                    bg.variant_pool[variant] for variant in allele.defining_variants
                ].count(Zygosity.HET)
                == 1
            ]
        if alleles_with_1_HET_var:
            for pair in bg.alleles[allele_type]:
                if pair.contains_reference and not pair.all_reference:
                    ref, allele = split_pair_by_ref(pair)
                    for allele2 in alleles_with_1_HET_var:
                        if allele.sub_type != allele2.sub_type:
                            continue
                        if allele.weight_geno > allele2.weight_geno:
                            to_remove.append(pair)
        if to_remove:
            bg.remove_pairs(to_remove, "cant_pair_with_ref_cuz_trumped", allele_type)

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

    """

    for allele_type in [AlleleState.CO]:
        if bg.alleles[allele_type] is None:
            continue
        to_remove = []
        for variant, zygo in bg.variant_pool.items():
            if zygo == Zygosity.HET:
                for pair in bg.alleles[allele_type]:
                    allele1_vars_plus_het_var = set(pair.allele1.defining_variants) | {
                        variant
                    }
                    allele2_vars_plus_het_var = set(pair.allele2.defining_variants) | {
                        variant
                    }
                    flattened = {
                        allele
                        for pair2 in bg.alleles[allele_type]
                        for allele in pair2
                        if pair2 != pair
                    }
                    hits = 0
                    for a in flattened:
                        if a not in pair:
                            flat_vars = set(a.defining_variants)
                            if (
                                flat_vars == allele1_vars_plus_het_var
                                or flat_vars == allele2_vars_plus_het_var
                            ):
                                hits += 1
                    if hits > 1:
                        to_remove.append(pair)
        bg.remove_pairs(to_remove, "ensure_co_existing_HET_SNP_used", allele_type)

    return bg


@apply_to_dict_values
def filter_HET_pairs_by_weight(bg: BloodGroup) -> BloodGroup:
    """This filter forced us to make decisions where there was not always
    a clearly correct answer. I have left the example long to show the evolution of
    thought and why we landed where we did.

    Args:
        bg (BloodGroup): The BloodGroup object containing allele pairs and variant pool
        information.

    Returns:
        BloodGroup: The updated BloodGroup object with inconsistent allele pairs
        removed.

    Example
    ----------
     FUT1/2
    All bar ref here is HET
    So while it's possible that any allele can exist in a pair,
    no pair can exist without most weighted allele, if it just has 1 defining SNP
    (or should it be just has 1 HET defining SNP, cuz if the rest HOM, same, right?)
    # update - only if alleles are same subtype! and if the 1 SNP of the most weighted
    # allele is 'in' one of the pair (meaning that most weighted allele doesn't exist?)
    . Ref has to be included (as all SNPS could be together)
    but FUT2*01/FUT2*01.03.01 isn't possible
    Further, FUT2*01/FUT2*01N.16 isn't possible as FUT2*01.03.01 trumps FUT2*01N.16
    [Allele(genotype='FUT2*01',
        genotype_alt='.',
        defining_variants=frozenset({'19:49206250_ref'}), hom
        weight_geno=LOW_WEIGHT,
        weight_pheno=2,
        reference=True,
        Sub_type='FUT2*01'),
    Allele(genotype='FUT2*01.03.01',
        genotype_alt='.',
        defining_variants=frozenset({'19:49206286_A_G'}),
        weight_geno=LOW_WEIGHT,
        weight_pheno=1,
        reference=False,
        Sub_type='FUT2*01'),
    Allele(genotype='FUT2*01N.02',
        genotype_alt='.',
        defining_variants=frozenset({'19:49206674_G_A'}),
        weight_geno=1,
        weight_pheno=5,
        reference=False,
        Sub_type='FUT2*01'),
    Allele(genotype='FUT2*01N.16',
        genotype_alt='.',
        defining_variants=frozenset({'19:49206985_G_A'}),
        weight_geno=8,
        weight_pheno=5,
        reference=False,
        Sub_type='FUT2*01')]

    res["FUT2"].variant_pool: {'19:49206250_ref': 'Homozygous',
                               '19:49206286_A_G': 'Heterozygous',
                               '19:49206674_G_A': 'Heterozygous',
                               '19:49206985_G_A': 'Heterozygous'}


    issue

    [Allele(genotype='JK*01N.19',
                        genotype_alt='.',
                        defining_variants=frozenset({'18:43319274_G_A'}),
                        weight_geno=7,
                        reference=False,
                        sub_type='JK*01',

    Allele(genotype='JK*02N.17',
                        genotype_alt='.',
                        defining_variants=frozenset({'18:43319274_G_A',
                                                    '18:43319519_G_A'}),
                        weight_geno=13,
                        reference=False,
                        sub_type='JK*02',

    Allele(genotype='JK*01W.06',
                        genotype_alt='.',
                        defining_variants=frozenset({'18:43310415_G_A',
                                                    '18:43316538_A_G'}), hom
                        weight_geno=1000,
                        reference=False,
                        sub_type='JK*01',

    JK*02	43319519_G_A

    pool: {'18:43310415_G_A': 'Heterozygous',
                '18:43316538_A_G': 'Homozygous',
                '18:43319274_G_A': 'Heterozygous',
                '18:43319519_G_A': 'Heterozygous'}

    current code removed JK*01W.06/JK*02N.17. At first glance I thought this was the
    correct behaviour
    because JK*01N.19 trumps JK*01W.06 and only has 1 SNP. but I noticed
    that the 1 defining SNP for JK*01N.19 (43319274_G_A) is part of the definition of
    JK*02N.17, which is further complicated by the fact its a different subtype.
    Based on my understanding I think that JK*01W.06/JK*02N.17 is right, because
    the trumpiness of JK*01N.19 is lost to JK*02N.17.

    Is there a clear rule to cover this problem? One the one hand, nulls are always
    the most weighted allele (they trump other options) but at the same time alleles
    whose defining variants are all in the list of defining variants of another alleles
    are no longer considered. Which of these rules is applied first? In this example
    the 1 defining SNP for JK*01N.19 (43319274_G_A) is part of the definition of
    JK*02N.17, so JK*01N.19 'doesn't exist' anymore. That makes sense to me, but if
    JK*01N.19 caused truncation or frameshift then I would say that this should
    be the most important consideration... unless the addition of the other variant to
    make JK*02N.17 rescues from being null? As I think this through I think I'm flipping
    and would think that any SNP that makes null would have to be to first consideration

    assumig we want to handle 'in' and subtype:
        each pair will have just 1 of each subtype (other filters will remove)
        so if pair has to include most weighted allele of subtype (sigh)
            if that allele has


    TODO
    need to also handle the case where the most weighted has multiple SNPs, but the
    2nd or 3rd most weighted has 1 SNP... ffs
    """

    flattened_alleles = flatten_alleles(bg.alleles[AlleleState.NORMAL])
    if not flattened_alleles:
        return bg
    if all(allele.weight_geno == LOW_WEIGHT for allele in flattened_alleles):
        return bg
    sub_types = {allele.sub_type for allele in flattened_alleles}

    keepers = []

    for sub_type in sub_types:
        max_weight = min(
            [
                allele.weight_geno
                for allele in flattened_alleles
                if allele.sub_type == sub_type
            ]
        )  # lowest weight, base 1

        weights_with_1_SNP = [
            allele.weight_geno
            for allele in flattened_alleles
            if one_HET_variant(allele, bg.variant_pool) and allele.sub_type == sub_type
        ]
        max_weight_with_1_SNP = (
            min(weights_with_1_SNP) if weights_with_1_SNP else LOW_WEIGHT
        )

        if max_weight_with_1_SNP == max_weight:
            alleles_with_max_weight_and_1_SNP = [
                allele
                for allele in flattened_alleles
                if one_HET_variant(allele, bg.variant_pool)
                and allele.sub_type == sub_type
                and allele.weight_geno == max_weight
            ]
            for pair in bg.alleles[AlleleState.NORMAL]:
                if all(
                    allele in pair.allele1 or allele in pair.allele2
                    for allele in alleles_with_max_weight_and_1_SNP
                ) or any(
                    allele.weight_geno == max_weight
                    for allele in pair
                    if allele.sub_type == sub_type
                ):
                    keepers.append(pair)
        else:  # if the allele with max weight does't have just 1 het SNP, need to keep all
            # TODO make good unit tests then refactor this mess
            for pair in bg.alleles[AlleleState.NORMAL]:
                keepers.append(pair)
    to_remove = [pair for pair in bg.alleles[AlleleState.NORMAL] if pair not in keepers]

    bg.remove_pairs(to_remove, "filter_HET_pairs_by_weight")

    return bg


### does this ever come up?TODO e2e failed when removed.
# issue is that these filters aren't
# being recorded properly!!! (still some might be unused)
@apply_to_dict_values
def filter_pairs_by_context(bg: BloodGroup) -> BloodGroup:
    """Filter allele pairs by context.

    This function removes allele pairs from the NORMAL state of the BloodGroup if
    the context indicates that the pair cannot exist. It checks whether the
    combination of remaining variant counts and an allele's defining variants
    match those in other alleles. This ensures that pairs like
    'A4GALT*01/A4GALT*02' are filtered out when a more comprehensive allele is
    implied by the available variants. This logic overlaps with
    cant_pair_with_ref_cuz_SNPs_must_be_on_other_side, but is still needed

    Example:
        Given allele definitions:
            - A4GALT*01.02: defining_variants = {'22:43089849_T_C'}
            - A4GALT*02:   defining_variants = {'22:43113793_C_A'}
            - A4GALT*02.02: defining_variants = {'22:43113793_C_A',
                                                 '22:43089849_T_C'}
            - A4GALT*01:   defining_variants = {'22:43113793_ref'}
        And a variant pool:
            {'22:43089849_T_C': 'Heterozygous',
             '22:43113793_C_A': 'Heterozygous'}
        The valid pairs are:
            'A4GALT*01.02/A4GALT*02' and
            'A4GALT*01/A4GALT*02.02',
        while 'A4GALT*01/A4GALT*02' is not possible.

    Args:
        bg (BloodGroup): A BloodGroup object containing allele pairs and a numeric
            variant pool.

    Returns:
        BloodGroup: The updated BloodGroup with contextually invalid pairs removed.

    """
    # need to catch; if remaing var and allele var combine to define another allele
    # and if remaining var defines or combines with other allele to define another
    # allele then pair is not possible, but check that the pair that rules it out
    # does already exist

    for allele_type in [AlleleState.NORMAL]:
        if bg.alleles[allele_type] is None:
            continue
        to_remove = []
        for pair in bg.alleles[allele_type]:
            def_vars = {
                a.defining_variants
                for pair2 in bg.alleles[allele_type]
                for a in pair2
                if a not in pair
            }
            for allele in pair:
                variant_pool_copy = bg.variant_pool_numeric.copy()
                if allele.reference:
                    continue

                for variant in allele.defining_variants:
                    variant_pool_copy[variant] -= 1
                left_over_vars = [k for k, v in variant_pool_copy.items() if v > 0]
                if any([len(left_over_vars) == 0, len(def_vars) < 2]):
                    continue
                remaining = [
                    tuple(sorted(set(left_over_vars))),
                    tuple(sorted(set(left_over_vars + list(allele.defining_variants)))),
                ]
                if all(variants in def_vars for variants in remaining):
                    to_remove.append(pair)
                    break
        bg.remove_pairs(to_remove, "filter_pairs_by_context", allele_type)

    return bg


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
    """
    if phased:
        to_remove = []
        for allele in bg.alleles[AlleleState.FILT]:
            phases_with_hom = _get_allele_phase_info(allele, bg.variant_pool_phase)
            phases = set(
                [phase_set for phase_set in phases_with_hom if phase_set != "."]
            )
            phase_sets_and_hom = _get_allele_phase_info(
                allele, bg.variant_pool_phase_set
            )
            if all(phase_set == "." for phase_set in phase_sets_and_hom):
                continue  # all hom
            phase_sets = set(
                [phase_set for phase_set in phase_sets_and_hom if phase_set != "."]
            )
            # ic(allele, phase_sets, phases, bg.variant_pool_phase, bg.variant_pool_phase_set, bg.variant_pool)
            if len(phase_sets) > 2:
                # logger.info(
                #     f"\n{bg.sample} : {allele} multiple phase sets: {phase_sets}"
                # )
                continue
            elif len(phase_sets) == 2:
                # logger.info(
                #     f"\n{bg.sample} : {allele} different or unknown phase set: {phase_sets}"
                # )
                continue
            elif any("/" in phase for phase in phases):
                # logger.info(
                #     f"\n{bg.sample} : {allele} unknown phase: {phases}"
                # )
                continue
            elif len(phase_sets) == 1 and len(phases) > 1:
                # logger.warning(
                #     f"\n{bg.sample} : {allele} is not phased: {phase_sets} & {phases} \n {bg.variant_pool} \n\n{bg.variant_pool_phase} \n\n{bg.variant_pool_phase_set}"
                # )
                to_remove.append(allele)
            elif len(phase_sets) == 1 and len(phases) == 1:
                continue
            else:
                raise ValueError(
                    f"beyond logic 121212 {allele} phased: {phase_sets} & {phases}"
                )
        if to_remove:
            # ic(44455566677744, bg.type, bg.alleles[AlleleState.FILT], to_remove)
            bg.remove_alleles(to_remove, "remove_unphased", AlleleState.FILT)
    return bg


# TODO tidy up all phase funcs once tested
def _get_allele_phase_info(allele, phase_dict):
    """ """
    return [phase_dict[variant] for variant in allele.defining_variants]


# TODO tidy up all phase funcs once tested
def _get_allele_phase_info2(allele, phase_dict):
    """ """
    return [
        phase_dict[variant]
        for variant in allele.defining_variants
        if "9:133257521" not in variant and "9:136132908" not in variant
    ]


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
    """

    if phased:
        to_remove = []
        for pair in bg.alleles[AlleleState.NORMAL]:
            if pair.contains_reference:
                continue
            p1_phases = set(
                _get_allele_phase_info2(pair.allele1, bg.variant_pool_phase)
            )
            p1_zygo = set(_get_allele_phase_info(pair.allele1, bg.variant_pool))
            p1_phase_sets = set(
                _get_allele_phase_info(pair.allele1, bg.variant_pool_phase_set)
            )
            p2_phases = set(
                _get_allele_phase_info2(pair.allele2, bg.variant_pool_phase)
            )
            p2_zygo = set(_get_allele_phase_info(pair.allele2, bg.variant_pool))
            p2_phase_sets = set(
                _get_allele_phase_info(pair.allele2, bg.variant_pool_phase_set)
            )

            phase_set = p1_phase_sets.union(p2_phase_sets)
            # ic(666777777, p1_phases, p2_phases)
            if len(phase_set) != 1:
                continue  # can't use phasing info
            # if bg.type == 'JK':
            #     print(11111, pair, p1_phases, p2_phases)
            if p1_zygo == {Zygosity.HOM} and p2_zygo == {Zygosity.HOM}:  # all hom
                continue
            # if bg.type == 'JK':
            #     print(22222, pair,p1_phases, p2_phases)
            elif p1_phases == p2_phases:
                to_remove.append(pair)
            # elif remove_homs(p1) != remove_homs(p2):
            #     continue
            # elif remove_homs(p1) == remove_homs(p2):
            #     to_remove.append(pair)
        if len(bg.alleles[AlleleState.NORMAL]) == len(to_remove):
            for pair in to_remove:
                bg.alleles[AlleleState.NORMAL].append(
                    Pair(reference_alleles[bg.type], pair.allele1)
                )
                bg.alleles[AlleleState.NORMAL].append(
                    Pair(reference_alleles[bg.type], pair.allele2)
                )
        # ic(9998887776655, bg.type, bg.alleles[AlleleState.NORMAL], to_remove)
        if to_remove:
            bg.remove_pairs(to_remove, "filter_pairs_by_phase")

    return bg


@apply_to_dict_values
def filter_impossible_alleles_phased(bg: BloodGroup, phased: bool) -> BloodGroup:
    """
    Filters alleles in a BloodGroup object based on phasing consistency and subsumption.

    When `phased` is True, this function attempts to simplify the set of possible
    alleles (`bg.alleles[AlleleState.POS]`) by looking for alleles that are impossibe
    due to being 'in' (defining variants are subset of) another allele:
        1. A1 is A2, all vars HOM (should be removed above - not phased dependant)
        2. A1 is A2, all vars HET and in same phase
        2. A1 is A2, vars are mix of HET and HOM, all HET are in same phase


    The function modifies `bg.alleles[AlleleState.POS]` in-place by removing the
    alleles deemed "impossible" under these phased conditions. Details of removed
    allele pairs due to subsumption are stored in
    `bg.filtered_out["allele_subsumed_by_other_phased"]`.

    Args:
        bg (BloodGroup): A BloodGroup object containing allele states, variant pool
            (with observed zygosities), and phasing information within Allele objects.
        phased (bool): A flag indicating whether phasing rules should be applied.
            If False, the function returns the BloodGroup object unmodified.

    Returns:
        BloodGroup: The modified BloodGroup object. The `alleles[AlleleState.POS]`
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
    to_test = "ABO"
    if phased:
        if len(bg.alleles[AlleleState.NORMAL]) == 1:
            return bg
        if len(bg.alleles[AlleleState.NORMAL]) == 0:
            # raise?TODO
            # assert bg.type == 'KN'
            return bg
        alleles = list(flatten_alleles(bg.alleles[AlleleState.NORMAL]))
        # ic(alleles)
        alleles_not_hom = [
            allele for allele in alleles if not _check_hom(bg.variant_pool, allele)
        ]
        # if bg.type == to_test:
        #     ic(22222333333, [a.genotype for a in alleles_not_hom])
        alleles_with_variants_in_same_phase_set = [
            allele
            for allele in alleles_not_hom
            if _check_phase_set(bg.variant_pool_phase_set, allele)
        ]
        # if bg.type == to_test:
        #     ic(
        #         3333333444444,
        #         [a.genotype for a in alleles_with_variants_in_same_phase_set],
        #     )
        alleles_with_variants_in_same_phase = [
            allele
            for allele in alleles_with_variants_in_same_phase_set
            if _check_phase(bg.variant_pool_phase, allele)
        ]
        # if bg.type == to_test:
        #     ic(444444455555, [a.genotype for a in alleles_with_variants_in_same_phase])
        l1, l2 = [], []
        # sort to rescue 9:133257521"/ "9:136132908"
        for allele in sorted(
            alleles_with_variants_in_same_phase,
            key=lambda allele: len(allele.defining_variants),
            reverse=True,
        ):
            # ic(222222444444, allele.genotype, len(allele.defining_variants))
            phases = set([])
            for variant in allele.defining_variants:
                if "9:133257521" not in variant and "9:136132908" not in variant:
                    phase = bg.variant_pool_phase[variant]
                    if phase != "1/1":
                        phases.add(bg.variant_pool_phase[variant])
            if not phases:
                if bg.type == 'ABO' and any(
                    [
                        "9:133257521" in variant or "9:136132908" in variant
                        for variant in allele.defining_variants
                    ]
                ):
                    if l1 and l2:
                        #ic(11111111199999, allele.sub_type)
                        if variant.endswith('_ref') and 'ABO*O.01.' in allele.genotype:
                            if 'ABO*O.01.' in l1[0].genotype:
                                l1.append(allele)
                            else:
                                l2.append(allele)
                        elif variant.endswith('_T_TC') and allele.sub_type != 'ABO*O':
                            if l1[0].sub_type in ['ABO*A', 'ABO*B']:
                                l1.append(allele)
                            else:
                                l2.append(allele)

                #ic("no phase", allele.genotype)
                continue
            assert len(phases) == 1
            phase = phases.pop()
            if phase == "1|0":
                l1.append(allele)
            elif phase == "0|1":
                l2.append(allele)
            else:
                assert phase in "unknown" or "/" in phase
                # ic(11111111, allele, phase)
        # if bg.type == to_test:
        #     ic(len([a.genotype for a in l1]), len([a.genotype for a in l2]))
        alleles_to_remove = []
        for allele in l1:
            for allele2 in l1:
                if allele in allele2:
                    alleles_to_remove.append(allele)
        for allele in l2:
            for allele2 in l2:
                if allele in allele2:
                    alleles_to_remove.append(allele)
        # if bg.type == to_test:
        #     ic(len([a.genotype for a in alleles_to_remove]))
        to_remove = []
        for pair in bg.alleles[AlleleState.NORMAL]:
            if pair.allele1 in alleles_to_remove or pair.allele2 in alleles_to_remove:
                to_remove.append(pair)
        try:
            assert len(to_remove) != len(bg.alleles[AlleleState.NORMAL])
        except AssertionError:
            # ic(
            #     11111222233333897,
            #     to_remove,
            #     bg.alleles[AlleleState.NORMAL],
            #     len(to_remove),
            #     len(bg.alleles[AlleleState.NORMAL]),
            # )
            raise
        if to_remove:
            bg.remove_pairs(to_remove, "filter_impossible_alleles_phased")
        try:
            assert bg.alleles[AlleleState.NORMAL]
        except AssertionError:
            # ic(bg.type)
            raise
    return bg


def _check_phase(variant_pool: dict[str, str], current_allele: Allele) -> bool:
    """
    True if all same phase set and or HOM
    """

    phase_sets = [
        phase
        for variant, phase in variant_pool.items()
        if variant in current_allele.defining_variants
        and phase != "1/1"
        and "9:133257521" not in variant
        and "9:136132908" not in variant
    ]  # 9:133257521 == 9:136132908 == c.261delG locus, is like an atithetical switch for most ABO-O
    # but it doesn't always get phased properly, so it's given special consideration (excluded)

    # if current_allele.genotype in ['ABO*O.01.05','ABO*O.01.24','ABO*O.01.41']:
    #     ic(
    #         456789,
    #         current_allele,
    #         variant_pool,
    #         phase_sets,
    #         set(phase_sets),
    #         len(set(phase_sets)) == 1,
    #     )
    if not phase_sets and "ABO" in current_allele.genotype:
        #ic("wtf", phase_sets, current_allele.genotype)
        return True
    return len(set(phase_sets)) == 1


def _check_phase_set(variant_pool: dict[str, str], current_allele: Allele) -> bool:
    """
    True if all same phase set and or HOM
    """
    phase_sets = [
        phase
        for variant, phase in variant_pool.items()
        if variant in current_allele.defining_variants and phase != "."
    ]

    return len(set(phase_sets)) == 1


def _check_hom(variant_pool: dict[str, str], current_allele: Allele) -> bool:
    """
    Checks if all elements in a list are hom.

    Args:
        current_allele: Allele

    Returns:
        True if all elements match the target_str, or if all elements
        except exactly one match the target_str. False otherwise.
    """
    variants = [
        zygo
        for variant, zygo in variant_pool.items()
        if variant in current_allele.defining_variants
    ]

    return all(item == Zygosity.HOM for item in variants)


@apply_to_dict_values
def filter_impossible_coexisting_alleles_phased(
    bg: BloodGroup, phased: bool
) -> BloodGroup:
    """
    Filters alleles in a BloodGroup object based on phasing consistency and subsumption.

    When `phased` is True, this function attempts to simplify the set of possible
    alleles (`bg.alleles[AlleleState.POS]`) by looking for alleles that are impossibe
    due to being 'in' (defining variants are subset of) another allele:
        1. A1 is A2, all vars HOM (should be removed above - not phased dependant)
        2. A1 is A2, all vars HET and in same phase
        2. A1 is A2, vars are mix of HET and HOM, all HET are in same phase


    The function modifies `bg.alleles[AlleleState.POS]` in-place by removing the
    alleles deemed "impossible" under these phased conditions. Details of removed
    allele pairs due to subsumption are stored in
    `bg.filtered_out["allele_subsumed_by_other_phased"]`.

    Args:
        bg (BloodGroup): A BloodGroup object containing allele states, variant pool
            (with observed zygosities), and phasing information within Allele objects.
        phased (bool): A flag indicating whether phasing rules should be applied.
            If False, the function returns the BloodGroup object unmodified.

    Returns:
        BloodGroup: The modified BloodGroup object. The `alleles[AlleleState.POS]`
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
    if phased:
        co_alleles = bg.alleles.get(AlleleState.CO)
        if co_alleles is None:
            return bg

        if len(co_alleles) == 1:
            return bg
        alleles = list(flatten_alleles(co_alleles))
        # ic(alleles)
        alleles_not_hom = [
            allele for allele in alleles if not _check_hom(bg.variant_pool, allele)
        ]
        # ic(bg.variant_pool, alleles_not_hom)
        alleles_with_variants_in_same_phase_set = [
            allele
            for allele in alleles_not_hom
            if _check_phase_set(bg.variant_pool_phase_set, allele)
        ]
        # ic(bg.variant_pool_phase_set, alleles_with_variants_in_same_phase_set)
        alleles_with_variants_in_same_phase = [
            allele
            for allele in alleles_with_variants_in_same_phase_set
            if _check_phase(bg.variant_pool_phase, allele)
        ]
        # ic(bg.variant_pool_phase, alleles_with_variants_in_same_phase)
        l1, l2 = [], []
        for allele in alleles_with_variants_in_same_phase:
            phase = bg.variant_pool_phase[list(allele.defining_variants)[0]]
            if phase == "1|0":
                l1.append(allele)
            elif phase == "0|1":
                l2.append(allele)
            else:
                try:
                    assert phase == '1/1'
                except AssertionError:
                    ic(11111111, allele, phase)
                    raise
        # ic(bg.variant_pool_phase, l1,l2)
        alleles_to_remove = []
        for allele in l1:
            for allele2 in l1:
                if allele in allele2:
                    alleles_to_remove.append(allele)
        for allele in l2:
            for allele2 in l2:
                if allele in allele2:
                    alleles_to_remove.append(allele)
        # ic(alleles_to_remove)
        to_remove = []
        for pair in co_alleles:
            if pair.allele1 in alleles_to_remove or pair.allele2 in alleles_to_remove:
                to_remove.append(pair)
        assert len(to_remove) != len(co_alleles)
        if to_remove:
            bg.remove_pairs(
                to_remove, "filter_impossible_co_alleles_phased", AlleleState.CO
            )
    return bg
    # def check_hom(current_allele: Allele) -> bool:
    #     """
    #     Checks if all elements in a list are hom.

    #     Args:
    #         current_allele: Allele

    #     Returns:
    #         True if all elements match the target_str, or if all elements
    #         except exactly one match the target_str. False otherwise.
    #     """
    #     variants = [
    #         zygo
    #         for variant, zygo in bg.variant_pool.items()
    #         if variant in current_allele.defining_variants
    #     ]

    #     return all(item == Zygosity.HOM for item in variants)

    # if bg.type != "KN":
    #     return bg
    # alleles_not_flat = bg.alleles.get(AlleleState.CO)
    # if alleles_not_flat is not None:
    #     alleles = list(flatten_alleles(alleles_not_flat))
    # else:
    #     return bg
    # to_test = "KN"
    # if phased:
    #     impossible_alleles = set([])
    #     impossible_pairs = set([])
    #     for allele in alleles:
    #         # if bg.type == to_test:
    #         #     ic(11111, allele)
    #         if allele in impossible_alleles:
    #             # if bg.type == to_test:
    #             #     ic(22222)
    #             continue

    #         if allele.genotype not in ["KN*01", "KN*02"]:
    #             continue
    #         if check_hom(allele):
    #             # if bg.type == to_test:
    #             #     ic(3333333)
    #             continue
    #         for allele2 in alleles:
    #             # if bg.type == to_test:
    #             #     ic(44444444, allele2)
    #             if allele == allele2:
    #                 # if bg.type == to_test:
    #                 #     ic(555555)
    #                 continue
    #             # if bg.type == 'KN' and not allele2.genotype not in ['KN*01','KN*02']:
    #             #     continue
    #             if check_hom(allele2):
    #                 # if bg.type == to_test:
    #                 #     ic(6666666)
    #                 continue
    #             if allele in allele2:
    #                 # if bg.type == to_test:
    #                 #     ic(7777777,allele.genotype, allele2.genotype)
    #                 impossible_alleles.add(allele)
    #                 impossible_pairs.add(Pair(allele1=allele2, allele2=allele))
    #                 continue
    #             if allele2 in allele:
    #                 # if bg.type == to_test:
    #                 #     ic(888888, allele.genotype, allele2.genotype)
    #                 impossible_alleles.add(allele2)
    #                 impossible_pairs.add(Pair(allele1=allele2, allele2=allele))
    #                 continue

    #     bg.filtered_out["allele1_not_possible_due_to_allele2_co_phased"] = list(
    #         impossible_pairs
    #     )
    #     # if bg.type == to_test:
    #     #     ic(bg.alleles[AlleleState.POS])
    #     to_remove = []

    #     co_alleles = bg.alleles.get(AlleleState.CO)
    #     for pair in co_alleles:
    #         if bg.type == to_test:
    #             ic(pair, impossible_alleles)
    #         if pair.allele1 in impossible_alleles or pair.allele2 in impossible_alleles:
    #             to_remove.append(pair)
    #     bg.remove_pairs(
    #         to_remove, "filter_impossible_alleles_co_phased", AlleleState.CO
    #     )

    # return bg


# @apply_to_dict_values
# def filter_pairs_by_phase_set(
#     bg: BloodGroup, phased: bool, reference_alleles
# ) -> BloodGroup:
#     """
#     Filters out allele pairs where both alleles are on the same phase strand.

#     This function is intended to remove allele pairs from a BloodGroup object when both
#     alleles in a pair are on the same phase strand, indicating they are on the same
#     chromosome and cannot be inherited together. The function operates under the
#     following logic:

#     - If `phased` is False, the function returns the BloodGroup object unchanged.
#     - For each allele pair in `bg.alleles[AlleleState.NORMAL]`:
#         - If the pair contains a reference allele, it is retained.
#         - Extract the phase sets (`p1` and `p2`) for each allele in the pair.
#         - If both alleles are homozygous (phase sets are {"."}), the pair is retained.
#         - If the phase sets are identical, the pair is removed.
#         - If the non-homozygous phase sets differ, the pair is retained.
#         - If the non-homozygous phase sets are identical, the pair is removed.
#         - If none of the above conditions are met, a ValueError is raised.

#     - If all pairs are removed and there were pairs with phase information, new pairs
#     are created by pairing each allele with the reference allele for the blood group
#     type.

#     Args:
#         bg (BloodGroup): The BloodGroup object containing allele pairs.
#         phased (bool): A flag indicating whether phase information is available.
#         reference_alleles (dict): A dictionary mapping blood group types to reference
#         alleles.

#     Returns:
#         BloodGroup: The updated BloodGroup object with inconsistent allele pairs removed.


#     Example:
#     ----------
#     Suppose you have allele pairs where both alleles are on the same phase strand.
#     This function will remove such pairs, ensuring that only valid allele combinations
#     are retained. If all pairs are removed and phase information is present, it will
#     create new pairs with the reference allele to represent possible allele
#     combinations.


#     Meant to remove pairs where both alleles are on the same strand ie
#     to_remove: [[Allele(genotype='FUT2*01N.16',
#                         defining_variants=frozenset({'19:48703728_G_A'}),
#                         weight_geno=8,
#                         weight_pheno=5,
#                         reference=False,
#                         sub_type='FUT2*01',
#                         phase=48649018),
#                  Allele(genotype='FUT2*01N.02',
#                         defining_variants=frozenset({'19:48703417_G_A'}),
#                         weight_geno=1,
#                         weight_pheno=5,
#                         reference=False,
#                         sub_type='FUT2*01',
#                         phase=48649018)]]

#     dont remove if ref in pair
#     if there is only 1 pair and they are phased then change to 2 pairs (or &) with ref
#     """

#     def remove_homs(phase_sets):
#         return [phase_set for phase_set in phase_sets if phase_set != "."]

#     if phased:
#         to_remove = []
#         for pair in bg.alleles[AlleleState.NORMAL]:
#             if pair.contains_reference:
#                 continue
#             ic(pair.allele1, pair.allele2)
#             p1 = set(pair.allele1.phases)
#             p2 = set(pair.allele2.phases)
#             if p1 == {"."} and p2 == {"."}:  # all hom
#                 continue
#             elif p1 == p2:
#                 to_remove.append(pair)
#             elif remove_homs(p1) != remove_homs(p2):
#                 continue
#             elif remove_homs(p1) == remove_homs(p2):
#                 to_remove.append(pair)
#         if len(bg.alleles[AlleleState.NORMAL]) == len(to_remove):
#             for pair in to_remove:
#                 bg.alleles[AlleleState.NORMAL].append(
#                     Pair(reference_alleles[bg.type], pair.allele1)
#                 )
#                 bg.alleles[AlleleState.NORMAL].append(
#                     Pair(reference_alleles[bg.type], pair.allele2)
#                 )
#         bg.remove_pairs(to_remove, "filter_pairs_by_phase")

#     return bg

# @apply_to_dict_values
# def filter_impossible_alleles_phased_old(bg: BloodGroup, phased: bool) -> BloodGroup:
#     """
#     Filters alleles in a BloodGroup object based on phasing consistency and subsumption.

#     When `phased` is True, this function attempts to simplify the set of possible
#     alleles (`bg.alleles[AlleleState.POS]`) by looking for alleles that are impossibe
#     due to being 'in' (defining variants are subset of) another allele:
#         1. A1 is A2, all vars HOM (should be removed above - not phased dependant)
#         2. A1 is A2, all vars HET and in same phase
#         2. A1 is A2, vars are mix of HET and HOM, all HET are in same phase


#     The function modifies `bg.alleles[AlleleState.POS]` in-place by removing the
#     alleles deemed "impossible" under these phased conditions. Details of removed
#     allele pairs due to subsumption are stored in
#     `bg.filtered_out["allele_subsumed_by_other_phased"]`.

#     Args:
#         bg (BloodGroup): A BloodGroup object containing allele states, variant pool
#             (with observed zygosities), and phasing information within Allele objects.
#         phased (bool): A flag indicating whether phasing rules should be applied.
#             If False, the function returns the BloodGroup object unmodified.

#     Returns:
#         BloodGroup: The modified BloodGroup object. The `alleles[AlleleState.POS]`
#             list may be reduced, and `filtered_out` may be updated.


#     Example (same phase mainly HET):
#     GYPB*03/GYPB*03N.03 - in GYPB*03N.04
#     GYPB*03/GYPB*03N.04 - is the only posibility as all vars in phase
#     GYPB*03/GYPB*06.02 - in GYPB*03N.04
#     Phenotypes (numeric): MNS:3,-4,5 | MNS:3,-4,6
#     Phenotypes (alphanumeric): S+,s-,He+ | S+,s-,U+

#     #Data:
#     Vars: {
#     '4:143999443_G_A': 'Homozygous',
#     '4:144001261_T_C': 'Heterozygous',
#     '4:144001254_T_A': 'Heterozygous',
#     '4:144001250_T_C': 'Heterozygous',
#     '4:144001249_C_A': 'Heterozygous',
#     '4:144001262_A_C': 'Heterozygous',
#     '4:143997535_C_A': 'Heterozygous'}

#     Filtered:
#     Allele
#     genotype: GYPB*03
#     defining_variants:
#             4:143999443_G_A #hom
#     weight_geno: 1000
#     phenotype: MNS:3,-4 or S+,s-
#     reference: False
#     phases: ('.',)

#     Allele
#     genotype: GYPB*06.02
#     defining_variants:
#             4:144001249_C_A
#             4:143999443_G_A #hom
#             4:144001254_T_A
#             4:144001262_A_C
#             4:144001261_T_C
#             4:144001250_T_C
#     weight_geno: 1000
#     phenotype: MNS:3,-4,6 or S+,s-,He+
#     reference: False
#     phases: ('143997535', '143997535', '143997535',
#       '143997535', '143997535', '.')

#     Allele
#     genotype: GYPB*03N.03
#     defining_variants:
#             4:143997535_C_A
#             4:143999443_G_A #hom
#     weight_geno: 1000
#     phenotype: MNS:-3,-4,5w or S-,s-,U+w
#     reference: False
#     phases: ('.', '143997535')

#     Allele
#     genotype: GYPB*03N.04
#     defining_variants:
#             4:143997535_C_A
#             4:144001249_C_A
#             4:143999443_G_A #hom
#             4:144001254_T_A
#             4:144001262_A_C
#             4:144001261_T_C
#             4:144001250_T_C
#     weight_geno: 1000
#     phenotype: MNS:-3,-4,5w or S-,s-,U+w
#     reference: False
#     phases: ('143997535', '143997535', '143997535',
#       '143997535', '143997535', '143997535', '.')


#     """

#     def check_hom(current_allele: Allele) -> bool:
#         """
#         Checks if all elements in a list are hom.

#         Args:
#             current_allele: Allele

#         Returns:
#             True if all elements match the target_str, or if all elements
#             except exactly one match the target_str. False otherwise.
#         """
#         variants = [
#             zygo
#             for variant, zygo in bg.variant_pool.items()
#             if variant in current_allele.defining_variants
#         ]

#         return all(item == Zygosity.HOM for item in variants)

#     def in_all(current_allele: Allele) -> bool:
#         """
#         check if an allele is in all other alleles that are either in
#         the same phase or all HOM, or a mix of HOM and same phase

#         example:

#                 Vars: {
#         '7:142957921_ref': 'Homozygous',
#         '7:142944339_C_T': 'Heterozygous',
#         '7:142943026_A_G': 'Heterozygous'}
#         Vars_phase: {
#         '7:142957921_ref': '.',
#         '7:142944339_C_T': '0|1',
#         '7:142943026_A_G': '1|0'}
#         Vars_phase_set: {
#         '7:142957921_ref': '.',
#         '7:142944339_C_T': '139403281',
#         '7:142943026_A_G': '139403281'}
#         Raw:
#             Allele
#             genotype: KEL*02 - is in both other options
#             defining_variants:
#                     7:142957921_ref
#             weight_geno: 1000
#             phenotype: KEL:-1,2 or K-,k+
#             reference: True

#             Allele
#             genotype: KEL*02.-19
#             defining_variants:
#                     7:142957921_ref
#                     7:142944339_C_T
#             weight_geno: 1000
#             phenotype: KEL:-1,2,-19 or K-,k+,K19-
#             reference: False

#             Allele
#             genotype: KEL*02.06
#             defining_variants:
#                     7:142943026_A_G
#                     7:142957921_ref
#             weight_geno: 1000
#             phenotype: KEL:-1,2,6,-7 or K-,k+,Js(a+),Js(b-)
#             reference: False
#         """
#         if len(alleles) <= 2:
#             return False  # need all
#         else:
#             counter = 0
#             for other_allele in alleles:
#                 if other_allele == current_allele:
#                     continue
#                 if current_allele in other_allele:
#                     counter += 1

#             return counter == len(alleles) - 1

#     #ic(789, bg.alleles[AlleleState.NORMAL])
#     alleles = list(flatten_alleles(bg.alleles[AlleleState.NORMAL]))

#     to_test = "JK"
#     if phased:
#         impossible_alleles = set([])
#         impossible_pairs = set([])
#         for allele in alleles:
#             if bg.type == to_test:
#                 ic(11111, allele)
#             if allele in impossible_alleles:
#                 if bg.type == to_test:
#                     ic(22222)
#                 continue
#             # if allele.reference:
#             #     if bg.type == to_test:
#             #         ic(22222333333333)
#             #     continue

#             if check_hom(allele) and not in_all(allele):
#                 if bg.type == to_test:
#                     ic(3333333)
#                 continue
#             for allele2 in alleles:
#                 if bg.type == to_test:
#                     ic(44444444, allele2)
#                 if allele == allele2:
#                     if bg.type == to_test:
#                         ic(555555)
#                     continue
#                 # if allele2.reference:
#                 #     continue
#                 if check_hom(allele2) and not in_all(allele2):
#                     if bg.type == to_test:
#                         ic(6666666)
#                     continue
#                 if allele in allele2:
#                     if bg.type == to_test:
#                         ic(7777777, allele.genotype, allele2.genotype)
#                     impossible_alleles.add(allele)
#                     impossible_pairs.add(Pair(allele1=allele2, allele2=allele))
#                     continue
#                 if allele2 in allele:
#                     if bg.type == to_test:
#                         ic(888888, allele.genotype, allele2.genotype)
#                     impossible_alleles.add(allele2)
#                     impossible_pairs.add(Pair(allele1=allele2, allele2=allele))
#                     continue
#         if bg.type == to_test:
#             ic(9999999, impossible_alleles)

#         bg.filtered_out["allele1_not_possible_due_to_allele2_phased"] = list(
#             impossible_pairs
#         )
#         if bg.type == to_test:
#             ic(111111122222, bg.alleles[AlleleState.POS])
#         to_remove = []
#         for pair in bg.alleles[AlleleState.NORMAL]:
#             if bg.type == to_test:
#                 ic(1111113333, pair, impossible_alleles)
#             if pair.allele1 in impossible_alleles or pair.allele2 in impossible_alleles:
#                 to_remove.append(pair)

#         bg.remove_pairs(to_remove, "filter_impossible_alleles_phased")

#     return bg
