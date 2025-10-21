import itertools
import operator
from collections import defaultdict
from dataclasses import dataclass
from functools import partial
from typing import Any, Protocol
from loguru import logger

from rbceq2.core_logic.alleles import Allele, BloodGroup, Pair
from rbceq2.core_logic.constants import AlleleState
from rbceq2.core_logic.utils import (
    Zygosity,
    apply_to_dict_values,
    check_available_variants,
    chunk_geno_list_by_rank,
    get_non_refs,
    collapse_variant,
)
from rbceq2.db.db import Db
from rbceq2.IO.vcf import VCF
import pandas as pd
from icecream import ic
from rbceq2.core_logic.large_variants import select_best_per_vcf, MatchResult


def raw_results(
    db: Db, vcf: VCF, exclude: list[str], var_map, matches
) -> dict[str, list[Allele]]:
    """Generate raw results from database alleles and VCF data based on phasing
    information.

    Args:
        db (Db): The database containing allele definitions and methods to generate
        them.
        vcf (VCF): The VCF data containing variants and possibly phasing information.

    Returns:
        Dict[str, List[Allele]]: A dictionary mapping blood groups to lists of Allele
        objects.

    best = select_best_per_vcf(matches, tie_tol=1e-9)
    best_var_map={}
    for match in best:
        best_var_map[f"{match.vcf.chrom}:{match.db.raw}"] = match.variant
    alleles_with_big_vars = []
    for allele in bg.alleles[AlleleState.RAW]:
        if any(var in best_var_map for var in allele.defining_variants):
            allele = allele.with_big_variants(best_var_map)
        alleles_with_big_vars.append(allele)
    bg.alleles[AlleleState.RAW] = alleles_with_big_vars
    """
    res: dict[str, list[Allele]] = defaultdict(list)
    for allele in db.make_alleles():
        if any(x in allele.genotype for x in exclude):
            continue

        if all(variant in vcf.variants for variant in allele.defining_variants):
            if any(variant in var_map for variant in allele.defining_variants):
                # best_match = [
                #     match for match in select_best_per_vcf(matches, tie_tol=1e-9) if
                #       allele.blood_group in match.db.id
                #       ]
                # #ic(allele.blood_group,best_match, select_best_per_vcf(matches, tie_tol=1e-9))
                # #ic(allele,allele.big_variants, var_map, best_match)
                # if len(best_match) == 0:
                #     continue
                # assert len(best_match) ==1
                # best = best_match[0]
                # best_var_map = {db_var: vcf_var for db_var, vcf_var in var_map.items() if
                #                 db_var == f"{best.vcf.chrom}:{best.db.raw}"}
                # #ic(allele, best_var_map,allele.big_variants, var_map, best_match)
                allele = allele.with_big_variants(var_map)
            res[allele.blood_group].append(allele)
            # ic(allele.blood_group, allele)
    return res


def make_blood_groups(
    res: dict[str, list[Allele]], sample: str
) -> dict[str, BloodGroup]:
    """Create a dictionary of BloodGroup objects from allele data.

    Iterates through the 'res' mapping of blood group identifiers to lists of Allele
    objects, and constructs a new dictionary where each key is a blood group name and
    each value is a BloodGroup instance.

    Args:
        res (dict[str, list[Allele]]):
            A dictionary mapping blood group names to a list of Allele objects.
        sample (str):
            The sample identifier to be associated with each BloodGroup.

    Returns:
        dict[str, BloodGroup]:
            A dictionary mapping blood group identifiers to BloodGroup instances.
    """
    new_dict: dict[str, BloodGroup] = {}
    for blood_group, alleles in res.items():
        new_dict[blood_group] = BloodGroup(
            type=blood_group, alleles={AlleleState.RAW: alleles}, sample=sample
        )

    return new_dict


@apply_to_dict_values
def add_phasing(
    bg: BloodGroup,
    phased: bool,
    variant_metrics: dict[str, dict[str, str]],
    phase_sets: dict[str, dict[int, tuple[int, int]]],
) -> BloodGroup:
    """Add phasing information to alleles in a BloodGroup.

    The Biological Scenario (gemini 2.5 pro)
    Imagine a gene called GENEX. In an individual, we find two different
    heterozygous variants within this gene:
    At position chr1:1,000,000, there's a A > G variant.
    At position chr1:1,002,500, there's a C > T variant.
    Since this person is diploid, they have two copies of GENEX (one
    maternal, one paternal). The critical question that phasing answers
    is: Are the two alternate alleles (G and T) on the same copy of the
    chromosome, or on different copies?
    There are two possibilities:
    Scenario A (Cis configuration): The alternate alleles are on the same chromosome.
    Maternal Chromosome: ---G---T---
    Paternal Chromosome: ---A---C--- (reference alleles)
    Scenario B (Trans configuration): The alternate alleles are on opposite chromosomes.
    Maternal Chromosome: ---G---C---
    Paternal Chromosome: ---A---T---

    How This Looks in the VCF
    A phasing algorithm (using long reads, parental data, or statistical inference)
    will try to determine whether we are in Scenario A or B. If it succeeds,
    it will report both variants in the same phase set.
    Let's assign the PS (Phase Set) identifier 12345 to this phased block.
    VCF for Scenario A (Cis):
    The two alternate alleles (1) are on the same haplotype.
    Generated vcf
    #CHROM  POS        ID  REF ALT ... FORMAT     SAMPLE1
    chr1    1000000    .   A   G   ... GT:PS      0|1:12345
    chr1    1002500    .   C   T   ... GT:PS      0|1:12345
    Vcf
    PS:12345: Both variants are in the same phase set.
    0|1 and 0|1: This tells us that the allele designated '1' (the ALT allele)
    for the first variant is on the same chromosome as the allele designated
    '1' for the second variant. Haplotype 1 is G-T, Haplotype 0 is A-C.
    VCF for Scenario B (Trans):
    The two alternate alleles (1) are on opposite haplotypes.
    Generated vcf
    #CHROM  POS        ID  REF ALT ... FORMAT     SAMPLE1
    chr1    1000000    .   A   G   ... GT:PS      0|1:12345
    chr1    1002500    .   C   T   ... GT:PS      1|0:12345
    Vcf

    PS:12345: They are still in the same phase set. This is the key point. The
    PS just tells you "the relationship between these variants is known."

    0|1 and 1|0: The GT fields now describe the trans relationship. For the first
    variant, the alternate allele (G) is on Haplotype 1. For the second variant,
    the alternate allele (T) is on Haplotype 0. Haplotype 1 is G-C, Haplotype
    0 is A-T.

    When Would They Be in Different Phase Sets?
    Variants within the same gene would only appear in different phase sets if
    the phasing process failed to connect them.


    If 'phased' is True, updates the alleles in the given 'bg' object by assigning
    phase sets (PS) from 'variant_metrics'. Alleles containing reference-only variants
    are ignored or reduced in count.

    Args:
        bg (BloodGroup):
            The BloodGroup object whose alleles are to be phased.
        phased (bool):
            A flag indicating whether phasing is enabled or not.
        variant_metrics (dict[str, dict[str, str]]):
            A nested dictionary of metrics for each variant. The inner dictionary
            should include the phase set (PS) or other variant-related metrics.

    Returns:
        BloodGroup:
            The updated BloodGroup object with phased alleles, if 'phased' is True.
            Otherwise, the original BloodGroup object.

    Raises:
        AssertionError:
            If the number of updated alleles does not match the number of original alleles,
            or if the computed phases do not match the expected length based on 'defining_variants'.
    """

    def assign_ref_phase(
        current_variant: str,
        # bg_variant_pool: dict[str, str],
        # current_phase_pool: dict[str, str],
    ) -> str:
        """
        Assigns the correct phase to a reference variant at a heterozygous site.
        """
        zygosity = bg.variant_pool.get(current_variant)
        if zygosity == Zygosity.HOM:
            return "1/1"
        elif "/" in phase_pool[current_variant]:
            return "unknown"
        return phase_pool[current_variant]
        # # 1. Validate that the input variant is a heterozygous reference variant
        # if not current_variant.endswith("_ref"):
        #     raise ValueError(
        #         f"Input variant '{current_variant}' must be a reference variant (ending in '_ref')."
        #     )

        # if zygosity != Zygosity.HET:
        #     raise ValueError(
        #         f"Variant '{current_variant}' must be present and 'Heterozygous' in bg_variant_pool. "
        #         f"Found: {zygosity}"
        #     )

        # # 2. Find the corresponding alternate allele variant at the same position
        # position = current_variant.split("_")[0]
        # partner_variant = None
        # for key in phase_pool:
        #     # Match keys that start with the same position but are not the ref variant itself
        #     if key.startswith(position + "_") and key != current_variant:
        #         partner_variant = key
        #         break  # Found the partner, no need to continue searching
        # if partner_variant is None:
        #     return phase_pool[current_variant]

        # # 3. Get the phase of the partner variant
        # partner_phase = phase_pool.get(partner_variant)
        # if partner_phase is None:
        #     raise ValueError(
        #         f"Partner variant '{partner_variant}' found but has no entry in phase_pool."
        #     )

        # # 4. Calculate the ref phase by "flipping" the partner's phase
        # phase_parts = partner_phase.split("|")
        # if len(phase_parts) != 2:
        #     raise ValueError(
        #         f"Invalid phase format for partner '{partner_variant}': '{partner_phase}'"
        #     )

        # # The ref allele is on the opposite haplotype from the alt allele
        # # e.g., if partner is 1|0, ref is 0|1
        # ref_phase = f"{phase_parts[1]}|{phase_parts[0]}"
        # ic(888888888888888888,current_variant, ref_phase)
        # return ref_phase

    def assign_ref_phase_set(current_variant):
        """ """

        def get_phase_set_for_loci() -> str:
            """Queries for a phase set ID given a locus.

            Args:
                loci (str): A locus in "CHROM:POS" format (e.g., "1:12345"
                            or "chr1:12345").

            Returns:
                int | None: The phase set ID if the locus falls within a known
                            phased block, otherwise None.
            """
            chrom, pos_str = current_variant.split("_")[0].split(":")
            pos = int(pos_str)

            # Normalize chromosome name to match internal representation (e.g. '1' not 'chr1')
            chrom = chrom.replace("chr", "")

            # Get all phase sets for the given chromosome
            chrom_phase_sets = phase_sets.get(chrom)
            if not chrom_phase_sets:
                return "unknown"

            # Check if the position falls within any of the phase set ranges
            for ps_id, (min_pos, max_pos) in chrom_phase_sets.items():
                if min_pos <= pos <= max_pos:
                    return str(ps_id)

            # If no matching phase set is found
            return "unknown"

        zygosity = bg.variant_pool.get(current_variant)
        if zygosity == Zygosity.HOM:
            return "."
        # 2. Find the corresponding alternate allele variant at the same position
        position = current_variant.split("_")[0]
        partner_variant = None
        for key in phase_set_pool:
            # Match keys that start with the same position but are not the ref variant itself
            if key.startswith(position + "_") and key != current_variant:
                partner_variant = key
                break  # Found the partner, no need to continue searching
        if partner_variant is None:
            return get_phase_set_for_loci()
        partner = phase_set_pool.get(partner_variant)
        if partner is not None:
            return partner
        else:
            return get_phase_set_for_loci()

    if phased:  # TODO enum for GT, PS etc
        phase_pool = {
            variant: variant_metrics[variant]["GT"] for variant in bg.variant_pool
        }
        phase_set_pool = {
            variant: variant_metrics[variant].get("PS") for variant in bg.variant_pool
        }

        phase_pool_ref_fixed = {}
        for variant, phase in phase_pool.items():
            if variant.endswith("_ref"):
                new_phase = assign_ref_phase(variant)
                phase_pool_ref_fixed[variant] = new_phase
            else:
                phase_pool_ref_fixed[variant] = phase
        bg.variant_pool_phase = phase_pool_ref_fixed

        phase_set_pool_ref_fixed = {}
        for variant, phase in phase_set_pool.items():
            if variant.endswith("_ref") or phase is None:
                new_phase = assign_ref_phase_set(variant)
                phase_set_pool_ref_fixed[variant] = new_phase
            else:
                phase_set_pool_ref_fixed[variant] = phase

        bg.variant_pool_phase_set = phase_set_pool_ref_fixed
    # if bg.type == 'RHD':
    #     ic(bg.variant_pool_phase, bg.variant_pool_phase_set)
    return bg


@apply_to_dict_values
def ABO_phasing(
    bg: BloodGroup,
    phased: bool,
) -> BloodGroup:
    # aboO_phases2 = set([]) I'm opting out of this path but if it ever gets revisited
    # do it at DB level #TODO look for vars that have only ever been observed in cis with
    # ABO*O, or never been observed with it and take phasing info from there. Best
    # to actually phase indels though
    """9:133257521(GRCh38) and 9:136132908 (GRCh37) _T_TC or _ref are pivotal for ABO
    calls but aren't always assigned a phase group by phasing algos

    some alleles are same/similar except for this locus

    ie HG04054 for 1kg data

    Allele
    genotype: ABO*B.01
    defining_variants:
                9:133255935_G_T 0|1
                9:133257521_T_TC 0/1
                9:133257486_T_C 1/1
                9:133255928_C_G 0|1
                9:133256074_G_A 0|1
                9:133256028_C_T 0|1
                9:133256205_G_C 0|1
                9:133255801_C_T 0|1
    weight_geno: 1000
    phenotype: . or B
    reference: False

    Allele
    genotype: ABO*O.01.41
    defining_variants:
                9:133257486_T_C 1/1
                9:133255928_C_G 0|1
                9:133256074_G_A 0|1
                9:133256028_C_T 0|1
                9:133256205_G_C 0|1
                9:133257521_ref unknown
                9:133255801_C_T 0|1

    this function will infer the phase of this locus from the phase of ABO*O.01
    specific variants
    """
    if not phased:
        return bg

    if bg.type != "ABO":
        return bg
    # 261delG
    c261delGs = [
        "9:133257521_ref",
        "9:133257521_T_TC",
        "9:136132908_ref",
        "9:136132908_T_TC",
    ]
    if any(bg.variant_pool.get(c261delG) == Zygosity.HOM for c261delG in c261delGs):
        return bg

    aboO = set([])
    other = set([])
    for allele in bg.alleles[AlleleState.FILT]:
        if allele.genotype.startswith("ABO*O.01"):
            for variant in allele.defining_variants:
                aboO.add(variant)
        else:
            for variant in allele.defining_variants:
                other.add(variant)
    aboO_phases = set([])
    for variant in aboO.difference(other):
        if not variant.startswith(("9:133257521", "9:136132908")):
            phase = bg.variant_pool_phase[variant]
            aboO_phases.add(phase)

    if len(aboO_phases) == 0:
        return bg  # can't rescue ABO
    if len(aboO_phases) > 1:
        return bg  # can't rescue ABO
    abo_phase = aboO_phases.pop()
    not_abo_phase = "1|0" if abo_phase == "0|1" else "0|1"
    new_phases = {}

    for variant, phase in bg.variant_pool_phase.items():
        if variant in c261delGs:
            if variant.endswith("_ref"):
                new_phases[variant] = abo_phase
            else:
                new_phases[variant] = not_abo_phase
        else:
            new_phases[variant] = phase
    bg.variant_pool_phase = new_phases

    return bg


@apply_to_dict_values
def make_variant_pool(bg: BloodGroup, vcf: VCF) -> BloodGroup:
    """Construct or update a variant pool for a BloodGroup from VCF data.

    This function traverses the alleles in the BloodGroup object, extracts reference
    information for each defining variant from the VCF, and combines these into a
    single dictionary (the variant pool).

    Args:
        bg (BloodGroup):
            The BloodGroup object to be updated with the new variant pool.
        vcf (VCF):
            The VCF object providing variant data.

    Returns:
        BloodGroup:
            The updated BloodGroup object, including the combined 'variant_pool' with
            reference data for each defining variant.

    Raises:
        KeyError:
            If a variant in 'bg.alleles[AlleleState.FILT]' is not found in 'vcf.variants'.
    """
    variant_pool = {}

    for allele in bg.alleles[AlleleState.FILT]:
        zygosity = {var: get_ref(vcf.variants[var]) for var in allele.defining_variants}
        variant_pool = variant_pool | zygosity
    bg.variant_pool = variant_pool

    return bg


def _modify_variant_pool_with_large_indel(
    variant_pool: dict, sample: str, bg_type: str, is_phase_pool: bool = False
) -> dict:
    """Internal helper to adjust variant zygosity when large deletions are present.

    Args:
        variant_pool (dict): Dictionary mapping variant strings to zygosity values.
        sample (str): Sample identifier for logging.
        bg_type (str): Blood group type for logging.
        is_phase_pool (bool): If True, expects phase notation ('1/1', '1'),
            otherwise Zygosity enums. Defaults to False.

    Returns:
        dict: Modified variant pool dictionary, or empty dict if no modifications needed.
    """

    def get_start_pos(current_variant):
        return int(current_variant.strip().split(":")[1].split("_")[0])

    # Identify large deletions
    big_dels = []
    for variant in variant_pool:
        no_seq_variant = collapse_variant(variant)
        if "DEL" in no_seq_variant.upper():
            big_dels.append(
                (variant, no_seq_variant)
            )  # these are teh db version of var,
            # should try get the VCF version TODO

    # Early exit if no large deletions or only deletions present
    if not big_dels or len(variant_pool) <= len(big_dels):
        return {}

    # Determine zygosity values based on pool type
    if is_phase_pool:
        hom_value = "1/1"
        het_value = ("1|0", "0|1")
        hem_value = "1"
    else:
        hom_value = Zygosity.HOM
        het_value = Zygosity.HET
        hem_value = Zygosity.HEM

    new_variant_pool = {}

    for big_del, no_seq_variant in big_dels:
        start = get_start_pos(no_seq_variant)
        length = no_seq_variant.split("_")[-1]
        length = int(length[:-2]) * 1000 if length.endswith("kb") else int(length)
        end = start + length

        # Check if deletion is homozygous
        big_del_is_hom = variant_pool.get(big_del) == hom_value

        for variant, zygosity in variant_pool.items():
            if variant == big_del:
                # Keep deletion as-is
                new_variant_pool[variant] = zygosity
            elif start < get_start_pos(variant) < end:
                # Variant falls within deletion range
                if big_del_is_hom:
                    assert variant.endswith("_ref")
                    new_variant_pool[variant] = zygosity
                else:
                    if zygosity == hom_value:
                        # Convert homozygous to hemizygous
                        new_variant_pool[variant] = hem_value
                    else:
                        # Warn about unexpected heterozygous variant
                        new_variant_pool[variant] = zygosity
                        logger.warning(
                            f"Heterozygous variant detected where hemizygousity expected "
                            f"{sample} {bg_type} {variant}"
                        )
                    if not is_phase_pool:
                        assert variant_pool[big_del] in het_value
            else:
                # Variant outside deletion range
                new_variant_pool[variant] = zygosity

    return new_variant_pool


@apply_to_dict_values
def modify_variant_pool_if_large_indel(bg: BloodGroup) -> BloodGroup:
    """Adjusts variant zygosity calls when large deletions are present.

    When a blood group contains large deletions alongside other variants, this function
    modifies the zygosity of variants that fall within the deletion boundaries. Variants
    located inside a deletion region are reset to hemizygous (HEM) status when homozygous,
    as they cannot reliably be called homozygous when overlapped by a heterozygous deletion.

    The function only processes variant pools that contain both large deletions and
    additional variants (i.e., pools with more than just deletions).

    Args:
        bg (BloodGroup): A BloodGroup object containing a variant_pool dictionary
            mapping variant strings to Zygosity enum values.

    Returns:
        BloodGroup: The modified BloodGroup object with updated variant_pool. If no
            large deletions are found or if the pool contains only deletions, returns
            the original BloodGroup unchanged.

    Note:
        - Variant strings are expected in format: "chr:pos_ref_alt" or "chr:pos_...DEL...kb"
        - Large deletions are identified by the presence of 'DEL' in the collapsed variant string
        - Deletion sizes ending in 'kb' are converted to base pairs (multiplied by 1000)
        - The function asserts that all large deletions in the pool are heterozygous, if there
          are other small HET variants that overlap the range

    Example:
        bg.variant_pool: {'4:143995187_del_103kb': 'Heterozygous',
                          '4:143999443_ref': 'Homozygous'}
        Becomes:
        bg.variant_pool: {'4:143995187_del_103kb': 'Heterozygous',
                          '4:143999443_ref': 'Hemizygous'}
    """
    new_variant_pool = _modify_variant_pool_with_large_indel(
        bg.variant_pool, bg.sample, bg.type, is_phase_pool=False
    )

    if new_variant_pool:
        bg.variant_pool = new_variant_pool

    return bg


@apply_to_dict_values
def modify_variant_phase_pool_if_large_indel(bg: BloodGroup) -> BloodGroup:
    """Adjusts variant phase calls when large deletions are present.

    When a blood group contains large deletions alongside other variants, this function
    modifies the phase notation of variants that fall within the deletion boundaries.
    Variants located inside a deletion region are converted from '1/1' to '1' (hemizygous)
    when overlapped by a heterozygous deletion.

    The function only processes variant pools that contain both large deletions and
    additional variants (i.e., pools with more than just deletions).

    Args:
        bg (BloodGroup): A BloodGroup object containing a variant_pool_phase dictionary
            mapping variant strings to phase notation strings.

    Returns:
        BloodGroup: The modified BloodGroup object with updated variant_pool_phase. If
            no large deletions are found or if the pool contains only deletions, returns
            the original BloodGroup unchanged.

    Note:
        - Variant strings are expected in format: "chr:pos_ref_alt" or "chr:pos_...DEL...kb"
        - Large deletions are identified by the presence of 'DEL' in the collapsed variant string
        - Deletion sizes ending in 'kb' are converted to base pairs (multiplied by 1000)

    Example:
        bg.variant_pool_phase: {'4:143995187_del_103kb': '0/1',
                                 '4:143999443_ref': '1/1'}
        Becomes:
        bg.variant_pool_phase: {'4:143995187_del_103kb': '0/1',
                                 '4:143999443_ref': '1'}
    """
    new_variant_pool = _modify_variant_pool_with_large_indel(
        bg.variant_pool_phase, bg.sample, bg.type, is_phase_pool=True
    )

    if new_variant_pool:
        bg.variant_pool_phase = new_variant_pool

    return bg


@apply_to_dict_values
def modify_allele_pool_if_large_indel(bg: BloodGroup) -> BloodGroup:
    """removes alleles that can't exist due to big indel

    Example - HOM del
    allele1: Allele
              genotype: RHD*01N.01
              defining_variants:
                        1:25272547_DEL_59419
              weight_geno: 500
              phenotype: . or D-
              reference: False
    Removed as '_ref' is over ridden by big del
    allele2: Allele
              genotype: RHD*10.00
              defining_variants:
                        1:25317062_ref
              weight_geno: 1000
              phenotype: . or DAU0
              reference: False

    """
    keepers = []
    for allele in bg.alleles[AlleleState.FILT]:
        for variant in allele.defining_variants:
            try:
                assert variant in bg.variant_pool
            except AssertionError:
                ic(bg.sample, allele, variant, bg.variant_pool)
        if all(variant in bg.variant_pool for variant in allele.defining_variants):
            keepers.append(allele)
    bg.alleles[AlleleState.FILT] = keepers
    return bg


@apply_to_dict_values
def modify_phase_of_large_indel(bg: BloodGroup, phased: bool) -> BloodGroup:
    """Infers and updates phase information for large deletions and reference variants.

    When a large deletion is unphased (has '/' separator) but overlaps with phased variants
    (using '|' separator) in the same phase set, this function infers the deletion's phase
    as the opposite haplotype of the overlapping variants. This is valid because a deletion
    on one haplotype means the other haplotype carries the reference allele where variants
    are called.

    Additionally, reference variants (ending in '_ref') have their phase inferred based on
    which deletion they overlap with. If a reference variant overlaps with a deletion
    that is phased as '0|1', the reference is phased as '1|0' (opposite haplotype).

    The function only modifies deletion phase if:
    1. The deletion overlaps with at least one phased variant
    2. All variants (including those ending in '_ref') share the same phase set
    3. The deletion is currently unphased (contains '/')

    Args:
        bg: A BloodGroup object containing variant_pool_phase (dict mapping variants to
            phase strings like '0|1' or '0/1') and variant_pool_phase_set (dict mapping
            variants to phase set IDs).
        phased: Boolean indicating whether phasing information is available. If False,
            returns the BloodGroup unchanged.

    Returns:
        The modified BloodGroup object with updated phase information for large deletions,
        reference variants, and phase sets. If no phasing is available or no large
        deletions need updating, returns the original BloodGroup unchanged.

    Raises:
        ValueError: If a reference variant overlaps with multiple deletions, which is
            biologically impossible (only 2 alleles exist per locus).

    Note:
        - Reference variants (ending in '_ref') ARE included in phase set validation
        - Phase is flipped: if overlapping variant is '1|0', deletion becomes '0|1'
        - Deletions are identified by 'del' or 'DEL' in the variant string
        - Only deletions with '/' in their phase string are considered for updating
        - Multiple non-overlapping deletions can have different phases (different chromosomes)
        - Any variant overlapping a deletion must be on the other chromosome

    Example:
        >>> bg = BloodGroup(
        ...     variant_pool_phase={
        ...         '4:143995187_del_103kb': '0/1',
        ...         '4:143997559_C_G': '1|0',
        ...         '4:143999443_ref': 'unknown'
        ...     },
        ...     variant_pool_phase_set={
        ...         '4:143995187_del_103kb': '143821229',
        ...         '4:143997559_C_G': '143821229',
        ...         '4:143999443_ref': 'unknown'
        ...     }
        ... )
        >>> modified_bg = modify_phase_if_large_indel(bg, phased=True)
        >>> modified_bg.variant_pool_phase['4:143995187_del_103kb']
        '0|1'
        >>> modified_bg.variant_pool_phase_set['4:143999443_ref']
        '143821229'
        >>> modified_bg.variant_pool_phase['4:143999443_ref']
        '1|0'
    """

    def get_start_pos(current_variant):
        return int(current_variant.strip().split(":")[1].split("_")[0])

    def flip_phase(phase_str):
        """Flip the phase: '1|0' becomes '0|1' and vice versa."""
        alleles = phase_str.split("|")
        return f"{alleles[1]}|{alleles[0]}"

    def get_deletion_boundaries(big_del):
        """Extract start and end positions for a deletion."""
        start = get_start_pos(big_del)
        length_str = big_del.split("_")[-1]
        length = (
            int(length_str[:-2]) * 1000
            if length_str.endswith("kb")
            else int(length_str)
        )
        end = start + length
        return start, end

    if not phased:
        return bg

    # Find large deletions that are unphased
    unphased_big_dels = []
    for variant, phase in bg.variant_pool_phase.items():
        if ("del" in variant.lower() or "DEL" in variant) and "/" in phase:
            unphased_big_dels.append(variant)

    if not unphased_big_dels:
        return bg

    # Collect all phase sets (excluding 'unknown' and '.')
    all_phase_sets = set()
    for variant, phase_set in bg.variant_pool_phase_set.items():
        if phase_set not in ["unknown", "."]:
            all_phase_sets.add(phase_set)

    # Only proceed if all variants share the same phase set
    if len(all_phase_sets) != 1:
        return bg

    common_phase_set = list(all_phase_sets)[0]

    # Process each unphased deletion
    for big_del in unphased_big_dels:
        start, end = get_deletion_boundaries(big_del)

        # Find overlapping phased variants
        overlapping_variants = []

        for variant, phase in bg.variant_pool_phase.items():
            if variant == big_del:
                continue

            if "|" in phase:  # Only consider phased variants
                variant_pos = get_start_pos(variant)
                if start < variant_pos < end:
                    overlapping_variants.append(variant)

        # Update deletion phase if we have overlapping variants
        if overlapping_variants:
            # Use the first overlapping variant's phase to infer deletion phase
            overlapping_phase = bg.variant_pool_phase[overlapping_variants[0]]
            inferred_del_phase = flip_phase(overlapping_phase)

            # Update the deletion's phase and phase set
            bg.variant_pool_phase[big_del] = inferred_del_phase
            bg.variant_pool_phase_set[big_del] = common_phase_set

    # Now infer phase for reference variants based on overlapping deletions
    # First, collect all deletions (both originally phased and newly phased)
    all_deletions = {}
    for variant, phase in bg.variant_pool_phase.items():
        if ("del" in variant.lower() or "DEL" in variant) and "|" in phase:
            start, end = get_deletion_boundaries(variant)
            all_deletions[variant] = {"phase": phase, "start": start, "end": end}

    # Update reference variants
    for variant in list(bg.variant_pool_phase.keys()):
        if variant.endswith("_ref"):
            variant_pos = get_start_pos(variant)

            # Find which deletion(s) this reference overlaps with
            overlapping_deletions = []
            for del_variant, del_info in all_deletions.items():
                if del_info["start"] < variant_pos < del_info["end"]:
                    overlapping_deletions.append(del_variant)

            # Check for invalid case: multiple deletions overlapping same locus
            if len(overlapping_deletions) > 1:
                raise ValueError(
                    f"Reference variant {variant} at position {variant_pos} overlaps with "
                    f"multiple deletions: {overlapping_deletions}. This is biologically "
                    f"impossible as there are only 2 alleles per locus."
                )

            # If the reference overlaps with a deletion, infer its phase
            if overlapping_deletions:
                # Reference is on the opposite haplotype from the deletion
                del_phase = all_deletions[overlapping_deletions[0]]["phase"]
                ref_phase = flip_phase(del_phase)
                bg.variant_pool_phase[variant] = ref_phase
                bg.variant_pool_phase_set[variant] = common_phase_set
            else:
                # If no overlap with deletions, just update phase set
                bg.variant_pool_phase_set[variant] = common_phase_set

    return bg


def get_ref(ref_dict: dict[str, str]) -> str:
    """Determine the zygosity from a reference dictionary containing genotype
    information.

    Args:
        ref_dict (Dict[str, str]): A dictionary containing the genotype ("GT") and
        possibly other information.

    Returns:
        str: A string indicating the zygosity as 'Homozygous' or 'Heterozygous'.

    Raises:
        ValueError: If the genotype string does not conform to the expected format.

    The genotype string is expected to be in the format '0/1', '0|1', etc., where
    the delimiter can be '/' or '|'.
    A genotype of '0/0' or '1/1', etc., where both alleles are the same, will return
    'Homozygous'.
    A genotype of '0/1', '1/0', etc., will return 'Heterozygous'.
    """
    # 0/1:41,47:88:99:1080,0,1068:0.534:99
    ref_str = ref_dict["GT"]
    assert len(ref_str) == 3

    ref_str = ref_str.replace(".", "0")
    if ref_str[0] == ref_str[2]:
        return Zygosity.HOM
    return Zygosity.HET


@apply_to_dict_values
def get_genotypes(bg: BloodGroup) -> BloodGroup:
    """Generate genotype combinations for a given blood group from allele pairs.

    Args:
        bg (BloodGroup): The blood group object containing alleles.

    Returns:
        BloodGroup: The blood group object with updated genotypes based on allele
        combinations.

    This function processes 'pairs' and 'co_existing' alleles to create sorted genotype
    strings.
    """

    def make_list_of_lists(alleles):
        return [pair.genotypes for pair in alleles]

    if bg.alleles[AlleleState.CO] is not None:
        bg.genotypes = [
            "/".join(sorted(co))
            for co in make_list_of_lists(bg.alleles[AlleleState.CO])
        ]
    else:
        bg.genotypes = [
            "/".join(sorted(normal_pair))
            for normal_pair in make_list_of_lists(bg.alleles[AlleleState.NORMAL])
        ]

    return bg


def filter_vcf_metrics(
    alleles: list[Allele],
    variant_metrics: dict[str, dict[str, str]],
    metric_name: str,
    metric_threshold: float,
    microarray: bool,
) -> tuple[defaultdict[str, list[Allele]], list[Allele]]:
    """Filter out alleles based on a specified read depth metric.

    Iterates through each allele's defining variants and compares the specified metric
    (e.g., "DP" for read depth) to a threshold value. For microarray data, the read depth
    is set to a constant value of 30.0. Alleles whose read depth falls below the threshold
    are collected in `filtered_out`; all others are placed in `passed_filtering`.

    Args:
        alleles (list[Allele]):
            A list of allele objects to be evaluated.
        variant_metrics (dict[str, dict[str, str]]):
            A nested dictionary where the key is a variant identifier and the value
            is a dictionary of metrics (e.g., {"DP": "45", ...}).
        metric_name (str):
            The name of the metric to evaluate (e.g., "DP" for read depth).
        metric_threshold (float):
            The threshold value for the chosen metric. Alleles below this value
            are excluded.
        microarray (bool):
            If True, overrides the chosen metric by setting read depth to 30.0.

    Returns:
        tuple[defaultdict[str, list[Allele]], list[Allele]]:
            A tuple containing two elements:
            1. `filtered_out`: A defaultdict where each key is
               "variant:read_depth" and each value is a list of alleles
               that did not meet the threshold.
            2. `passed_filtering`: A list of alleles that passed the threshold.

    Raises:
        KeyError:
            If a required variant or metric is missing in `variant_metrics`.
    """
    # TODO large dels will have depth zero
    filtered_out = defaultdict(list)
    passed_filtering = []
    metric_threshold = float(metric_threshold)
    for allele in alleles:
        if allele.big_variants:
            passed_filtering.append(allele)
            continue
        keep = True
        for variant in allele.defining_variants:
            read_depth = float(variant_metrics[variant][metric_name])
            if microarray:  # TODO !!
                read_depth = 30.0  # for microarray
            else:
                read_depth = float(variant_metrics[variant][metric_name])
            if read_depth < metric_threshold:
                filtered_out[f"{variant}:{str(read_depth)}"].append(allele)
                keep = False
        if keep:
            passed_filtering.append(allele)

    return filtered_out, passed_filtering


@apply_to_dict_values
def remove_alleles_with_low_read_depth(
    bg: BloodGroup,
    variant_metrics: dict[str, str],
    min_read_depth: int,
    microarray: bool,
) -> BloodGroup:
    """
    Remove alleles from a BloodGroup object that have defining variants with read depth
    below a specified minimum threshold.

    Args:
        bg (BloodGroup): The BloodGroup object containing alleles to filter.
        variant_metrics (dict[str, dict[str, int]]): A dictionary containing variant
        metrics with read depth information.
        min_read_depth (int): The minimum read depth threshold.

    Returns:
        BloodGroup: The BloodGroup object with alleles filtered based on read depth.
    """

    filtered_out, passed_filtering = filter_vcf_metrics(
        bg.alleles[AlleleState.FILT], variant_metrics, "DP", min_read_depth, microarray
    )
    if filtered_out:
        vars_affected = ",".join(filtered_out.keys())
        message = f"Read Depth. Sample: {bg.sample}, BG: {bg.type}, variant/s: {vars_affected}"
        logger.warning(message)
    bg.filtered_out["insufficient_read_depth"] = filtered_out
    bg.alleles[AlleleState.FILT] = passed_filtering
    return bg


@apply_to_dict_values
def only_keep_alleles_if_FILTER_PASS(
    bg: BloodGroup, df: pd.DataFrame, no_filter: bool
) -> BloodGroup:
    """
    Remove alleles from a BloodGroup object that have defining variants with read depth
    below a specified minimum threshold.

    Args:
        bg (BloodGroup): The BloodGroup object containing alleles to filter.
        variant_metrics (dict[str, dict[str, int]]): A dictionary containing variant
        metrics with read depth information.
        min_read_depth (int): The minimum read depth threshold.

    Returns:
        BloodGroup: The BloodGroup object with alleles filtered based on read depth.
    """
    if no_filter:
        bg.alleles[AlleleState.FILT] = bg.alleles[AlleleState.RAW]
        return bg
    passed_filtering = []
    for allele in bg.alleles[AlleleState.RAW]:
        keeper = True
        for variant in allele.defining_variants:
            if "_ref" in variant:
                continue
            vcf_var = allele.big_variants.get(variant, variant)
            loci = vcf_var.split("_")[0]
            try:
                filter_value = df.query("variant.str.contains(@vcf_var)")[
                    "FILTER"
                ].iloc[0]
            except IndexError:
                message = f"FILTER parsing failed. Sample: {bg.sample}, BG: {bg.type}, variant/s: {variant}"
                logger.error(message)
                ic(
                    bg.sample,
                    bg.type,
                    loci,
                    allele,
                    allele.big_variants,
                    vcf_var,
                    df,
                    df.query("loci == '@loci'"),
                    df.query("CHROM == '4'"),
                )
                raise
            if filter_value != "PASS":
                keeper = False
                break
        if keeper:
            passed_filtering.append(allele)

    bg.filtered_out["FILTER_not_PASS"] = [
        allele
        for allele in bg.alleles[AlleleState.RAW]
        if allele not in passed_filtering
    ]
    bg.alleles[AlleleState.FILT] = passed_filtering

    return bg


@apply_to_dict_values
def remove_alleles_with_low_base_quality(
    bg: BloodGroup,
    variant_metrics: dict[str, str],
    min_base_quality: int,
    microarray: bool,
) -> BloodGroup:
    """
    Remove alleles from a BloodGroup object that have defining variants with base
    quality below a specified minimum threshold.

    Args:
        bg (BloodGroup): The BloodGroup object containing alleles to filter.
        variant_metrics (dict[str, dict[str, int]]): A dictionary containing variant
        metrics with read depth information.
        min_base_quality (int): The minimum base_quality threshold.

    Returns:
        BloodGroup: The BloodGroup object with alleles filtered based on read depth.
    """

    filtered_out, passed_filtering = filter_vcf_metrics(
        bg.alleles[AlleleState.FILT],
        variant_metrics,
        "GQ",
        min_base_quality,
        microarray,
    )
    if filtered_out:
        vars_affected = ",".join(filtered_out.keys())
        message = f"Base Quality. Sample: {bg.sample}, BG: {bg.type}, variant/s: {vars_affected}"
        logger.warning(message)
    bg.filtered_out["insufficient_min_base_quality"] = filtered_out
    bg.alleles[AlleleState.FILT] = passed_filtering

    return bg


def get_fully_homozygous_alleles(
    ranked_chunks: list[list[Allele]], variant_pool: dict[str, Any]
) -> list[list[Allele]]:
    """Filter out alleles that are not fully homozygous from a list of ranked allele chunks.

    Uses a partial function to check each allele's variants in the provided `variant_pool`.
    Only alleles where every relevant variant equals the required homozygous genotype (2)
    are included in the result.

    Args:
        ranked_chunks (list[list[Allele]]):
            A list of lists (chunks), where each chunk contains ranked Allele objects.
        variant_pool (dict[str, Any]):
            A dictionary containing variant data used for assessing homozygosity.
            The exact structure depends on the `check_available_variants` function.

    Returns:
        list[list[Allele]]:
            A list of lists, each mirroring the structure of `ranked_chunks`
            but including only alleles that are fully homozygous in every variant.

    Raises:
        KeyError:
            If a variant key is missing in `variant_pool`.
    """
    check_hom = partial(check_available_variants, 2, variant_pool, operator.eq)
    homs = [[] for _ in ranked_chunks]

    for i, chunk in enumerate(ranked_chunks):
        for allele in chunk:
            if all(check_hom(allele)):
                homs[i].append(allele)
    return homs


def unique_in_order(lst: list) -> list:
    """
    Return a list of unique elements from 'lst' in the order they first appear,
    without using a set or other unordered data structure.

    Args:
        lst: The input list (possibly with duplicates).

    Returns:
        A list of items from 'lst' with duplicates removed in order.

    Example:
        >>> unique_in_order([3, 3, 1, 2, 1, 3])
        [3, 1, 2]
    """
    unique_items = []
    for item in lst:
        # Append item only if it's not already in the unique list
        if item not in unique_items:
            unique_items.append(item)
    return unique_items


# -----------------------------------------------------------
# Protocol for structural subtyping
# -----------------------------------------------------------
class GeneticProcessingProtocol(Protocol):
    """Protocol defining a process method for genetic data."""

    def process(
        self, bg: BloodGroup, reference_alleles: dict[str, Allele]
    ) -> list[Pair]:
        """Process a BloodGroup and return `AlleleState.NORMAL` pairs."""
        ...


# -----------------------------------------------------------
# Concrete strategies
# -----------------------------------------------------------
@dataclass
class NoVariantStrategy:
    """Handles the case where there are no variants."""

    def process(
        self, bg: BloodGroup, reference_alleles: dict[str, Allele]
    ) -> list[Pair]:
        ref_allele = reference_alleles[bg.type]
        return [Pair(ref_allele, ref_allele)]


@dataclass
class SingleVariantStrategy:
    """Handles the case where there is a single variant."""

    def process(
        self, bg: BloodGroup, reference_alleles: dict[str, Allele]
    ) -> list[Pair]:
        return [
            make_pair(
                reference_alleles, bg.variant_pool_numeric, bg.alleles[AlleleState.FILT]
            )
        ]


@dataclass
class MultipleVariantDispatcher:
    """Chooses a sub-strategy when multiple variants are present."""

    def process(
        self, bg: BloodGroup, reference_alleles: dict[str, Allele]
    ) -> list[Pair]:
        options = unique_in_order(bg.alleles[AlleleState.FILT])
        non_ref_options = get_non_refs(options)
        ranked_chunks = chunk_geno_list_by_rank(non_ref_options)
        homs = get_fully_homozygous_alleles(ranked_chunks, bg.variant_pool_numeric)

        first_chunk = ranked_chunks[0]
        weight_first_chunk = first_chunk[0].weight_geno
        trumpiest_homs = homs[0]
        weight_trumpiest_homs = (
            trumpiest_homs[0].weight_geno if trumpiest_homs else 1000
        )

        # Sub-strategy selection
        if len(trumpiest_homs) == 1:
            return SingleHomMultiVariantStrategy(
                hom_allele=trumpiest_homs[0], first_chunk=first_chunk
            ).process(bg, reference_alleles)
        elif len(trumpiest_homs) > 1 and weight_first_chunk == weight_trumpiest_homs:
            return MultipleHomMultiVariantStrategy(
                homs=trumpiest_homs, first_chunk=first_chunk
            ).process(bg, reference_alleles)
        elif any(len(hom_chunk) > 0 for hom_chunk in homs):
            return SomeHomMultiVariantStrategy(ranked_chunks=ranked_chunks).process(
                bg, reference_alleles
            )
        else:
            return NoHomMultiVariantStrategy(non_ref_options=non_ref_options).process(
                bg, reference_alleles
            )


@dataclass
class SingleHomMultiVariantStrategy:
    hom_allele: Allele
    first_chunk: list[Allele]

    def process(
        self, bg: BloodGroup, reference_alleles: dict[str, Allele]
    ) -> list[Pair]:
        hom_pair = [Pair(self.hom_allele, self.hom_allele)]
        if len(self.first_chunk) == 1:
            return hom_pair
        elif any(self.hom_allele in allele for allele in self.first_chunk):
            return combine_all(self.first_chunk, bg.variant_pool_numeric)
        else:
            return hom_pair + combine_all(self.first_chunk, bg.variant_pool_numeric)


@dataclass
class MultipleHomMultiVariantStrategy:
    homs: list[Allele]
    first_chunk: list[Allele]

    def process(
        self, bg: BloodGroup, reference_alleles: dict[str, Allele]
    ) -> list[Pair]:
        new_pairs = [Pair(h, h) for h in self.homs]
        if len(self.first_chunk) > len(self.homs):
            return new_pairs + combine_all(self.first_chunk, bg.variant_pool_numeric)
        return new_pairs


@dataclass
class SomeHomMultiVariantStrategy:
    ranked_chunks: list[list[Allele]]

    def process(
        self, bg: BloodGroup, reference_alleles: dict[str, Allele]
    ) -> list[Pair]:
        homs = get_fully_homozygous_alleles(self.ranked_chunks, bg.variant_pool_numeric)
        if len(homs) > 2 and len(homs[0]) == 0 and len(homs[1]) == 0:
            flat = [item for sublist in self.ranked_chunks for item in sublist]
            return combine_all(
                flat, bg.variant_pool_numeric
            )  # pass to filters (low_weight_hom)

        first_chunk = self.ranked_chunks[0]
        if len(first_chunk) == 1 and len(self.ranked_chunks) == 1:
            return [
                make_pair(
                    reference_alleles,
                    bg.variant_pool_numeric.copy(),
                    first_chunk,
                )
            ]
        return combine_all(
            self.ranked_chunks[0] + self.ranked_chunks[1], bg.variant_pool_numeric
        )


@dataclass
class NoHomMultiVariantStrategy:
    non_ref_options: list[Allele]

    def process(
        self, bg: BloodGroup, reference_alleles: dict[str, Allele]
    ) -> list[Pair]:
        ref_allele = reference_alleles[bg.type]
        ref_options = self.non_ref_options + [ref_allele]

        return combine_all(ref_options, bg.variant_pool_numeric)


# -----------------------------------------------------------
# Picking the right protocol-based strategy
# -----------------------------------------------------------
def _pick_strategy(bg: BloodGroup) -> GeneticProcessingProtocol:
    """Decide which strategy (protocol implementer) to use."""
    options = unique_in_order(bg.alleles[AlleleState.FILT])
    if len(options) == 0:
        return NoVariantStrategy()
    elif len(options) == 1:
        return SingleVariantStrategy()
    else:
        return MultipleVariantDispatcher()


@apply_to_dict_values
def process_genetic_data(
    bg: BloodGroup, reference_alleles: dict[str, Allele]
) -> BloodGroup:
    """Process genetic data to identify alleles and genotypes.

    Args:
        bg (BloodGroup):
            The blood group data that contains alleles (POS, NORMAL, etc.).
        reference_alleles (dict[str, Allele]):
            A dictionary mapping blood group types to reference Allele objects.

    Returns:
        An updated BloodGroup with `AlleleState.NORMAL` pairs set appropriately.

    Raises:
        ValueError: When constraints in the multiple-variant scenario are violated.
    """

    strategy: GeneticProcessingProtocol = _pick_strategy(
        bg
    )  # Returns a Protocol implementer
    normal_pairs = strategy.process(bg, reference_alleles)
    bg.alleles[AlleleState.NORMAL] = normal_pairs

    return bg


@apply_to_dict_values
def find_what_was_excluded_due_to_rank(
    bg: BloodGroup, reference_alleles: dict[str, Allele]
) -> BloodGroup:
    """Find all possible allele pairs based on genetic data.

    If the pairs are not present, list them in
    filtered_out["excluded_due_to_rank*"].

    Args:
        bg (BloodGroup): A BloodGroup object containing alleles, the variant pool,
            and other genetic data.
        reference_alleles (dict[str, Allele]): A dictionary mapping blood group types
            to their reference Allele.

    Returns:
        BloodGroup: The updated BloodGroup with pairs excluded due to rank added to
            the filtered_out collections.
    """

    options = set(bg.alleles[AlleleState.FILT])
    non_ref_options = get_non_refs(options)
    if non_ref_options:
        for pair in combine_all(non_ref_options, bg.variant_pool_numeric):
            if pair not in bg.alleles[AlleleState.NORMAL]:
                bg.filtered_out["excluded_due_to_rank"].append(pair)
        ref_options = non_ref_options + [
            reference_alleles[non_ref_options[0].blood_group]
        ]
        for pair in combine_all(ref_options, bg.variant_pool_numeric):
            if pair not in bg.alleles[AlleleState.NORMAL]:
                bg.filtered_out["excluded_due_to_rank_ref"].append(pair)
        ranked_chunks = chunk_geno_list_by_rank(non_ref_options)
        homs = get_fully_homozygous_alleles(ranked_chunks, bg.variant_pool_numeric)
        for ranked_homs in homs:
            for hom in ranked_homs:
                pair = Pair(allele1=hom, allele2=hom)
                if pair not in bg.alleles[AlleleState.NORMAL]:
                    bg.filtered_out["excluded_due_to_rank_hom"].append(pair)

    return bg


def make_pair(
    reference_alleles: dict[str, str], variant_pool: list[str], sub_results: list[str]
) -> list[str]:
    """Creates a pair of alleles based on the given parameters.

    Args:
        reference_alleles (Dict[str, str]): A mapping from blood group to reference
        allele.
        variant_pool (List[str]): A list of available variants.
        sub_results (List[str]): A list containing the initial results, expected to be
        of length 1.

    Returns:
        List[str]: A list containing the original results and an additional allele,
        either a duplicate of the first (if checks pass) or a corresponding
        reference allele.

    Raises:
        AssertionError: If the length of `sub_results` is not 1.
    """
    sub_results = list(sub_results)
    check_vars = partial(check_available_variants, 2, variant_pool, operator.eq)
    assert len(sub_results) == 1
    if all(check_vars(sub_results[0])):  # this is essentially fully_hom (func)
        sub_results.append(sub_results[0])
    else:
        sub_results.append(reference_alleles[sub_results[0].blood_group])
    return Pair(*sub_results)


def pair_can_exist(
    pair: tuple[Allele, Allele], variant_pool_copy: dict[str, int]
) -> bool:
    """Check if a pair of alleles can exist based on the variant pool.

    NB: This is a bit of a misnomer, as it only subtracts in more complex cases,
    like "009Kenya A4GALT": A4GALT*01/A4GALT*02 is not possible because if
    'A4GALT*02' then 22:43089849_T_C is on the other side so it has to be
    'A4GALT*01.02' and not reference.

    Args:
        pair (tuple[Allele, Allele]): A tuple containing two Allele objects.
        variant_pool_copy (dict[str, int]): A dictionary mapping variant identifiers
            to their available counts.

    Returns:
        bool: True if the pair can exist based on the variant pool, False otherwise.
    """
    allele1, allele2 = pair
    if allele1.reference or allele2.reference:
        return True
    for variant in allele1.defining_variants:
        variant_pool_copy[variant] -= 1
    return all(variant_pool_copy[variant] >= 1 for variant in allele2.defining_variants)


def combine_all(alleles: list[Allele], variant_pool: dict[str, int]) -> list[Pair]:
    """Combine all alleles into pairs, if possible.

    Args:
        alleles (list[Allele]): A list of Allele objects to be paired.
        variant_pool (dict[str, int]): A dictionary mapping variant identifiers
            to their available counts.

    Returns:
        list[Pair]: A list of Pair objects where each pair satisfies the variant
            pool constraints.
    """
    ranked = []
    for pair in itertools.combinations(alleles, 2):
        if pair_can_exist(pair, variant_pool.copy()):
            ranked.append(Pair(*pair))
    return ranked


@apply_to_dict_values
def add_CD_to_XG(bg: BloodGroup) -> BloodGroup:
    """TODO why not just use the CD99 vars??
    adds CD to XG blood group.

    Args:
        bg (BloodGroup): The BloodGroup object to be processed.

    Returns:
        BloodGroup: The processed BloodGroup object.
    """
    if bg.genotypes == ["XG*01/XG*01"]:
        bg.genotypes = ["XG*01/XG*01", "CD99*01/CD99*01"]
    return bg


def add_refs(db: Db, res: dict[str, BloodGroup], exclude) -> dict[str, BloodGroup]:
    """Add reference genotypes to existing results or create new entries for them.

    Args:
        db (Db): The database object containing reference alleles.
        res (Dict[str, BloodGroup]): Dictionary of BloodGroup objects to be updated
        with reference data.

    Returns:
        Dict[str, BloodGroup]: The updated dictionary of BloodGroup objects with added
          reference genotypes.

    This function checks for existing blood groups in the results dictionary and adds
    the reference genotype from the database if not present. It initializes a new
    BloodGroup object for any blood group type not already included in the results with
    the reference genotype as both a 'raw' and 'paired' allele.
    """
    for blood_group, reference in db.reference_alleles.items():
        if blood_group in exclude:
            continue
        if blood_group not in res:
            res[blood_group] = BloodGroup(
                type=blood_group,
                alleles={
                    AlleleState.RAW: [reference],
                    AlleleState.FILT: [reference],
                    AlleleState.NORMAL: [Pair(*[reference] * 2)],
                },
                sample="ref",
                genotypes=[f"{reference.genotype}/{reference.genotype}"],
            )
    return res


@apply_to_dict_values
def only_keep_alleles_if_best_big_del(
    bg: BloodGroup, matches: list[MatchResult]
) -> BloodGroup:
    """ """

    best = select_best_per_vcf(matches, tie_tol=1e-9)
    best_var_map = {}
    for match in best:
        best_var_map[f"{match.vcf.chrom}:{match.db.raw}"] = match.variant
    alleles_with_big_vars = []
    for allele in bg.alleles[AlleleState.RAW]:
        if any(var in best_var_map for var in allele.defining_variants):
            allele = allele.with_big_variants(best_var_map)
        alleles_with_big_vars.append(allele)
    bg.alleles[AlleleState.RAW] = alleles_with_big_vars

    return bg
