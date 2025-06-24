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
)
from rbceq2.db.db import Db
from rbceq2.IO.vcf import VCF
from icecream import ic

# def add_phase(instance: Allele, phases: list[str]) -> Allele:
#     """Add phase information to an Allele instance and return a new Allele with the
#     updated data.

#     Args:
#         instance (Allele): The original Allele instance to update.
#         phases (List[str]): A list of strings representing phase information.

#     Returns:
#         Allele: A new Allele instance with updated phase information.
#     """
#     new_instance_dict = instance.__dict__.copy()

#     new_instance_dict["phases"] = phases
#     del new_instance_dict["number_of_defining_variants"]
#     return Allele(**new_instance_dict)


def raw_results(db: Db, vcf: VCF, exclude: list[str]) -> dict[str, list[Allele]]:
    """Generate raw results from database alleles and VCF data based on phasing
    information.

    Args:
        db (Db): The database containing allele definitions and methods to generate
        them.
        vcf (VCF): The VCF data containing variants and possibly phasing information.

    Returns:
        Dict[str, List[Allele]]: A dictionary mapping blood groups to lists of Allele
        objects.
    """
    res: dict[str, list[Allele]] = defaultdict(list)

    for allele in db.make_alleles():
        if any(x in allele.genotype for x in exclude):
            continue
        if all(var in vcf.variants for var in allele.defining_variants):
            res[allele.blood_group].append(allele)

    return res


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

        This function finds the corresponding alternate allele at the same genomic
        position, retrieves its phase, and returns the inverse phase for the
        reference allele.

        Args:
            variant: The identifier of the reference variant (e.g., '4:144120554_ref').
            bg_variant_pool: A dictionary mapping variant IDs to their zygosity.
            phase_pool: A dictionary mapping variant IDs to their phase string (e.g., '1|0').

        Returns:
            The calculated phase string for the reference variant (e.g., '0|1').

        Raises:
            ValueError: If the input variant is not a valid heterozygous reference variant,
                        if a corresponding alternate variant cannot be found, or if the
                        phase format is invalid.
        """
        zygosity = bg.variant_pool.get(current_variant)
        if zygosity == Zygosity.HOM:
            return "."
        if "/" in phase_pool[current_variant]:
            return "unknown"
        # 1. Validate that the input variant is a heterozygous reference variant
        if not current_variant.endswith("_ref"):
            raise ValueError(
                f"Input variant '{current_variant}' must be a reference variant (ending in '_ref')."
            )

        if zygosity != Zygosity.HET:
            # ic(phase_pool)
            raise ValueError(
                f"Variant '{current_variant}' must be present and 'Heterozygous' in bg_variant_pool. "
                f"Found: {zygosity}"
            )

        # 2. Find the corresponding alternate allele variant at the same position
        position = current_variant.split("_")[0]
        partner_variant = None
        for key in phase_pool:
            # Match keys that start with the same position but are not the ref variant itself
            if key.startswith(position + "_") and key != current_variant:
                partner_variant = key
                break  # Found the partner, no need to continue searching
        # ic(current_variant, phase_pool)
        if partner_variant is None:
            return phase_pool[current_variant]

        # 3. Get the phase of the partner variant
        partner_phase = phase_pool.get(partner_variant)
        if partner_phase is None:
            raise ValueError(
                f"Partner variant '{partner_variant}' found but has no entry in phase_pool."
            )

        # 4. Calculate the ref phase by "flipping" the partner's phase
        phase_parts = partner_phase.split("|")
        if len(phase_parts) != 2:
            raise ValueError(
                f"Invalid phase format for partner '{partner_variant}': '{partner_phase}'"
            )

        # The ref allele is on the opposite haplotype from the alt allele
        # e.g., if partner is 1|0, ref is 0|1
        ref_phase = f"{phase_parts[1]}|{phase_parts[0]}"

        return ref_phase

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
            # try:
            chrom, pos_str = current_variant.split("_")[0].split(':')
            pos = int(pos_str)
            # except (ValueError, AttributeError):
            #     logger.warning(f"Invalid loci format for query: {current_variant}")
            #     return None

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
        # ic(current_variant, phase_set_pool)
        if partner_variant is None:
            return get_phase_set_for_loci()
        return phase_set_pool.get(partner_variant)

    if phased:  # TODO enum for GT, PS etc
        phase_pool = {
            variant: variant_metrics[variant]["GT"] for variant in bg.variant_pool
        }
        phase_set_pool = {
            variant: variant_metrics[variant].get("PS", "None")
            for variant in bg.variant_pool
        }
        phase_pool_ref_fixed = {}
        # if bg.type == 'ABO':
        #     ic(1111,phase_pool, phase_set_pool)
        for variant, phase in phase_pool.items():
            if variant.endswith("_ref"):
                new_phase = assign_ref_phase(variant)
                # ic(variant,new_phase, phase)
                phase_pool_ref_fixed[variant] = new_phase
            else:
                phase_pool_ref_fixed[variant] = phase
        phase_set_pool_ref_fixed = {}

        for variant, phase in phase_set_pool.items():
            #ic(variant)
            if variant.endswith("_ref"):
                new_phase = assign_ref_phase_set(variant)
                # ic(variant,new_phase, phase)
                phase_set_pool_ref_fixed[variant] = new_phase
            else:
                phase_set_pool_ref_fixed[variant] = phase
        bg.variant_pool_phase = phase_pool_ref_fixed
        bg.variant_pool_phase_set = phase_set_pool_ref_fixed
        # if bg.type == 'ABO':
        #     ic(22222,phase_pool_ref_fixed, phase_set_pool_ref_fixed)
        # raw_phased = []
        # for allele in bg.alleles[AlleleState.FILT]:
        #     phases = defaultdict(dict)
        #     # phase_sets = {}
        #     for variant in allele.defining_variants:
        #         if variant.endswith("_ref"):
        #             # if bg.type == "GYPA":
        #             #     ic(
        #             #         variant,
        #             #         bg.variant_pool,
        #             #         phase_pool,
        #             #         phase_set_pool,
        #             #         assign_ref_phase(variant, bg.variant_pool, phase_pool),
        #             #     )
        #             phase_of_variant = assign_ref_phase(variant, bg.variant_pool, phase_pool)
        #             phase_set = phase_set_pool[variant]
        #         else:
        #             phase_of_variant = variant_metrics[variant].get("GT", ".")
        #             phase_set = variant_metrics[variant].get("PS", ".")
        #         # TODO GQ etc like this too
        #         # GT:GQ:DP:AD:AF:PS
        #         # OR
        #         # GT:GQ:DP:AF:PS
        #         # or??? TODO

        #         # TODO is PS called IPS sometimes? ie GIAB
        #         phases[variant]["phase_set"] = phase_set
        #         phases[variant]["phase"] = phase_of_variant
        #         ic(phases)

        #     allele_phased = add_phase(allele, phases)
        #     raw_phased.append(allele_phased)
        #     block_numbers = [phase for phase in phases if phase != "."]
        #     if len(set(block_numbers)) > 1:
        #         logger.warning(
        #             f"first example of unphased vars!!! {allele.defining_variants} {phases}"
        #         )
        # assert len(bg.alleles[AlleleState.FILT]) == len(raw_phased)
        # bg.alleles[AlleleState.FILT] = raw_phased
        # # ic(11111,bg.alleles[AlleleState.FILT])

    return bg


# @apply_to_dict_values
# def add_phasing_for_antitheticals(
#     bg: BloodGroup, phased: bool, antitheticals: dict[str, str]
# ) -> BloodGroup:
#     """There are actually 5 SNPs that define JK*01N.20, the 4 shown plus it assumes
#     that 18:45739554 is ref (right??). But we know that '18:45739554_G_A' is on the
#     same strand. So, this is useful to rule out alleles. But it makes it very hard
#     to call alleles with phased data as there's no phase data for SNPs positions
#     that are meant to be reference...
#     chr18	JK*01			c.838G

#     chr18	JK*01N.20		c.226G>A;c.28G>A;c.303G>A;c.588A>G

#     Steps:

#     1. Check the ref pos
#         if not in VCF:
#             hom - phased
#         elif hom alt:
#             phased
#         else:
#             get phase of SNP - can rule out that, then just have to assume that it
#             is the other phase (if only 1) or unknown phase if multiple options

#     TODO - this is covered (almost???) by filter_pairs_by_phase
#     just by checking if phase of pair1 == pair2
#     should be enough, right? cuz the dif subtypes are based on antitheticals and they're
#     internally phased (even if they skip the ref (no change ie c.838G)).

#     hmm TODO manually validate a bunch (get Alexis to maybe??)

#     FYI its GM19473 that has:
#      all_phases: {'18:45730348_G_A': '20911299',
#                  '18:45731089_G_A': '20911299',
#                  '18:45731166_G_A': '20911299',
#                  '18:45736573_A_G': '.',
#                  '18:45739554_G_A': '20911299'} (c.838G)
#     """
#     if phased:
#         if bg.type == "JK":
#             pass

#     return bg


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
        keep = True
        for variant in allele.defining_variants:
            read_depth = float(variant_metrics[variant][metric_name])
            if microarray:
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
        bg.alleles[AlleleState.RAW], variant_metrics, "DP", min_read_depth, microarray
    )
    if filtered_out:
        vars_affected = ",".join(filtered_out.keys())
        message = f"Read Depth. Sample: {bg.sample}, BG: {bg.type}, variant/s: {vars_affected}"
        logger.warning(message)
    bg.filtered_out["insufficient_read_depth"] = filtered_out
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


@apply_to_dict_values
def rule_out_impossible_alleles(bg: BloodGroup) -> BloodGroup:
    """
    Rule out impossible alleles based on Homozygous variants.

    Args:
        bg (BloodGroup): A BloodGroup object containing allele states and phasing
            information.

    Returns:
        bg (BloodGroup): A BloodGroup object containing allele states and phasing
            information.

    Example:

    Allele(genotype='JK*02W.03',
            defining_variants={'18:43310415_G_A',
                                18:43316538_A_G',
                               '18:43319519_G_A'},
    Allele(genotype='JK*02W.04',
           defining_variants={'18:43310415_G_A',
                              '18:43319519_G_A'},

    '18:43310415_G_A': 'Heterozygous',
    '18:43316538_A_G': 'Homozygous',
    '18:43319519_G_A': 'Heterozygous',

    JK*02W.04 is impossible because '18:43316538_A_G' is Homozygous
    """
    # TODO mv this to filters.geno and apply after 'process_genetic_data', like the phased version
    # too scared to do it now in the middle of a validation round
    homs = {
        variant for variant, zygo in bg.variant_pool.items() if zygo == Zygosity.HOM
    }
    impossible_alleles = []
    impossible_pairs = []
    for allele in bg.alleles[AlleleState.FILT]:
        for allele2 in bg.alleles[AlleleState.FILT]:
            if allele.genotype != allele2.genotype and allele2 in allele:
                dif = allele.defining_variants.difference(allele2.defining_variants)
                if all(variant in homs for variant in dif):
                    impossible_alleles.append(allele2)
                    impossible_pairs.append(Pair(allele1=allele2, allele2=allele))

    bg.alleles[AlleleState.POS] = [
        allele
        for allele in bg.alleles[AlleleState.FILT]
        if allele not in impossible_alleles
    ]
    bg.filtered_out["allele1_not_possible_due_to_allele2"] = impossible_pairs

    return bg


# @apply_to_dict_values
# def rule_out_impossible_alleles_phased(bg: BloodGroup, phased: bool) -> BloodGroup:
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

#     to_test = 'YT'
#     if phased:
#         impossible_alleles = set([])
#         impossible_pairs = set([])
#         for allele in bg.alleles[AlleleState.POS]:
#             if bg.type == to_test:
#                 ic(11111, allele)
#             if allele in impossible_alleles:
#                 if bg.type == to_test:
#                     ic(22222)
#                 continue

#             if bg.type == 'KN' and allele.genotype not in ['KN*01','KN*02']:
#                 continue

#             if check_hom(allele):
#                 if bg.type == to_test:
#                     ic(3333333)
#                 continue
#             for allele2 in bg.alleles[AlleleState.POS]:
#                 if bg.type == to_test:
#                     ic(44444444, allele2)
#                 if allele == allele2:
#                     if bg.type == to_test:
#                         ic(555555)
#                     continue
#                 if bg.type == 'KN' and not allele2.genotype not in ['KN*01','KN*02']:
#                     continue
#                 if check_hom(allele2):
#                     if bg.type == to_test:
#                         ic(6666666)
#                     continue
#                 if allele in allele2:
#                     if bg.type == to_test:
#                         ic(7777777,allele.genotype, allele2.genotype)
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
#             ic(impossible_alleles)
#         bg.alleles[AlleleState.POS] = [
#             allele
#             for allele in bg.alleles[AlleleState.POS]
#             if allele not in impossible_alleles
#         ]
#         bg.filtered_out["allele1_not_possible_due_to_allele2_phased"] = list(impossible_pairs)
#     if bg.type == to_test:
#         ic(bg.alleles[AlleleState.POS])
#     return bg


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
                reference_alleles, bg.variant_pool_numeric, bg.alleles[AlleleState.POS]
            )
        ]


@dataclass
class MultipleVariantDispatcher:
    """Chooses a sub-strategy when multiple variants are present."""

    def process(
        self, bg: BloodGroup, reference_alleles: dict[str, Allele]
    ) -> list[Pair]:
        options = unique_in_order(bg.alleles[AlleleState.POS])
        non_ref_options = get_non_refs(options)
        ranked_chunks = chunk_geno_list_by_rank(non_ref_options)
        homs = get_fully_homozygous_alleles(ranked_chunks, bg.variant_pool_numeric)

        first_chunk = ranked_chunks[0]
        trumpiest_homs = homs[0]
        # if bg.type == 'YT':
        #     ic(options, non_ref_options, ranked_chunks, homs, first_chunk, trumpiest_homs)
        # Sub-strategy selection
        if len(trumpiest_homs) == 1:
            return SingleHomMultiVariantStrategy(
                hom_allele=trumpiest_homs[0], first_chunk=first_chunk
            ).process(bg, reference_alleles)
        elif len(trumpiest_homs) > 1:
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
            raise ValueError("No homs in the first two rank tiers.")

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
        ref_allele = reference_alleles[
            bg.type
        ]  # TODO - this is putting refs back in that have already been discarded as impossible
        ref_options = self.non_ref_options + [ref_allele]
        return combine_all(ref_options, bg.variant_pool_numeric)


# -----------------------------------------------------------
# Picking the right protocol-based strategy
# -----------------------------------------------------------
def _pick_strategy(bg: BloodGroup) -> GeneticProcessingProtocol:
    """Decide which strategy (protocol implementer) to use."""
    options = unique_in_order(bg.alleles[AlleleState.POS])
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
    # ic(strategy)
    normal_pairs = strategy.process(bg, reference_alleles)
    bg.alleles[AlleleState.NORMAL] = normal_pairs
    # if bg.type == 'YT':
    #     ic(bg.alleles[AlleleState.NORMAL], bg.alleles[AlleleState.POS])
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

    options = set(bg.alleles[AlleleState.POS])
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
    """
    adds CD to XG blood group.

    Args:
        bg (BloodGroup): The BloodGroup object to be processed.

    Returns:
        BloodGroup: The processed BloodGroup object.
    """
    if bg.genotypes == ["XG*01/XG*01"]:
        bg.genotypes = ["XG*01/XG*01", "CD99*01/CD99*01"]
    # ic(22222,bg.alleles[AlleleState.FILT])
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
