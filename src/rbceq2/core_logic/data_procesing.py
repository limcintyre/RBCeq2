import itertools
import operator
import re
from collections import defaultdict
from dataclasses import dataclass
from functools import partial
from typing import Any, Protocol
from loguru import logger

from rbceq2.core_logic.alleles import Allele, BloodGroup, Pair
from rbceq2.core_logic.constants import (
    AlleleState,
    CRITICAL_VARIANTS,
    HAPLOID_SECOND_SLOT,
    SYNTHESISED_HOM_REF_GT,
    UNNAMED_SECOND_SLOT,
)
from rbceq2.core_logic.filter_semantics import filter_excludes_allele
from rbceq2.core_logic.utils import (
    BeyondLogicError,
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
                allele = allele.with_big_variants(var_map)
            res[allele.blood_group].append(allele)
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
    ) -> str:
        """
        Assigns the correct phase to a reference variant at a heterozygous site.
        """
        zygosity = bg.variant_pool.get(current_variant)
        if zygosity == Zygosity.HOM:
            return "1/1"
        if zygosity == Zygosity.HEM:
            # One copy of the locus, so the token is on the one chromosome that locus
            # has - which is what a haploid GT says, and '1' is what a real one leaves
            # here. Without this the fall-through returns the raw GT, and for a
            # synthesised lane row that is the SYNTHESISED_HOM_REF_GT sentinel, same
            # leak the NO_COPIES/NO_DATA guard below exists to stop.
            return "1"
        if zygosity in (Zygosity.NO_COPIES, Zygosity.NO_DATA):
            # No chromosome to be phased, or nothing measured to phase. Without this the
            # fall-through returns the raw GT, which for a synthesised lane row is the
            # SYNTHESISED_HOM_REF_GT sentinel - a GT leaking into a phase field.
            return "unknown"
        elif "/" in phase_pool[current_variant]:
            return "unknown"
        return phase_pool[current_variant]

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

    return bg


# @apply_to_dict_values
# def ABO_phasing(
#     bg: BloodGroup,
#     phased: bool,
# ) -> BloodGroup:
#     # aboO_phases2 = set([]) I'm opting out of this path but if it ever gets revisited
#     # do it at DB level #TODO look for vars that have only ever been observed in cis with
#     # ABO*O, or never been observed with it and take phasing info from there. Best
#     # to actually phase indels though
#     """9:133257521(GRCh38) and 9:136132908 (GRCh37) _T_TC or _ref are pivotal for ABO
#     calls but aren't always assigned a phase group by phasing algos

#     some alleles are same/similar except for this locus

#     ie HG04054 for 1kg data

#     Allele
#     genotype: ABO*B.01
#     defining_variants:
#                 9:133255935_G_T 0|1
#                 9:133257521_T_TC 0/1
#                 9:133257486_T_C 1/1
#                 9:133255928_C_G 0|1
#                 9:133256074_G_A 0|1
#                 9:133256028_C_T 0|1
#                 9:133256205_G_C 0|1
#                 9:133255801_C_T 0|1
#     weight_geno: 1000
#     phenotype: . or B
#     reference: False

#     Allele
#     genotype: ABO*O.01.41
#     defining_variants:
#                 9:133257486_T_C 1/1
#                 9:133255928_C_G 0|1
#                 9:133256074_G_A 0|1
#                 9:133256028_C_T 0|1
#                 9:133256205_G_C 0|1
#                 9:133257521_ref unknown
#                 9:133255801_C_T 0|1

#     this function will infer the phase of this locus from the phase of ABO*O.01
#     specific variants
#     """
#     if not phased:
#         return bg

#     if bg.type != "ABO":
#         return bg
#     # 261delG
#     c261delGs = [
#         "9:133257521_ref",
#         "9:133257521_T_TC",
#         "9:136132908_ref",
#         "9:136132908_T_TC",
#     ]
#     if any(bg.variant_pool.get(c261delG) == Zygosity.HOM for c261delG in c261delGs):
#         return bg

#     aboO = set([])
#     other = set([])
#     for allele in bg.alleles[AlleleState.FILT]:
#         if allele.genotype.startswith("ABO*O.01"):
#             for variant in allele.defining_variants:
#                 aboO.add(variant)
#         else:
#             for variant in allele.defining_variants:
#                 other.add(variant)
#     aboO_phases = set([])
#     for variant in aboO.difference(other):
#         if not variant.startswith(("9:133257521", "9:136132908")):
#             phase = bg.variant_pool_phase[variant]
#             aboO_phases.add(phase)

#     if len(aboO_phases) == 0:
#         return bg  # can't rescue ABO
#     if len(aboO_phases) > 1:
#         return bg  # can't rescue ABO
#     abo_phase = aboO_phases.pop()
#     not_abo_phase = "1|0" if abo_phase == "0|1" else "0|1"
#     new_phases = {}

#     for variant, phase in bg.variant_pool_phase.items():
#         if variant in c261delGs:
#             if variant.endswith("_ref"):
#                 new_phases[variant] = abo_phase
#             else:
#                 new_phases[variant] = not_abo_phase
#         else:
#             new_phases[variant] = phase
#     bg.variant_pool_phase = new_phases

#     return bg


def chrom_copies_for_bg(
    bg: BloodGroup, vcf: VCF, single_copy_types: dict[str, str]
) -> int:
    """How many allele slots this blood group's result has for this sample.

    1 only where both halves agree: the blood group lies outside PAR on X/Y
    (Db.get_single_copy_types, from the database), and this sample's caller emitted
    haploid GTs on that chromosome outside PAR (VCF._infer_haploid_chroms, from the GTs).
    2 otherwise - every autosomal blood group, every female sample, and every VCF that
    gave no ploidy evidence. Five of the nine e2e datasets are in that last group and
    return 2 in every cell. The four short read ones are the only e2e coverage this
    branch has: 482 of 967 samples and 1,603 of 3,209 get a single slot, always in XK,
    GATA1 and ATP11C, which are the three blood groups outside PAR on X.

    Both halves are needed and neither is sufficient. The database half alone would make
    XK haploid for everyone including females; the sample half alone has no way to answer
    for a blood group whose only variant was dropped as hom ref, which is state table row
    B2 - a male hemizygous for reference XK, where the pool is empty and the answer is
    still one slot.

    Args:
        bg (BloodGroup): The blood group being built.
        vcf (VCF): The VCF this sample came from, carrying haploid_chroms.
        single_copy_types (dict[str, str]): Blood group type -> chromosome, from Db.

    Returns:
        int: 1 or 2.
    """
    chrom = single_copy_types.get(bg.type)
    if chrom is not None and chrom in vcf.haploid_chroms:
        return 1
    return 2


# How many dissenting positions locus_copies_for_bg names before it summarises the rest.
# Enough to go and look at every one of them, few enough to stay one log line.
MAX_DISSENTERS_LOGGED = 6


def locus_copies_for_bg(
    bg: BloodGroup, vcf: VCF, loci_by_type: dict[str, dict[str, frozenset[int]]]
) -> int | None:
    """How many copies of this gene the caller says are present, or None if it did not say.

    Some callers encode gene copy number as GT ploidy, writing a haploid GT across a gene
    the sample carries once. That is a statement about the *gene*, not about the sample's
    chromosomes, and the two have different answers: one copy of RHD is still two allele
    slots, one of them holding no gene at all. Reading it as chromosome ploidy would report
    'RHD*01/-' and lose a chromosome the sample has.

    The evidence required is agreement across the whole gene: every database locus the VCF
    actually reported for this blood group has to be haploid. One haploid GT among diploid
    neighbours is a caller quirk, not a copy number, and is left to be rejected by name in
    get_ref - the same treatment as before. Requiring agreement is what keeps this from
    firing on stray malformed rows, and it needs no flag and no per-gene list: a caller that
    encodes copy number this way does it consistently, and one that does not never trips it.

    Only a *disagreement* is worth a warning, and disagreement means the two readings
    are genuinely competing. A gene the caller wrote diploid with a couple of haploid
    rows in it is not competing with anything - it is a diploid gene and a caller having
    a bad day, the ordinary case, and it takes the same silent None as a gene nobody
    typed. The warning is kept for the other shape: the gene is mostly haploid, so the
    caller does look like it was reporting a copy number, and unanimity is the only thing
    standing in the way. That is the case where the conservative answer may be costing a
    real call, so it names the dissenting positions - they are the whole diagnosis, and
    there are few of them by construction.

    A tie says neither, and stays quiet. On a sparsely called input a tie is usually one
    locus against one, which is not evidence of a copy number.

    Note this is a decision about what to *say*, not about what to return. Unanimity
    still governs the result, unchanged, and the quiet case returns exactly what it
    returned when it was loud.

    Deliberately returns None rather than 2 for 'no evidence'. Nothing here can distinguish
    a gene at two copies from a gene nobody typed, and recording the second as the first
    would be the same class of mistake as reading './.' as wildtype.

    Note this asks only about the gene. Whether the *sample* has one chromosome there is
    chrom_copies_for_bg's question, and the two are read together at output: a gene at one
    copy on one chromosome is ordinary hemizygosity, on two chromosomes it is a missing
    copy that needs naming.

    Args:
        bg (BloodGroup): The blood group being built.
        vcf (VCF): The VCF this sample came from, carrying the per locus ploidy.
        loci_by_type (dict[str, dict[str, frozenset[int]]]): Blood group -> chromosome ->
        the positions that can vote on its copy number, from Db. Small variant positions
        only - see Db.get_loci_by_type for why a breakpoint must not be one of them.

    Returns:
        int | None: 1 where the caller reported the gene consistently single copy, None
        where it said nothing or contradicted itself.
    """
    by_chrom = loci_by_type.get(bg.type)
    if not by_chrom:
        return None

    haploid: list[tuple[str, int]] = []
    diploid: list[tuple[str, int]] = []
    for chrom, positions in by_chrom.items():
        for pos in positions & vcf.haploid_loci.get(chrom, frozenset()):
            haploid.append((chrom, pos))
        for pos in positions & vcf.diploid_loci.get(chrom, frozenset()):
            diploid.append((chrom, pos))

    if not haploid:
        return None
    if diploid:
        if len(haploid) > len(diploid):
            dissenters = sorted(diploid)[:MAX_DISSENTERS_LOGGED]
            shown = ", ".join(f"{chrom}:{pos}" for chrom, pos in dissenters)
            if len(diploid) > len(dissenters):
                shown += f" and {len(diploid) - len(dissenters)} more"
            logger.warning(
                f"{bg.sample}: {bg.type} has {len(haploid)} haploid and "
                f"{len(diploid)} diploid genotypes, so the copy number is "
                f"contradictory and is not being read as one. The gene reads as "
                f"single copy apart from {shown}. See issue #40"
            )
        return None
    return 1


def variant_was_discarded(vcf_var: str, df: pd.DataFrame) -> bool:
    """Did this token's own VCF row fail the FILTER classification?

    Not the same question as 'is this token absent from the variant pool', and conflating
    the two is what made a sample read as homozygous reference at a locus the caller called
    heterozygous with a PASS. A token leaves the pool whenever the only allele carrying it
    is removed, and that removal may have been over an entirely different variant - GYPA*01
    needs three positions, so one doubtful call at the third takes two good calls with it.
    Only a token the caller itself flagged is evidence that a heterozygous partner is gone.

    Returns False for a token with no row rather than raising. The caller is deciding
    whether to promote a '_ref' token to homozygous, and no row is no evidence for that -
    unlike only_keep_alleles_if_FILTER_PASS, where a missing row means the allele should
    never have been built and the raise is the point.

    Args:
        vcf_var (str): The variant token, ie the name used in the pool.
        df (pd.DataFrame): The sample's VCF rows.

    Returns:
        bool: True where the row exists and its FILTER value means the call is doubtful.
    """
    filter_value = filter_value_for(vcf_var, df)
    if filter_value is None:
        return False
    excludes, _ = filter_excludes_allele(filter_value)

    return excludes


def filter_values_for(vcf_var: str, df: pd.DataFrame) -> list[str]:
    """Every FILTER value carried by the rows this token matched, in file order.

    The one place the lookup itself lives. It was written out twice, and the two copies
    differed in what they did with no row at all, which is a difference that has to
    survive - see filter_value_for and only_keep_alleles_if_FILTER_PASS.

    A list rather than a value because the answer is genuinely plural: a file can carry
    more than one row for the same variant, and a row's variant cell can hold several
    tokens comma joined, which is why the match is by substring rather than by equality.
    Callers that want one value take the first, which is what both of them did before
    and what they still do - the plural form is here so that the choice is visible at
    the point it is made rather than hidden inside an index.

    The str() is a no op on every input form to hand - 141,729 FILTER cells across short
    read, long read and a jointly called cohort are all strings, lane synthesised rows
    included. It is there to make the annotation true rather than to convert anything. A
    cell with no value at all has never been seen and still has no defined behaviour: it
    failed before this and it fails after it, differently.

    Args:
        vcf_var (str): The variant token.
        df (pd.DataFrame): The sample's VCF rows.

    Returns:
        list[str]: The FILTER fields verbatim, empty where the token has no row.
    """
    matched = df.query("variant.str.contains(@vcf_var)")["FILTER"]

    return [str(value) for value in matched]


def filter_value_for(vcf_var: str, df: pd.DataFrame) -> str | None:
    """The FILTER field of the row this token came from, or None if it has no row.

    Args:
        vcf_var (str): The variant token.
        df (pd.DataFrame): The sample's VCF rows.

    Returns:
        str | None: The FILTER field verbatim, or None where the token has no row.
    """
    values = filter_values_for(vcf_var, df)

    return values[0] if values else None


def rows_disagree_about_exclusion(filter_values: list[str]) -> bool:
    """Would the verdict change depending on which of these rows was read?

    Not 'do the values differ'. Two rows can carry different values and still classify
    the same way, and on the one input form where this happens today they always do -
    'PASS' beside a value that marks which of two callers' rows is the one in conflict,
    both of which mean the call itself is not in doubt. Reporting those would be
    reporting nothing.

    What matters is a pair the classification splits, because then the answer depends on
    which row the file happens to list first. Nothing in the code holds that steady;
    what holds it steady is which side of filter_values.tsv each value sits on, so a
    value moving from one column to the other is enough to make the order load bearing
    with no code change anywhere.

    Args:
        filter_values (list[str]): The FILTER fields of every row one token matched.

    Returns:
        bool: True where the rows do not agree about whether to exclude.
    """
    if len(filter_values) < 2:
        return False

    return len({filter_excludes_allele(value)[0] for value in filter_values}) > 1


# How many order dependent tokens only_keep_alleles_if_FILTER_PASS names before it
# summarises the rest. One line per blood group either way.
MAX_ORDER_DEPENDENT_TOKENS_LOGGED = 4


def warn_if_the_row_order_decided_it(
    bg: BloodGroup, order_dependent: dict[str, list[str]]
) -> None:
    """Say out loud where the first row was taken and the second would have differed.

    The lookup reads the first row a token matched. Where every matching row classifies
    the same way that is a choice with no consequence, and this stays quiet. Where they
    split, the allele was kept or dropped on the strength of which row the file lists
    first, and nothing downstream can see that a second row said otherwise.

    Expected to be silent. It was measured at zero across every input form to hand, and
    that is the point rather than a reason to leave it out: the thing it watches is held
    steady by a classification table rather than by any code, so it goes from never to
    routinely the moment a value moves between columns of filter_values.tsv. Silence
    here is the evidence that the order still does not matter.

    One line per blood group, naming the tokens and both verdicts' values, because the
    fix is to look at the rows.

    Args:
        bg (BloodGroup): The blood group whose alleles were just filtered.
        order_dependent (dict[str, list[str]]): Token -> the FILTER values of every row
        it matched, for the tokens whose rows did not agree.
    """
    if not order_dependent:
        return
    tokens = sorted(order_dependent)
    shown = ", ".join(
        f"{token} ({', '.join(order_dependent[token])})"
        for token in tokens[:MAX_ORDER_DEPENDENT_TOKENS_LOGGED]
    )
    if len(tokens) > MAX_ORDER_DEPENDENT_TOKENS_LOGGED:
        shown += f" and {len(tokens) - MAX_ORDER_DEPENDENT_TOKENS_LOGGED} more"
    logger.warning(
        f"Sample: {bg.sample}, BG: {bg.type}. {len(tokens)} variant/s were reported on "
        f"more than one row, and the rows do not agree about whether the call is "
        f"doubtful. The first row in the file was used, so the answer here depends on "
        f"the order of the file: {shown}"
    )


def warn_if_critical_variant_not_trusted(bg: BloodGroup, df: pd.DataFrame) -> None:
    """Say out loud when a locus that decides a whole blood group was not trusted.

    Every exclusion is already in the log, but the log is per allele and this is a case
    where one row removes most of a blood group's definitions at once - the ABO c.261delG
    insertion is needed by 163 of them, and without it only the 43 that rest on its absence
    remain, so the sample reads as group O. Three samples in a densely called long read
    cohort were reported group O for exactly this reason, on rows the caller had scored
    QUAL 0 to 1.4 with the reference and alternate reads roughly even. A user reading the
    genotype has no way to see that from the output.

    Warned rather than raised, and the allele is still excluded - this does not change the
    call. The point is that the user knows to look, since 'A' and 'O' are equally plausible
    readings of the underlying data and only they can decide whether to rerun without the
    filter.

    One warning per locus per blood group, naming what was lost rather than how many rows
    matched, because the same locus can back several excluded alleles.

    Args:
        bg (BloodGroup): A BloodGroup whose FILTER exclusions have been recorded.
        df (pd.DataFrame): The sample's VCF rows.
    """
    for variant in sorted(
        {
            variant
            for allele in bg.filtered_out["FILTER_not_PASS"]
            for variant in allele.defining_variants
            if variant in CRITICAL_VARIANTS
        }
    ):
        filter_value = filter_value_for(variant, df)
        if filter_value is None:
            continue
        excludes, _ = filter_excludes_allele(filter_value)
        if not excludes:
            # The allele reached filtered_out over one of its other defining variants, so
            # this locus is not what was lost and blaming it would send the user to the
            # wrong row.
            continue
        dropped = sum(
            1
            for allele in bg.filtered_out["FILTER_not_PASS"]
            if variant in allele.defining_variants
        )
        excluded = (
            "1 allele needing it was excluded"
            if dropped == 1
            else f"{dropped} alleles needing it were excluded"
        )
        logger.warning(
            f"FILTER '{filter_value}' at {variant} - {CRITICAL_VARIANTS[variant]}. "
            f"Sample: {bg.sample}, BG: {bg.type}. {excluded}, so the call falls back to "
            f"whatever the rest support. This one row decides the answer for the whole "
            f"blood group, so check it before trusting this sample. Keeping it needs "
            f"--no_filter, which keeps every other flagged variant in the run too."
        )


@apply_to_dict_values
def make_variant_pool(
    bg: BloodGroup,
    vcf: VCF,
    single_copy_types: dict[str, str] | None = None,
    loci_by_type: dict[str, dict[str, frozenset[int]]] | None = None,
) -> BloodGroup:
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

    def find_matching_keys(dict_keys, ref_key):
        """to protect against instances like below
        ie when filter not passed, ref needs to be homozygous
        #Results:
        Genotypes count: 1
        Genotypes:
        Phenotypes (numeric):
        Phenotypes (alphanumeric):

        #Data:
        Vars:
        1:159205564_ref : Heterozygous
        Vars_phase:
        1:159205564_ref : 0|1
        Vars_phase_set:
        1:159205564_ref : 159155563
        Raw:
        Allele
        genotype: FY*01
        defining_variants:
                1:159205564_ref
        weight_geno: 1000
        phenotype: FY:1,-2 or Fy(a+),Fy(b-)
        reference: True

        Allele
        genotype: FY*02
        defining_variants:
                1:159205564_G_A
        weight_geno: 1000
        phenotype: FY:-1,2 or Fy(a-),Fy(b+)
        reference: False


        2025-11-12 09:39:51.322 | DEBUG    | ### Filters applied ###:

        2025-11-12 09:39:51.322 | DEBUG    |
        FILTER_not_PASS: Allele
        genotype: FY*02
        defining_variants:
                1:159205564_G_A
        weight_geno: 1000
        phenotype: FY:-1,2 or Fy(a-),Fy(b+)
        reference: False"""
        # Parse the reference key
        chrom, pos_ref = ref_key.split(":")
        pos = pos_ref.split("_")[0]

        # Find matching keys
        matches = [
            key
            for key in dict_keys
            if key.startswith(f"{chrom}:{pos}_") and key != ref_key
        ]

        return matches

    bg.chrom_copies = chrom_copies_for_bg(bg, vcf, single_copy_types or {})
    bg.locus_copies = locus_copies_for_bg(bg, vcf, loci_by_type or {})

    variant_pool = {}

    for allele in bg.alleles[AlleleState.FILT]:
        zygosity = {
            var: get_ref(vcf.variants[var], var, bg.chrom_copies, bg.locus_copies)
            for var in allele.defining_variants
        }
        variant_pool = variant_pool | zygosity

    for variant, zygo in variant_pool.items():
        if variant.endswith("_ref") and zygo == Zygosity.HET:
            matching = find_matching_keys(variant_pool.keys(), variant)
            if not matching:
                failed = []
                for allele in bg.filtered_out["FILTER_not_PASS"]:
                    for dropped in allele.defining_variants:
                        if dropped not in variant_pool and variant_was_discarded(
                            dropped, vcf.df
                        ):
                            failed.append(dropped)
                matching2 = find_matching_keys(failed, variant)
                if matching2:  # het pair gone
                    variant_pool[variant] = Zygosity.HOM
    bg.variant_pool = variant_pool
    check_token_copies_fit_chrom_copies(bg)

    return bg


def check_token_copies_fit_chrom_copies(bg: BloodGroup) -> None:
    """No token may claim more copies than the sample has chromosomes.

    The per-token form of the weak invariant in ploidy_state_table.md section 4. Vacuous
    while every region is two copies - Zygosity tops out at HOM, which is 2 - so this only
    ever fires on a single-copy region, where it catches one specific thing: a caller that
    mixes ploidy codings within one sample, writing '1' at one non-PAR X locus and '1/1'
    at another. Then the file says both one chromosome and two and there is no way to tell
    which is meant.

    Raised rather than warned. A wrong answer here is a wrong blood group, and unlike the
    contradictions in _modify_variant_pool_with_large_indel there is no defensible
    fallback - picking either reading invents a chromosome count. A BeyondLogicError
    rather than a bare assert so it survives 'python -O'.

    Args:
        bg (BloodGroup): A BloodGroup whose variant_pool has just been built.

    Raises:
        BeyondLogicError: If any token claims more copies than bg.chrom_copies.
    """
    for variant, copies in bg.variant_pool_numeric.items():
        if copies > bg.chrom_copies:
            raise BeyondLogicError(
                message=(
                    "A variant claims more copies than the sample has chromosomes "
                    "there. This happens when one sample mixes ploidy codings, ie a "
                    "haploid '1' at one non-PAR X locus and a diploid '1/1' at another. "
                    "The VCF is saying both one chromosome and two - see issue #40."
                ),
                context=(
                    f"sample: {bg.sample}, BG: {bg.type}, variant: {variant}, "
                    f"zygosity: {bg.variant_pool[variant]} ({copies} copies), "
                    f"chrom_copies: {bg.chrom_copies}"
                ),
            )


@apply_to_dict_values
def record_unused_variants(
    bg: BloodGroup,
    vcf: VCF,
    loci_by_type: dict[str, dict[str, frozenset[int]]],
    df: pd.DataFrame,
) -> BloodGroup:
    """Collect the variants at this blood group's loci that the pool does not hold.

    Diagnostic only. Nothing downstream reads what this writes; it exists so the debug
    trace can answer a question it currently cannot. When every allele a blood group had
    is discarded the variant pool comes out empty, and the trace then shows no variants
    at all - indistinguishable from a sample that was wildtype at every locus. Reading
    HG01174 GYPA meant opening the VCF by hand to discover the sample was 1/1 at all
    three loci with two of the three calls passing.

    The loci come from the database rather than from the pool, for the same reason
    locus_copies_for_bg uses them: by the time a pool exists it only holds what survived,
    so it cannot testify about what did not.

    Zygosity is derived the same way the pool's is, but failure to derive it is not
    allowed to end the run - this is a debug aid, and a genotype get_ref refuses is
    exactly the kind of thing worth seeing in the trace rather than dying on. Such a
    variant is recorded with its raw genotype instead.

    Phase and phase set are the caller's own, unprocessed: the reference-token fixups
    add_phasing applies are about deciding which chromosome an allele sits on, and
    nothing here sits on a chromosome.

    Args:
        bg (BloodGroup): The blood group, after make_variant_pool.
        vcf (VCF): The sample's VCF, for its variants and their metrics.
        loci_by_type (dict): Database positions per blood group, per chromosome.
        df (pd.DataFrame): The sample's rows, for each variant's FILTER.

    Returns:
        BloodGroup: With unused_pool and its three companions filled in.
    """
    positions = loci_by_type.get(bg.type)
    if not positions:
        return bg

    # Built once rather than queried per variant. The df.query idiom used elsewhere is
    # the slow path, and a row can carry several comma joined tokens from one
    # multi-allelic site, so each is mapped to its row's FILTER separately.
    filters: dict[str, str] = {}
    if df is not None and not df.empty:
        for raw, filter_value in zip(df["variant"], df["FILTER"]):
            for token in str(raw).split(","):
                filters.setdefault(token.strip(), filter_value)

    for variant, metrics in vcf.variants.items():
        if variant in bg.variant_pool:
            continue
        locus, _, _ = variant.partition("_")
        chrom, _, position = locus.partition(":")
        if not position.isdigit():
            continue
        if int(position) not in positions.get(chrom, frozenset()):
            continue
        try:
            bg.unused_pool[variant] = get_ref(
                metrics, variant, bg.chrom_copies, bg.locus_copies
            )
        except Exception:
            bg.unused_pool[variant] = metrics.get("GT", "?")
        bg.unused_pool_phase[variant] = metrics.get("GT", ".")
        bg.unused_pool_phase_set[variant] = metrics.get("PS") or "."
        # A '_ref' token asserts the reference and has no FILTER of its own, which is
        # why only_keep_alleles_if_FILTER_PASS skips it. Saying so beats printing the
        # row's value: a lane locus with no called rows is synthesised from COMMON_COLS,
        # so its FILTER field holds the literal string 'FILTER'.
        if "_ref" in variant:
            bg.unused_pool_filters[variant] = "not read - reference token"
        else:
            bg.unused_pool_filters[variant] = filters.get(variant, ".")

    return bg


@apply_to_dict_values
def remove_alleles_with_no_call_variants(bg: BloodGroup) -> BloodGroup:
    """Remove alleles that depend on a locus the caller did not call.

    An allele is only reported if every one of its defining variants was measured. If any
    of them is Zygosity.NO_DATA the allele cannot be confirmed, so it is excluded here
    with a named reason rather than being left in play - every one of the 18 sites that
    compares a zygosity does so with '== Zygosity.HOM' or '!= Zygosity.HET', all of which
    evaluate the wrong way for NO_DATA and would quietly stop a filter firing with nothing
    recorded in filtered_out.

    The NO_DATA token is deliberately left in bg.variant_pool. It is the evidence for the
    exclusion, and the Vars: block of the debug log is where a user goes to see why an
    allele went away. Dropping the token instead would remove the exclusion's own
    explanation from the trail.

    Note this is not the same check as only_keep_alleles_if_FILTER_PASS, and neither
    subsumes the other: FILTER is a per-site column and can read PASS on a row where this
    sample's GT is './.'. Measured on the 1kg microarray set, all 94 rows at DB loci are
    FILTER PASS while 5,343 sample genotypes at those loci are no-calls.

    Args:
        bg (BloodGroup): A BloodGroup whose variant_pool has been built.

    Returns:
        BloodGroup: The BloodGroup with unconfirmable alleles moved to filtered_out under
            'no_call_at_defining_variant'.
    """
    to_remove = [
        allele
        for allele in bg.alleles[AlleleState.FILT]
        if any(
            bg.variant_pool.get(variant) == Zygosity.NO_DATA
            for variant in allele.defining_variants
        )
    ]
    bg.remove_alleles(to_remove, "no_call_at_defining_variant", AlleleState.FILT)

    return bg


def _modify_variant_pool_with_large_indel(
    variant_pool: dict,
    sample: str,
    bg_type: str,
    is_phase_pool: bool = False,
    chrom_copies: int = 2,
) -> dict:
    """Internal helper to adjust variant zygosity when large deletions are present.

    A deletion reaches here labelled HOM, HET or HEM, and HEM is the one that means two
    opposite things. It is one copy of the locus carrying one copy of the token, which
    says nothing on its own about how many copies there were to start with, so the region
    has to supply that:

        chrom_copies == 1   the only chromosome carries the deletion, so everything
                            inside it is gone - the 'hom' reading. An XK whole gene
                            deletion in a male (McLeod, and the reason XK*N.01 exists).
        chrom_copies == 2   one chromosome of two carries it, so a variant inside can sit
                            on the other - the 'het' reading. This is what a copy number
                            aware caller writes for a heterozygous deletion, a haploid
                            '1' rather than '0/1'.

    Both are named explicitly. Until v2.4.6 only the first was, and the second tripped a
    bare assert - so an ordinary heterozygous deletion from a copy number aware caller
    crashed the sample, or, under 'python -O', silently did not. See issue #40.

    What remains genuinely beyond logic is a deletion that is *neither*: it removed some
    copies but not all, so it has to be on one chromosome of two, and if it is not het by
    either reading then the number of surviving copies of the variants inside it cannot be
    worked out. That raises, and it raises by name.

    Args:
        variant_pool (dict): Dictionary mapping variant strings to zygosity values.
        sample (str): Sample identifier for logging.
        bg_type (str): Blood group type for logging.
        is_phase_pool (bool): If True, expects phase notation ('1/1', '1'),
            otherwise Zygosity enums. Defaults to False.
        chrom_copies (int): Copies of the region the sample was born with. Defaults to 2.

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
        # A tuple in both branches, so the membership test below means membership. It was
        # a bare string in the zygosity branch, where 'in' is substring matching, and it
        # gave the right answer only because 'Heterozygous' contains itself.
        het_values = ("1|0", "0|1")
        hem_value = "1"
        no_data_values = ("./.", ".|.", ".")
        # The phase pool has no equivalent of 'no copies'. Phase is a property of a
        # chromosome that exists, so there is nothing honest to write inside a hom
        # deletion. None means 'leave the phase alone' - the zygosity pool is what records
        # the absence, and by the time any phased filter runs, the alleles that would read
        # this phase have been excluded.
        no_copies_value = None
    else:
        hom_value = Zygosity.HOM
        het_values = (Zygosity.HET,)
        hem_value = Zygosity.HEM
        no_data_values = (Zygosity.NO_DATA,)
        no_copies_value = Zygosity.NO_COPIES

    new_variant_pool = {}

    for big_del, no_seq_variant in big_dels:
        # An uncalled deletion is not evidence of anything, so nothing inside it is
        # adjusted. This is what a jointly called cohort produces and a per sample file
        # does not: the cohort carries a row for every structural variant *any* sample
        # had, so a sample without this one gets './.' rather than no row at all. Both
        # encodings describe the same sample and have to reach the same answer, and
        # skipping is what makes them agree - it leaves the pool exactly as the per
        # sample file would have left it by simply not having the row.
        #
        # Deliberately not read as 'the deletion is absent' either. './.' is the caller
        # declining to call, which is a different claim from a confident reference call,
        # and the difference is the whole reason NO_DATA exists.
        if variant_pool.get(big_del) in no_data_values:
            continue
        start = get_start_pos(no_seq_variant)
        length = no_seq_variant.split("_")[-1]
        length = int(length[:-2]) * 1000 if length.endswith("kb") else int(length)
        end = start + length

        # Check if deletion is homozygous - ie removes every copy of what it spans. On a
        # single-copy region hemizygous means exactly that, so it counts here too.
        del_zygosity = variant_pool.get(big_del)
        big_del_is_hom = del_zygosity == hom_value or (
            chrom_copies == 1 and del_zygosity == hem_value
        )
        # The mirror of the line above, and the other half of what a haploid GT on a
        # deletion means. HEM is one copy of the locus carrying one copy of the token, so
        # what it says about a *deletion* depends entirely on how many copies there were
        # to start with: on one chromosome the only copy is gone, which is the hom
        # reading above; on two, one of two is gone, which is what HET says. The label is
        # the same and the two readings are opposite, so both need naming.
        big_del_is_het = del_zygosity in het_values or (
            chrom_copies == 2 and del_zygosity == hem_value
        )

        for variant, zygosity in variant_pool.items():
            if variant == big_del:
                # Keep deletion as-is
                new_variant_pool[variant] = zygosity
            elif start < get_start_pos(variant) < end:
                # Variant falls within deletion range
                if big_del_is_hom:
                    if not variant.endswith("_ref"):
                        # An ALT call inside a deletion that removed both copies is a
                        # contradiction - one of the two calls is wrong. Warned rather
                        # than raised because it is an input problem, not a logic error
                        # here, and because the SV calls this is judged against are thin
                        # (see the fuzzy matching note in CLAUDE.md). Either way the
                        # locus is untrustworthy, so it is marked the same way.
                        logger.warning(
                            f"Non-reference variant inside a homozygous deletion, so "
                            f"both the deletion and the variant cannot be right "
                            f"{sample} {bg_type} {variant}"
                        )
                    # Was: kept as-is, ie Homozygous - a claim of two wildtype copies at a
                    # locus with no chromosomes at all. Nothing can be carried here.
                    new_variant_pool[variant] = (
                        zygosity if no_copies_value is None else no_copies_value
                    )
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
                    if not is_phase_pool and not big_del_is_het:
                        raise BeyondLogicError(
                            message=(
                                "A variant inside a deletion that is neither homozygous "
                                "nor heterozygous. Reaching here means the deletion "
                                "removed some copies but not all of them, so it has to "
                                "be on one chromosome of two, and this one is neither - "
                                "so how many copies of the variant survive cannot be "
                                "worked out. An uncalled deletion is skipped before this "
                                "point, so the remaining way to get here is a deletion "
                                "sitting inside another one."
                            ),
                            context=(
                                f"sample: {sample}, blood group: {bg_type}, deletion: "
                                f"{big_del}, deletion zygosity: {del_zygosity}, "
                                f"variant: {variant}, variant zygosity: {zygosity}, "
                                f"chrom_copies: {chrom_copies}"
                            ),
                            raised_by=(
                                "_modify_variant_pool_with_large_indel/"
                                "deletion_neither_hom_nor_het"
                            ),
                        )
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
        bg.variant_pool,
        bg.sample,
        bg.type,
        is_phase_pool=False,
        chrom_copies=bg.chrom_copies,
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
        bg.variant_pool_phase,
        bg.sample,
        bg.type,
        is_phase_pool=True,
        chrom_copies=bg.chrom_copies,
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

    The example above is what this function was written for, but until the C4 fix it could
    not happen: _modify_variant_pool_with_large_indel kept the '_ref' token as Homozygous,
    so every defining variant was still in the pool and nothing was ever removed. It is
    now marked Zygosity.NO_COPIES instead, which is what makes this reachable.

    Two separate reasons are recorded, because they are not the same problem:

    - hom_deletion_at_defining_variant - the allele needs a locus that has no chromosomes
      under it. Expected, and the point of this function.
    - defining_variant_missing_from_pool - the allele needs a variant that is not in the
      pool at all. Not expected; previously an `ic` to stdout and a bare assert that
      vanishes under `python -O`, and the alleles were dropped with nothing recorded.

    Args:
        bg (BloodGroup): A BloodGroup whose variant_pool has had large indels applied.

    Returns:
        BloodGroup: The BloodGroup with impossible alleles moved to filtered_out.
    """
    in_hom_del = []
    missing_from_pool = []

    for allele in bg.alleles[AlleleState.FILT]:
        if any(
            bg.variant_pool.get(variant) == Zygosity.NO_COPIES
            for variant in allele.defining_variants
        ):
            in_hom_del.append(allele)
        elif not all(
            variant in bg.variant_pool for variant in allele.defining_variants
        ):
            missing = sorted(
                variant
                for variant in allele.defining_variants
                if variant not in bg.variant_pool
            )
            logger.warning(
                f"Allele needs a variant that is not in the variant pool, which should "
                f"not be reachable - please report: {bg.sample} {bg.type} "
                f"{allele.genotype} {missing}"
            )
            missing_from_pool.append(allele)

    bg.remove_alleles(
        in_hom_del, "hom_deletion_at_defining_variant", AlleleState.FILT
    )
    bg.remove_alleles(
        missing_from_pool, "defining_variant_missing_from_pool", AlleleState.FILT
    )

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


def parse_GT(GT: str) -> tuple[str, ...]:
    """Split a genotype string into its allele indices.

    Splits on either separator, so phased and unphased genotypes parse the same way.
    No interpretation is applied - '.' is returned as '.', and multi-digit indices are
    returned whole.

    Args:
        GT (str): A genotype string, e.g. '0/1', '0|1', '1', './.'.

    Returns:
        tuple[str, ...]: The allele indices, e.g. ('0', '1') for '0/1', ('1',) for '1'.
    """
    return tuple(re.split(r"[/|]", GT))


def get_ref(
    ref_dict: dict[str, str],
    variant: str = "",
    chrom_copies: int = 2,
    locus_copies: int | None = None,
) -> str:
    """Determine the zygosity from a reference dictionary containing genotype
    information.

    The genotype string is expected to be in the format '0/1', '0|1', etc., where
    the delimiter can be '/' or '|'.
    A genotype of '0/0' or '1/1', etc., where both alleles are the same, will return
    'Homozygous'.
    A genotype of '0/1', '1/0', etc., will return 'Heterozygous'.

    A GT containing '.' returns 'No_data'. Until v2.4.3 it was passed through
    .replace(".", "0"), which turned './.' into '0/0' and therefore into 'Homozygous' -
    converting the absence of a measurement into a positive claim of wildtype. That is
    wrong in both directions: for an ALT token it asserted the sample carries the variant
    on both chromosomes, and it is exactly the case the caller declined to call. In a
    microarray call set '0/0' and './.' both appear in the same file with different
    meanings (a confident hom ref call, and a failed probe), so './.' cannot be read as
    hom ref. Half-calls ('0/.', '1/.') are No_data too - one known allele is not a
    confirmed genotype.

    The synthesised lane row is the one legitimate assertion of wildtype and carries
    SYNTHESISED_HOM_REF_GT rather than a real GT, so it is recognised here and is not
    affected.

    A haploid GT is interpreted where either count is 1, and they are different claims
    reaching the same zygosity. chrom_copies is 1 outside PAR on X/Y, where the sample
    carries one chromosome. locus_copies is 1 where a caller reported one copy of the gene
    consistently across it, which since v2.4.5 covers autosomes - state table rows D2 and
    D3. Either way '1' is Hemizygous and '.' is No_data; '0' never arrives, because
    remove_home_ref drops it the way it drops '0/0'.

    Hemizygous is the honest label - one copy of the locus carrying one copy of the token -
    but it is not what the pairing machinery needs, so the *numeric* value comes from
    BloodGroup.len_dict via chrom_copies rather than from the enum. See make_variant_pool,
    which is where the normalisation is applied and explained.

    What the two cases differ in is the *shape of the result*, not the zygosity, and that
    difference is settled at output rather than here. One chromosome gives one allele slot;
    one copy of a gene on two chromosomes gives two, the second holding no gene at all. See
    get_genotypes.

    A haploid GT with neither count at 1 is still rejected. That is a locus whose
    neighbours in the same gene were diploid, so the file is claiming one copy and two at
    once and there is nothing to prefer between them.

    Multi-allelic genotypes are rejected rather than guessed at - the VCF needs splitting
    with 'bcftools norm -m -both' first. The '.' check runs before the biallelic check, so
    a nonsense GT like '2/.' reports as No_data rather than as multi-allelic.

    Args:
        ref_dict (Dict[str, str]): A dictionary containing the genotype ("GT") and
        possibly other information.
        variant (str): The variant the genotype belongs to, used to make errors
        traceable back to a VCF row. Optional so existing callers keep working.
        chrom_copies (int): Copies of the region the sample was born with. 2 unless the
        VCF layer found haploid GTs outside PAR on this chromosome. Optional so existing
        callers keep working.
        locus_copies (int | None): Copies of the gene the caller reported, from
        locus_copies_for_bg. None means it did not say, which is the common case.
        Optional so existing callers keep working.

    Returns:
        str: A string indicating the zygosity as 'Homozygous', 'Heterozygous',
        'Hemizygous' or 'No_data'.

    Raises:
        BeyondLogicError: If the genotype is haploid where neither the region nor the gene
        has one copy, or is multi-allelic.
    """
    # 0/1:41,47:88:99:1080,0,1068:0.534:99
    GT = ref_dict["GT"]

    if GT == SYNTHESISED_HOM_REF_GT:
        return Zygosity.HOM

    alleles = parse_GT(GT)

    if len(alleles) == 1 and 1 in (chrom_copies, locus_copies):
        if alleles[0] == ".":
            return Zygosity.NO_DATA
        if alleles[0] == "0":
            # remove_home_ref drops these, so reaching here means the pool was built from
            # a df that never went through it. Absence is still the right encoding.
            return Zygosity.NO_DATA
        if alleles[0] != "1":
            raise BeyondLogicError(
                message=(
                    "Multi-allelic genotypes are not supported. A haploid genotype "
                    "names an alternate allele this row does not have, so the row was "
                    "not split into one alternate per row. Please split the VCF first, "
                    "ie 'bcftools norm -m -both'."
                ),
                context=(
                    f"variant: {variant}, GT: {GT}, chrom_copies: {chrom_copies}, "
                    f"locus_copies: {locus_copies}"
                ),
                raised_by="get_ref/multi_allelic_haploid_GT",
            )
        return Zygosity.HEM

    if len(alleles) != 2:
        raise BeyondLogicError(
            message=(
                "A haploid genotype needs either one chromosome or one copy of the gene, "
                "and this locus has neither. Outside the pseudoautosomal regions of X/Y "
                "the sample carries one chromosome; where a caller encodes gene copy "
                "number as GT ploidy it writes haploid genotypes consistently across the "
                "whole gene. A single haploid genotype among diploid ones in the same "
                "gene is claiming one copy and two at once, so it is rejected rather "
                "than guessed at - see issue #40."
            ),
            context=(
                f"variant: {variant}, GT: {GT}, chrom_copies: {chrom_copies}, "
                f"locus_copies: {locus_copies}"
            ),
            raised_by="get_ref/haploid_GT_where_neither_count_is_one",
        )

    if "." in alleles:
        return Zygosity.NO_DATA

    if not set(alleles) <= {"0", "1"}:
        raise BeyondLogicError(
            message=(
                "Multi-allelic genotypes are not supported. A diploid genotype names an "
                "alternate allele this row does not have, so the row was not split into "
                "one alternate per row. Please split the VCF first, ie "
                "'bcftools norm -m -both'."
            ),
            context=f"variant: {variant}, GT: {GT}",
            raised_by="get_ref/multi_allelic_diploid_GT",
        )

    allele1, allele2 = alleles
    if allele1 == allele2:
        return Zygosity.HOM
    return Zygosity.HET


@apply_to_dict_values
def get_genotypes(
    bg: BloodGroup,
    reference_alleles: dict[str, Any] | None = None,
    gene_absent_subtypes: dict[str, str] | None = None,
) -> BloodGroup:
    """Generate genotype combinations for a given blood group from allele pairs.

    Args:
        bg (BloodGroup): The blood group object containing alleles.
        reference_alleles (dict[str, Any] | None): Blood group -> reference Allele, from
        Db. Needed only to recognise the reference slot when a gene copy is missing.
        Optional so existing callers keep working.
        gene_absent_subtypes (dict[str, str] | None): Blood group -> the subtype naming a
        missing gene copy, from Db. Optional so existing callers keep working.

    Returns:
        BloodGroup: The blood group object with updated genotypes based on allele
        combinations.

    This function processes 'pairs' and 'co_existing' alleles to create sorted genotype
    strings.

    Two different single copy states are rendered here, and keeping them apart is the whole
    point of there being two numbers. Both are *reporting* decisions: the pair itself is
    left alone in each case, so phenotype, filters and the exclusion trail all still see an
    ordinary two-allele pair.

    Where bg.chrom_copies is 1 the two slots hold the same allele - one chromosome carrying
    it, which the pairing machinery represents as a duplicate - and the second slot is
    written as HAPLOID_SECOND_SLOT. 'XK*N.03/-' says there is no second chromosome;
    'XK*N.03/XK*N.03' would be indistinguishable in the TSV from a female homozygote, and a
    bare 'XK*N.03' would break every consumer that splits a genotype on '/'.

    A third second-slot value is written elsewhere and only passed through here. Where
    cant_name_second_slot_cuz_ref_impossible identified one chromosome and refused the
    other, it leaves the rendered string in bg.single_slot_genotypes and removes the
    pair, so there is nothing left to render and the string is used as it stands.
    Unlike the two below it is not a copy number statement: both chromosomes are there
    and one of them
    carries an allele the database cannot name.

    Where bg.locus_copies is 1 but chrom_copies is 2 the sample has two chromosomes and one
    of them carries no copy of the gene. The pairing machinery has nothing to put there -
    an array reports copy number without a deletion record, so no deletion allele was ever
    built - and pairs the real allele with the reference instead, which asserts wildtype on
    a chromosome there is evidence against. So the reference slot is replaced: by the
    database's subtype for a missing copy where it has one, and by UNNAMED_SECOND_SLOT
    where it does not. Subtype rather than allele because a copy number carries no
    breakpoints and cannot say which deletion it was.

    Example:
        XK, male, X:37686068 G>A called '1':
            pair       -> Pair(XK*N.16, XK*N.16)
            genotypes  -> ['XK*N.16/-']

        RHD, one gene copy, 1:25272598 G>A called '1':
            pair       -> Pair(RHD*01, RHD*09.01)
            genotypes  -> ['RHD*09.01/RHD*01N']
    """
    reference = (reference_alleles or {}).get(bg.type)
    reference_genotype = getattr(reference, "genotype", None)
    absent_slot = (gene_absent_subtypes or {}).get(bg.type, UNNAMED_SECOND_SLOT)
    missing_copy = bg.locus_copies == 1 and bg.chrom_copies == 2

    def make_list_of_lists(alleles):
        return [pair.genotypes for pair in alleles]

    def render(genotypes: list[str], co_existing: bool = False) -> str:
        """Join a pair's genotypes, collapsing a slot the sample does not have.

        pair.genotypes is already sorted by Pair._ordered, so it is not re-sorted here.

        co_existing says whether this pair came from the Knops path, and it is only
        consulted for the two-non-reference-alleles case, where the same shape means
        opposite things. Co-existing alleles sit on one chromosome by definition, so two
        of them at one gene copy is the ordinary result. Off that path it is impossible,
        and cant_have_2_non_ref_alleles_cuz_only_1_gene_copy has already excluded it by
        name - so reaching here means that filter did not run or did not catch it, which
        is worth saying rather than quietly writing the pair out as this used to.
        """
        if bg.chrom_copies == 1 and len(set(genotypes)) == 1:
            return f"{genotypes[0]}/{HAPLOID_SECOND_SLOT}"
        if missing_copy:
            # dict.fromkeys rather than a set: the duplicate has to go, but the order of
            # what is left has to stay reproducible. Pair.alleles is a frozenset, so the
            # only ordering guarantee here is the one Pair._ordered already applied.
            present = list(
                dict.fromkeys(
                    geno for geno in genotypes if geno != reference_genotype
                )
            )
            if len(present) < 2:
                # The ordinary shape: one real allele plus the reference slot that the
                # missing copy displaces, or the duplicate the pairing machinery writes
                # when one copy carries the allele and there is nothing to pair it with.
                return f"{(present or genotypes)[0]}/{absent_slot}"
            if not co_existing:
                logger.warning(
                    f"{bg.sample}: {bg.type} has one gene copy but a pair of two "
                    f"non-reference alleles ({'/'.join(present)}), which needs both on "
                    f"the one copy. It should have been excluded by "
                    f"cant_have_2_non_ref_alleles_cuz_only_1_gene_copy and was not, so "
                    f"it is being reported as called; see issue #40"
                )
        return "/".join(sorted(genotypes))

    if bg.alleles[AlleleState.CO] is not None:
        bg.genotypes = [
            render(co, co_existing=True)
            for co in make_list_of_lists(bg.alleles[AlleleState.CO])
        ]
    elif bg.alleles[AlleleState.NORMAL]:
        bg.genotypes = [
            render(normal_pair)
            for normal_pair in make_list_of_lists(bg.alleles[AlleleState.NORMAL])
        ]
    else:
        # Nothing paired. Where one slot was named and the other refused the strings are
        # already rendered and are the answer; where they are not this is the empty list
        # it has always been, which becomes 'Undetermined/Undetermined' downstream.
        bg.genotypes = list(bg.single_slot_genotypes)

    return bg


@apply_to_dict_values
def cant_revert_to_ref_cuz_a_passing_call_denies_it(
    bg: BloodGroup,
    vcf: VCF,
    df: pd.DataFrame,
    reference_alleles: dict[str, Allele],
) -> BloodGroup:
    """Decline to name a blood group whose reference a trusted call rules out.

    Rule 3 makes the reference the default when nothing is buildable, even where it
    needed a discarded variant, because that is the ISBT convention and what a lab
    scientist does by hand. This is the one shape where that default asserts something
    the caller positively denies.

    A reference allele at a lane locus is defined by '_ref' tokens, and a '_ref' token
    only exists where the sample has a reference copy. Where the alternate is homozygous
    there is no reference copy and no token, so the reference cannot be built at all -
    it is absent from the raw alleles rather than filtered out of them. Reporting it
    anyway claims wildtype at a position the caller called homozygous variant.

    The gate is that the denying call *passed*. Where the only thing standing between
    the sample and the reference is a call the caller doubted, rule 3 is exactly right
    and this does nothing: HG03600 GYPA is homozygous alternate at 144120554 with
    LowQual and heterozygous at the other two loci, so its reference is denied only by
    a doubted call and it keeps GYPA*02/GYPA*02.

    A reference can also be denied the other way round. Four of the database's 88
    references are defined partly by an *alternate* - KN*01, RHD*01, ABO*A1.01 and
    RHCE*01 - because at a lane locus the transcript reference differs from the genome
    reference. There the reference needs the alternate and a homozygous reference call
    denies it. That denial rests on the '_ref' token being present rather than on any
    FILTER value, because a '_ref' token has no FILTER of its own; absence of the token
    is no data and is never read as wildtype. HG04183 RHCE is the only instance across
    all nine datasets: homozygous reference at 25420739 where RHCE*01 needs 25420739_G_C,
    with its other route - 25408711 - denied only by a LowQual call.

    Deliberately narrow. It does not fire where the reference was *built and then
    struck* by FILTER - ABO, KN and RHD in the long read set - because there the
    reference needs the doubted variant itself and nothing contradicts it. That is
    messy but it is rule 3 working as agreed, not a defect.

    Args:
        bg (BloodGroup): The blood group, after process_genetic_data has supplied the
            reference pair.
        vcf (VCF): The sample's VCF, for which tokens exist at each locus.
        df (pd.DataFrame): The sample's rows, for each variant's FILTER.
        reference_alleles (dict[str, Allele]): The reference allele per blood group.

    Returns:
        BloodGroup: With the reference pair removed and recorded, or unchanged.

    Example:
        HG01871 GYPA. GYPA*01 needs three positions and is discarded over the LowQual
        call at 144120567. Nothing else is buildable, so the pool is empty and rule 3
        offers GYPA*02/GYPA*02 - which needs a reference copy at 144120554 and 144120555,
        both of them called 1/1 with PASS:

        4:144120554_C_A : Homozygous : PASS
        4:144120555_T_C : Homozygous : PASS
        4:144120567_A_G : Heterozygous : LowQual

        Note 144120567, the position that failed FILTER, is heterozygous and has its
        '_ref' token, so it is not what denies the reference. The two that do are both
        trusted.

        genotypes -> 'Undetermined/Undetermined', where it was 'GYPA*02/GYPA*02' (M-,N+).
    """
    pairs = bg.alleles.get(AlleleState.NORMAL) or []
    if len(pairs) != 1 or not pairs[0].all_reference:
        return bg
    # An empty pool is the rule 3 case. A pool with anything in it means some allele
    # survived and the reference pair was a real choice among others.
    if bg.variant_pool:
        return bg
    reference = reference_alleles.get(bg.type)
    if reference is None:
        return bg
    # Built and then struck is the other shape, and is left alone.
    if any(
        allele.genotype == reference.genotype
        for allele in bg.alleles.get(AlleleState.RAW) or []
    ):
        return bg

    for token in reference.defining_variants:
        if token in vcf.variants:
            continue
        locus = token.split("_")[0]
        if token.endswith("_ref"):
            # The reference wants wildtype here and the sample has none, so something
            # else is homozygous at this locus. It only counts if the caller vouched
            # for it.
            denied = any(
                called != token
                and called.split("_")[0] == locus
                and not variant_was_discarded(called, df)
                for called in vcf.variants
            )
        else:
            # The reference wants the alternate here - a lane locus, where the
            # transcript reference differs from the genome reference - and the sample
            # does not carry it. The '_ref' token has to actually be present to say so:
            # a locus nobody typed is no data, and reading absence as wildtype is the
            # mistake this whole filter exists to avoid. Where it is present the sample
            # is homozygous reference, since a heterozygote would carry both tokens,
            # and no FILTER is consulted because a '_ref' token has none of its own.
            denied = f"{locus}_ref" in vcf.variants
        if denied:
            bg.remove_pairs(
                list(pairs),
                "cant_revert_to_ref_cuz_a_passing_call_denies_it",
                reverts_to_reference=False,
            )
            return bg
    return bg


@apply_to_dict_values
def only_keep_alleles_if_FILTER_PASS(
    bg: BloodGroup, df: pd.DataFrame, no_filter: bool
) -> BloodGroup:
    """Keep only alleles whose every defining variant the caller vouched for.

    The one place upstream quality is acted on. Everything else is delegated upstream
    deliberately, so an allele needing a variant the caller doubted is dropped here and the
    blood group reverts to whatever the remaining alleles support - usually the reference
    allele. '_ref' tokens are skipped, having no FILTER of their own.

    Not every non-PASS value is a doubt about the call, which is why the test is a lookup
    rather than a comparison with 'PASS'. Some values report a statistic computed across a
    whole cohort, so on a jointly called file they say nothing about this sample; some mark
    which of several rows for a marker is the recommended one, which is probeset selection
    rather than call quality. filter_semantics carries the classification and the caller's
    own description of each value; an unrecognised value still excludes.

    Args:
        bg (BloodGroup): The BloodGroup object containing alleles to filter.
        df (pd.DataFrame): The sample's VCF rows, used to look up each variant's FILTER.
        no_filter (bool): Skip the check entirely and promote every raw allele, ie the
        --no_filter flag.

    Returns:
        BloodGroup: The BloodGroup with alleles[FILT] set, and anything dropped recorded
        under filtered_out['FILTER_not_PASS'].
    """
    if no_filter:
        bg.alleles[AlleleState.FILT] = list(bg.alleles[AlleleState.RAW])
        return bg
    passed_filtering = []
    unclassified: set[str] = set()
    order_dependent: dict[str, list[str]] = {}
    for allele in bg.alleles[AlleleState.RAW]:
        keeper = True
        for variant in allele.defining_variants:
            if "_ref" in variant:
                continue
            vcf_var = allele.big_variants.get(variant, variant)
            # loci = vcf_var.split("_")[0]
            filter_values = filter_values_for(vcf_var, df)
            if not filter_values:
                message = f"FILTER parsing failed. Sample: {bg.sample}, BG: {bg.type}, variant/s: {variant}"
                logger.error(message)
                raise IndexError(message)
            if rows_disagree_about_exclusion(filter_values):
                order_dependent[vcf_var] = filter_values
            filter_value = filter_values[0]
            excludes, unknown = filter_excludes_allele(filter_value)
            unclassified.update(unknown)
            if excludes:
                keeper = False
                break
        if keeper:
            passed_filtering.append(allele)

    for value in sorted(unclassified):
        logger.warning(
            f"FILTER value '{value}' is not in filter_values.tsv, so it was read as a "
            f"reason to exclude. Sample: {bg.sample}, BG: {bg.type}. If it does not mean "
            "the call is doubtful, classify it there."
        )
    warn_if_the_row_order_decided_it(bg, order_dependent)

    bg.filtered_out["FILTER_not_PASS"] = [
        allele
        for allele in bg.alleles[AlleleState.RAW]
        if allele not in passed_filtering
    ]
    bg.alleles[AlleleState.FILT] = passed_filtering
    warn_if_critical_variant_not_trusted(bg, df)

    return bg


def get_fully_homozygous_alleles(
    ranked_chunks: list[list[Allele]],
    variant_pool: dict[str, Any],
    chrom_copies: int = 2,
) -> list[list[Allele]]:
    """Filter out alleles that are not fully homozygous from a list of ranked allele chunks.

    Uses a partial function to check each allele's variants in the provided `variant_pool`.
    Only alleles where every relevant variant equals the required homozygous genotype (2)
    are included in the result.

    The question being asked is 'is this token on every chromosome the sample has', which
    was 2 while every region was assumed diploid. Where chrom_copies is 1 the answer is 1
    - a Hemizygous token on a single-copy region is on every chromosome present, and is
    fully homozygous in the only sense available. The comparison is against chrom_copies
    rather than a literal so that variant_pool_numeric stays truthful: Zygosity.HEM keeps
    scoring 1, because one chromosome carries it.

    Args:
        ranked_chunks (list[list[Allele]]):
            A list of lists (chunks), where each chunk contains ranked Allele objects.
        variant_pool (dict[str, Any]):
            A dictionary containing variant data used for assessing homozygosity.
            The exact structure depends on the `check_available_variants` function.
        chrom_copies (int):
            Copies of the region the sample was born with. Defaults to 2.

    Returns:
        list[list[Allele]]:
            A list of lists, each mirroring the structure of `ranked_chunks`
            but including only alleles that are fully homozygous in every variant.

    Raises:
        KeyError:
            If a variant key is missing in `variant_pool`.
    """
    check_hom = partial(
        check_available_variants, chrom_copies, variant_pool, operator.eq
    )
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
                reference_alleles,
                bg.variant_pool_numeric,
                bg.alleles[AlleleState.FILT],
                bg.chrom_copies,
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
        homs = get_fully_homozygous_alleles(
            ranked_chunks, bg.variant_pool_numeric, bg.chrom_copies
        )

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
        homs = get_fully_homozygous_alleles(
            self.ranked_chunks, bg.variant_pool_numeric, bg.chrom_copies
        )
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
                    bg.chrom_copies,
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
        homs = get_fully_homozygous_alleles(
            ranked_chunks, bg.variant_pool_numeric, bg.chrom_copies
        )
        for ranked_homs in homs:
            for hom in ranked_homs:
                pair = Pair(allele1=hom, allele2=hom)
                if pair not in bg.alleles[AlleleState.NORMAL]:
                    bg.filtered_out["excluded_due_to_rank_hom"].append(pair)

    return bg


def make_pair(
    reference_alleles: dict[str, str],
    variant_pool: list[str],
    sub_results: list[str],
    chrom_copies: int = 2,
) -> list[str]:
    """Creates a pair of alleles based on the given parameters.

    Where chrom_copies is 1 the allele is always duplicated rather than paired with the
    reference, because the single chromosome carries it and there is no second chromosome
    to be reference. That duplicate is what get_genotypes renders as 'ALLELE/-'. Taking
    the reference branch instead is what made a hemizygous male read as a heterozygous
    female before v2.4.4 - see issue #40.

    Args:
        reference_alleles (Dict[str, str]): A mapping from blood group to reference
        allele.
        variant_pool (List[str]): A list of available variants.
        sub_results (List[str]): A list containing the initial results, expected to be
        of length 1.
        chrom_copies (int): Copies of the region the sample was born with. Defaults to 2.

    Returns:
        List[str]: A list containing the original results and an additional allele,
        either a duplicate of the first (if checks pass) or a corresponding
        reference allele.

    Raises:
        AssertionError: If the length of `sub_results` is not 1.
    """
    sub_results = list(sub_results)
    check_vars = partial(
        check_available_variants, chrom_copies, variant_pool, operator.eq
    )
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


# @apply_to_dict_values
# def add_CD_to_XG(bg: BloodGroup) -> BloodGroup:
#     """TODO why not just use the CD99 vars??
#     adds CD to XG blood group.

#     Args:
#         bg (BloodGroup): The BloodGroup object to be processed.

#     Returns:
#         BloodGroup: The processed BloodGroup object.
#     """
#     if bg.genotypes == ["XG*01/XG*01"]:
#         bg.genotypes = ["XG*01/XG*01", "CD99*01/CD99*01"]
#     return bg


def add_refs(
    db: Db,
    res: dict[str, BloodGroup],
    exclude,
    vcf: VCF | None = None,
) -> dict[str, BloodGroup]:
    """Add reference genotypes to existing results or create new entries for them.

    A blood group only reaches here with no alleles of its own, so the genotype string is
    built directly rather than by get_genotypes. That means both single copy renderings
    have to be repeated here, and this is the case with no evidence left in the variant
    pool to work either out from - remove_home_ref dropped the only rows.

    A male with no XK variant at all is hemizygous *reference*, which is 'XK*01/-' and not
    'XK*01/XK*01'. State table row B2.

    A gene reported at one copy whose every locus was wildtype is 'RHD*01/RHD*01N' and not
    'RHD*01/RHD*01'. State table row D3, and the common shape rather than an edge case: a
    sample missing one copy of a gene is usually ordinary at the loci that remain, so every
    row is a haploid '0' and every one of them gets dropped.

    Args:
        db (Db): The database object containing reference alleles.
        res (Dict[str, BloodGroup]): Dictionary of BloodGroup objects to be updated
        with reference data.
        exclude: Blood groups to skip.
        vcf (VCF | None): The VCF this sample came from, carrying the ploidy evidence.
        None means every blood group gets two slots, which is the behaviour before v2.4.4.

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
            bg = BloodGroup(
                type=blood_group,
                alleles={
                    AlleleState.RAW: [reference],
                    AlleleState.FILT: [reference],
                    AlleleState.NORMAL: [Pair(*[reference] * 2)],
                },
                sample="ref",
            )
            if vcf is not None:
                bg.chrom_copies = chrom_copies_for_bg(bg, vcf, db.single_copy_types)
                bg.locus_copies = locus_copies_for_bg(bg, vcf, db.loci_by_type)

            if bg.chrom_copies == 1:
                second_slot = HAPLOID_SECOND_SLOT
            elif bg.locus_copies == 1:
                second_slot = db.gene_absent_subtypes.get(
                    blood_group, UNNAMED_SECOND_SLOT
                )
            else:
                second_slot = reference.genotype

            bg.genotypes = [f"{reference.genotype}/{second_slot}"]
            res[blood_group] = bg
    return res


@apply_to_dict_values
def only_keep_alleles_if_best_big_del(
    bg: BloodGroup, matches: list[MatchResult]
) -> BloodGroup:
    """Keep alleles only if they contain best-matching large deletion variants.

    Filters alleles to retain only those containing large deletion variants that
    are among the best matches from variant calling. Updates alleles with their
    corresponding best-match large variant information.

    Args:
        bg: BloodGroup object containing alleles to filter and update.
        matches: List of MatchResult objects from variant matching, containing
            variant call information and database references.

    Returns:
        The modified BloodGroup object with alleles updated to include large
        variant information for best matches.

    Note:
        Uses a very tight tolerance (1e-9) for selecting best matches to ensure
        only the highest quality variant calls are considered.
    """

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
