from __future__ import annotations
from collections.abc import Iterable
from rbceq2.core_logic.constants import AlleleState, SYNTHESISED_HOM_REF_GT

from rbceq2.core_logic.alleles import Allele, BloodGroup, Pair
from rbceq2.core_logic.utils import Zygosity


def flatten_alleles(pairs: list[Pair]) -> set[Allele]:
    """Flatten the pairs into a set of alleles.

    Args:
        pairs (list[Pair]): A list of Pair objects, where each Pair is an
            iterable containing Allele objects.

    Returns:
        set[Allele]: A set containing all Allele objects from the given pairs.
    """
    return {allele for pair in pairs for allele in pair}


def all_hom(variant_pool: dict[str, str], current_allele: Allele) -> bool:
    """All variants are homozygous"""

    return all(
        variant_pool.get(allele_var) == Zygosity.HOM
        for allele_var in current_allele.defining_variants
    )


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


def _het_phase_evidence(bg: BloodGroup, variant: str) -> tuple[str, str, str] | None:
    """Return chromosome, known phase set and orientation for a phased HET token.

    Homozygous and hemizygous tokens do not locate a heterozygous allele on one
    of two chromosomes. Unknown phase sets cannot connect two orientations.
    This reads existing normalized evidence without assigning new phase sets.

    Args:
        bg (BloodGroup): Existing variant and phase pools.
        variant (str): Defining token with chromosome and position.

    Returns:
        tuple[str, str, str] | None: Chromosome, phase set and HET orientation,
        or None when the token supplies no usable block evidence.
    """
    if bg.variant_pool.get(variant) != Zygosity.HET:
        return None
    phase = bg.variant_pool_phase.get(variant)
    phase_set = bg.variant_pool_phase_set.get(variant)
    if phase not in {"0|1", "1|0"}:
        return None
    if not isinstance(phase_set, str) or not carries_phase(phase_set):
        return None
    chrom, separator, _ = variant.partition(":")
    if not separator or not chrom:
        return None
    return chrom.removeprefix("chr"), phase_set, phase


def _het_phase_summary(
    bg: BloodGroup, variants: Iterable[str]
) -> tuple[str, str, str] | None:
    """Return one shared block and orientation supported by every HET definition.

    Each heterozygous defining token must carry usable evidence. Homozygous and
    hemizygous definitions are ignored only for orientation; copy counts stay intact.
    A collection without HET definitions supplies no such placement evidence.

    Args:
        bg (BloodGroup): Existing variant and phase pools.
        variants (Iterable[str]): Defining tokens whose HET evidence must agree.

    Returns:
        tuple[str, str, str] | None: Shared evidence, or None when any HET token
        is unknown, the tokens disagree, or there are no HET definitions.
    """
    evidence = set()
    for variant in variants:
        if bg.variant_pool.get(variant) != Zygosity.HET:
            continue
        current = _het_phase_evidence(bg, variant)
        if current is None:
            return None
        evidence.add(current)
    return evidence.pop() if len(evidence) == 1 else None


def identify_unphased(bg: BloodGroup, alleles: list[Allele]) -> list[Allele]:
    """Find definitions contradicted by opposite HET sides within one known block.

    Args:
        bg (BloodGroup): Existing zygosity, orientation and phase-set evidence.
        alleles (list[Allele]): Candidate definitions to inspect.

    Returns:
        list[Allele]: Definitions requiring variants on opposite haplotypes in
        the same chromosome/phase set. Callers record their named exclusions.
    """
    to_remove = []
    for allele in alleles:
        sides_by_block = {}
        for variant in allele.defining_variants:
            evidence = _het_phase_evidence(bg, variant)
            if evidence is None:
                continue
            block, side = evidence[:2], evidence[2]
            sides_by_block.setdefault(block, set()).add(side)
        if any(len(sides) > 1 for sides in sides_by_block.values()):
            to_remove.append(allele)
    return to_remove


def proceed(bg: BloodGroup, allele_state: AlleleState) -> bool:
    """Some filter functions iterate over AlleleStates
    ie for allele_state in [AlleleState.NORMAL, AlleleState.CO]:

    if CO, we only want to process knops and then only if there are
    instances of co_exsisting alleles"""

    if allele_state != AlleleState.CO:
        return True
    if bg.type != "KN" or bg.alleles[allele_state] is None:
        return False
    return True


def check_var(
    bg: BloodGroup, pair: Pair, allele_state: AlleleState, variant: str
) -> int:
    """checks that HET vars are used"""
    allele1_vars_plus_het_var = set(pair.allele1.defining_variants) | {variant}
    allele2_vars_plus_het_var = set(pair.allele2.defining_variants) | {variant}
    flattened = {
        allele
        for pair2 in bg.alleles[allele_state]
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
    return hits
