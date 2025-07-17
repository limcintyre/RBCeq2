from __future__ import annotations

from rbceq2.core_logic.alleles import Allele, BloodGroup, Pair
from rbceq2.core_logic.utils import (
    Zygosity
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

def all_hom(bg: BloodGroup, current_allele: Allele) -> bool:
        """if > 1 het vars in an alleles defiing variant set, this shouldnt apply
        unless phased, which is handled elsewhere - TODO ensure its handled elsewhere"""

        return all(
            bg.variant_pool.get(allele_var) == Zygosity.HOM
            for allele_var in current_allele.defining_variants
        )


def check_hom(variant_pool: dict[str, str], current_allele: Allele) -> bool:
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



def check_phase(variant_pool: dict[str, str], current_allele: Allele) -> bool:
    """
    True if all same phase set and or HOM
    """

    phase_sets = [
        phase
        for variant, phase in variant_pool.items()
        if variant in current_allele.defining_variants and phase != "1/1"
    ]

    if not phase_sets and "ABO" in current_allele.genotype:
        return True  # TODO still needed?
    return len(set(phase_sets)) == 1


def check_phase_set(variant_pool: dict[str, str], current_allele: Allele) -> bool:
    """
    True if all same phase set and or HOM
    """
    phase_sets = [
        phase
        for variant, phase in variant_pool.items()
        if variant in current_allele.defining_variants and phase != "."
    ]

    return len(set(phase_sets)) == 1


# def _check_hom(variant_pool: dict[str, str], current_allele: Allele) -> bool:
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
#         for variant, zygo in variant_pool.items()
#         if variant in current_allele.defining_variants
#     ]

#     return all(item == Zygosity.HOM for item in variants)

