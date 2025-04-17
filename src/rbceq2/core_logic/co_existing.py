import operator
from functools import partial
from itertools import combinations, combinations_with_replacement
from typing import Iterable

from rbceq2.core_logic.alleles import Allele, BloodGroup, Pair
from rbceq2.core_logic.constants import LOW_WEIGHT, AlleleState
from rbceq2.core_logic.utils import (
    apply_to_dict_values,
    check_available_variants,
    chunk_geno_list_by_rank,
)


def sub_alleles(lst: tuple[Allele], allele_relationship: dict[str, bool]) -> bool:
    """look up a precomputed relationship with a special key - KN only

    Args:
        lst (list): A list of strings to be checked.

    Returns:
        list: A list of booleans indicating whether each element is in another element.
    """
    for allele1 in lst:
        for allele2 in lst:
            key = f"{allele1.genotype}_isin_{allele2.genotype}"
            if allele1 == allele2:
                continue
            if allele_relationship[key]:
                return True
    return False


def all_hom_variants(all_homs: list[Allele], current_combo: tuple[Allele, ...]) -> bool:
    """
    Check if all homozygous alleles are present in the current combination.

    This function verifies whether all alleles designated as homozygous are
    present in the given combination of alleles.

    Args:
        all_homs (list[Allele]): A list of homozygous alleles to be checked. All_homs
        is a set of Alleles where all defining vars are HOM
        current_combo (tuple[Allele, ...]): A tuple of alleles representing the current
        combination.

    Returns:
        bool: True if the count DOES NOT match.

    Example:

    'possible': [
        Allele(genotype='KN*01.06',
                genotype_alt='.',
                defining_variants=frozenset({'1:207782856_A_G',
                                            '1:207782889_A_G',
                                            '1:207782916_A_T', hom
                                            '1:207782931_A_G'}), hom
                weight_geno=1000,
                reference=False,
        Allele(genotype='KN*01.07',
                genotype_alt='.',
                defining_variants=frozenset({'1:207782889_A_G',
                                            '1:207782916_A_T'}), hom
                weight_geno=1000,
                reference=False,
        Allele(genotype='KN*01.10',
                genotype_alt='.',
                defining_variants=frozenset({'1:207782916_A_T', hom
                                            '1:207782931_A_G'}), hom
                weight_geno=1000,
                reference=False,

    sample='140',
    variant_pool={'1:207782856_A_G': 'Heterozygous',
                '1:207782889_A_G': 'Heterozygous',
                '1:207782916_A_T': 'Homozygous',
                '1:207782931_A_G': 'Homozygous'},
    genotypes=['KN*01.06/KN*01.10'],
    phenotypes=[],
    """
    hom_count = 0
    for hom_allele in all_homs:
        for allele in current_combo:
            if hom_allele == allele:
                hom_count += 1
            if hom_allele in allele:
                hom_count += 1
    return hom_count != len(all_homs)  # not less OR MORE


def mushed_vars(mushed_combo: Iterable[Allele]) -> set[str]:
    """Combine the defining_variants of all alleles in mushed_combo into a
    single set.

    Args:
        mushed_combo (Iterable[Allele]): An iterable of Allele objects whose
            defining_variants will be combined.

    Returns:
        set[str]: A set containing all unique defining_variants from the given
            Allele objects.
    """
    return set().union(
        *[mushed_allele.defining_variants for mushed_allele in mushed_combo]
    )


def can_co_exist(mushed_combo1: set[str], tmp_pool23: dict[str, int]) -> bool:
    """Check if the defining variants in 'mushed_combo1' can be satisfied by the
    'tmp_pool23' pool of available variants (count>0).

    Args:
        mushed_combo1 (set[str]): A set of defining variants.
        tmp_pool23 (dict[str, int]): A dictionary mapping variant names to their
            available counts.

    Returns:
        bool: True if all variants in mushed_combo1 are available in tmp_pool23
            with count > 0, False otherwise.
    """
    available_variants: set[str] = {
        variant for variant, count in tmp_pool23.items() if count > 0
    }
    intersection = mushed_combo1.intersection(available_variants)

    return len(intersection) == len(mushed_combo1)


def decide_if_co_existing(
    tmp_pool2: dict[str, int],
    combo1: tuple[Allele, ...],
    combo2: tuple[Allele, ...],
    mushed_combo1: set[str],
    co_existing: list[tuple[tuple[Allele, ...], tuple[Allele, ...]]],
) -> list[tuple[tuple[Allele, ...], tuple[Allele, ...]]]:
    """
    Decide if combo1 and combo2 can co-exist based on the variant pool and
    already discovered co-existing combos.

    Args:
        variant_pool: A dictionary of variant -> numeric count.
        combo1: A tuple of Alleles forming one combination.
        combo2: A tuple of Alleles forming the second combination.
        combined_variants: A set of defining variants from combo1 (used for can_co_exist).
        existing_pairs: The current list of co-existing combos to which we may add.

    Returns:
        The possibly updated list of existing_pairs.
    """
    check_vars_other_strand = partial(
        check_available_variants, 0, tmp_pool2, operator.gt
    )
    if all(
        all(check_vars_other_strand(allele2)) for allele2 in combo1
    ) and can_co_exist(mushed_combo1, tmp_pool2):
        # ie can they exist given other strand
        # if can_co_exist(mushed_combo1, tmp_pool2):
        combo_pair = (combo1, combo2)
        if combo_pair not in co_existing:
            co_existing.append(combo_pair)  # TODO - order for testing?

        # TODO sort on order in ISBT
    return co_existing


@apply_to_dict_values
def homs(
    bg: BloodGroup,
) -> BloodGroup:
    """Check all combinations of alleles against all other combinations to see if
    they can co-exist once the variant pool is reduced by the respective
    defining variants.

    Args:
        bg (BloodGroup): A BloodGroup object containing alleles, type, variant pool,
            and misc information.

    Returns:
        BloodGroup: The updated BloodGroup with modified alleles and misc values.
    """
    if bg.type != "KN":
        bg.alleles[AlleleState.CO] = None
        return bg
    unique_alleles = set(bg.alleles[AlleleState.POS])

    check_hom_vars = partial(
        check_available_variants, 2, bg.variant_pool_numeric, operator.eq
    )
    bg.misc = {}  # TODO misc is a bit promiscuous
    bg.misc["homs"] = {
        allele
        for allele in unique_alleles
        if all(check_hom_vars(allele))
        if not allele.reference
    }

    return bg


@apply_to_dict_values
def max_rank(
    bg: BloodGroup,
) -> BloodGroup:
    """Check all combinations of alleles against all other combinations to see if
    they can co-exist once the variant pool is reduced by the respective defining
    variants.

    Args:
        bg (BloodGroup): A BloodGroup object containing alleles, type, variant pool,
            and misc information.

    Returns:
        BloodGroup: The updated BloodGroup with modified alleles and misc values.
    """
    if bg.type != "KN":
        bg.alleles[AlleleState.CO] = None
        return bg
    unique_alleles = set(bg.alleles[AlleleState.POS])

    trumpiest = [
        allele.weight_geno
        for allele in list(bg.misc["homs"])
        + [a for a in unique_alleles if a.number_of_defining_variants == 1]
    ]

    bg.misc["max_rank"] = min(trumpiest) if trumpiest else LOW_WEIGHT

    return bg


@apply_to_dict_values
def prep_co_putative_combos(
    bg: BloodGroup,
    allele_relationships: dict[str, dict[str, bool]],
) -> BloodGroup:
    """Check all combinations of alleles against all other combinations
    to see if can co-exist , once the variant pool is reduced
    by the respective defining variants"""

    def make_allele_combos(
        flattened_alleles: list[Allele],
        homs: list[Allele],
        allele_relationship: dict[str, bool],
    ) -> list[tuple[Allele, ...]]:
        """Generate combinations of alleles, filtering out those that do not meet
        criteria.

        This function generates combinations of alleles from the provided list,
        filtering out combinations based on the presence of homozygous alleles,
        sub-alleles, and reference alleles.

        Args:
            flattened_alleles (List[Allele]): List of all alleles to consider.
            homs (List[Allele]): List of homozygous alleles.
            allele_relationship (Dict[str, bool]): Dictionary of allele relationships.
            bg_type (str): The blood group type.

        Returns:
            Set[Tuple[Allele, ...]]: Set of allele combinations that meet the criteria.
        """
        combos = []
        for i in range(1, len(flattened_alleles) + 1)[::-1]:
            for combo1 in combinations(flattened_alleles, i):
                if all_hom_variants(homs, combo1):
                    continue
                if sub_alleles(combo1, allele_relationship):
                    continue
                if any(a.reference for a in combo1):
                    continue
                if len({a.sub_type for a in combo1}) > 1:
                    continue  # TODO can't combine subtypes - have asked Eileen
                ranked_combo1 = tuple(chunk_geno_list_by_rank(combo1)[0])

                assert ranked_combo1 not in combos
                combos.append(ranked_combo1)
        return combos

    if bg.type != "KN":
        bg.alleles[AlleleState.CO] = None
        return bg

    unique_alleles = sorted(
        set(bg.alleles[AlleleState.POS]), key=lambda allele: allele.genotype
    )
    bg.misc["combos"] = make_allele_combos(
        unique_alleles, bg.misc["homs"], allele_relationships[bg.type]
    )

    return bg


@apply_to_dict_values
def add_co_existing_alleles(
    bg: BloodGroup,
) -> BloodGroup:
    """Check all combinations of alleles against all other combinations
    to see if can co-exist , once the variant pool is reduced
    by the respective defining variants

    Identify and store combinations of alleles that can co-exist in a BloodGroup.

    This function is specifically designed for blood groups of type "KN". It examines
    all pairwise combinations (including self-combinations) in `bg.misc["combos"]`.
    For each pair, the function reduces the available variant pool by the defining
    variants of the second combination, then decides whether both sets of alleles
    can co-exist without conflict. If so, the combination is added to the
    `AlleleState.CO` list.

    Args:
        bg (BloodGroup):
            The BloodGroup object whose alleles are being checked for co-existence.

    Returns:
        BloodGroup:
            The input BloodGroup object with its `AlleleState.CO` attribute updated
            to either a list of co-existing allele combinations (if any) or `None`
            (if no viable combinations are found).

    Raises:
        KeyError:
            If an expected variant is missing from `bg.variant_pool_numeric`.
    """

    if bg.type != "KN":
        bg.alleles[AlleleState.CO] = None
        return bg
    co_existing = []
    for combo1, combo2 in combinations_with_replacement(bg.misc["combos"], 2):
        # This includes pairs like (A, A), (A, B), (B, B), but not both (A, B) and (B, A)
        mushed_combo1 = mushed_vars(combo1)
        tmp_pool2 = bg.variant_pool_numeric
        for variant_on_other_strand in mushed_vars(combo2):
            tmp_pool2[variant_on_other_strand] -= 1
        co_existing = decide_if_co_existing(
            tmp_pool2, combo1, combo2, mushed_combo1, co_existing
        )
        # TODO else - record why combo failed in bg
    if co_existing:
        bg.alleles[AlleleState.CO] = co_existing
    else:
        bg.alleles[AlleleState.CO] = None
    return bg


@apply_to_dict_values
def add_co_existing_allele_and_ref(
    bg: BloodGroup,
    reference_alleles: dict[str, Allele],
) -> BloodGroup:
    """Check all combinations of alleles against all other combinations
    to see if can co-exist , once the variant pool is reduced
    by the respective defining variants

    Evaluate allele combinations against a reference allele for co-existence.

    This function also targets blood groups of type "KN". For each combination
    in `bg.misc["combos"]`, it checks whether the allele set (combo1) can coexist
    with the specified reference allele (drawn from `reference_alleles`). It does so
    by reducing the `variant_pool_numeric` according to the reference allele's defining
    variants, then verifying if both sets remain viable.

    Args:
        bg (BloodGroup):
            The BloodGroup object whose alleles are being tested for co-existence
            with a reference allele.
        reference_alleles (dict[str, Allele]):
            A dictionary mapping blood group types to their corresponding reference
            Allele object.

    Returns:
        BloodGroup:
            The updated BloodGroup object. The `AlleleState.CO` attribute is assigned
            a list of (combination, reference) pairs that can coexist, or `None` if
            no valid pairs are found.

    Raises:
        KeyError:
            If an expected variant is missing from `bg.variant_pool_numeric` or
            `reference_alleles`.
    """
    if bg.type != "KN":
        bg.alleles[AlleleState.CO] = None
        return bg

    if bg.alleles[AlleleState.CO] is not None:
        co_existing = bg.alleles[AlleleState.CO]
    else:
        co_existing = []
    for combo1 in bg.misc["combos"]:
        mushed_combo1 = mushed_vars(combo1)
        if bg.misc["homs"]:
            continue
        if bg.misc["max_rank"] < min([allele.weight_geno for allele in combo1]):
            continue
        tmp_pool2 = bg.variant_pool_numeric
        ref = reference_alleles[bg.type]
        # some refs have defining vars (ie ABO*A1.01 has
        # 261delG/136132908_T_TC)
        for variant_on_other_strand in ref.defining_variants:
            if variant_on_other_strand in tmp_pool2:
                tmp_pool2[variant_on_other_strand] -= 1
        co_existing = decide_if_co_existing(
            tmp_pool2, combo1, (ref,), mushed_combo1, co_existing
        )
    if co_existing:
        bg.alleles[AlleleState.CO] = co_existing
    else:
        bg.alleles[AlleleState.CO] = None
    return bg


@apply_to_dict_values
def filter_redundant_pairs(bg: BloodGroup) -> BloodGroup:
    """
    Filter out duplicate pairs, considering both original and reversed forms.

    Args:
        co_existing (list[str]): List of string pairs separated by '/'.

    Returns:
        List of unique combo pairs.
    """

    if bg.alleles[AlleleState.CO] is None:
        return bg
    seen = set()
    result = []
    for combo1, combo2 in bg.alleles[AlleleState.CO]:
        geno1, geno2 = geno_str(combo1), geno_str(combo2)
        fwd = f"{geno1}/{geno2}"
        rev = f"{geno2}/{geno1}"
        if fwd not in seen and rev not in seen:
            seen.update({fwd, rev})
            result.append((combo1, combo2))
    bg.alleles[AlleleState.CO] = result

    return bg


@apply_to_dict_values
def mush(bg: BloodGroup) -> BloodGroup:
    """
    Mushes co-existing alleles into a new allele:
      - Iterates over each 'combo' in bg.alleles[AlleleState.CO]
      - For each sub-combo in combos:
          if len(combo) > 1 -> we unify numeric + alphanumeric phenotypes
                               and produce a new "mushed" Allele
          if len(combo) == 1 -> we keep the single allele as-is
      - Then each mushed_pair must have exactly two Alleles to form a Pair.
      - Raises ValueError if there's a conflict or if the pair doesn't have exactly two.

    Args:
        bg: A BloodGroup object whose 'CO' (co-existing) combos we want to merge.

    Returns:
        The same BloodGroup, with bg.alleles[AlleleState.CO] replaced by
        new Pairs of Alleles or set to None if none remain.

    Raises:
        ValueError: if any conflict arises or a pair isn't exactly two Alleles.
    """

    def _check_conflict(
        expressions: dict[str, str], name: str, expr: str
    ) -> dict[str, str]:
        """Raise ValueError if 'name' is already in 'expressions' with a different expr."""
        if name in expressions and expressions[name] != expr:
            raise ValueError(
                f"Conflicting expressions for antigen '{name}': "
                f"'{expressions[name]}' and '{expr}'"
            )
        expressions[name] = expr
        return expressions

    result = []
    if bg.alleles[AlleleState.CO] is None:
        return bg
    for combos in bg.alleles[AlleleState.CO]:
        mushed_pair = []
        for combo in combos:
            if len(combo) > 1:
                # Process multiple alleles in a combo
                geno_numeric = geno_str(combo)

                # Process pheno_numeric
                phenotypes = [allele.phenotype.split(":")[1] for allele in combo]
                antigens = []
                for phenotype in phenotypes:
                    antigens.extend(phenotype.split(","))

                antigen_expressions: dict[str, str] = {}
                for antigen in antigens:  # TODO should I use Ant objects here!!!!????
                    if antigen.startswith("-"):
                        expression, name = "-", antigen[1:]
                    else:
                        expression, name = "+", antigen
                    antigen_expressions = _check_conflict(
                        antigen_expressions, name, expression
                    )

                unique_antigens = []
                for name, expression in antigen_expressions.items():
                    if expression == "-":
                        unique_antigens.append(f"-{name}")
                    else:
                        unique_antigens.append(name)
                unique_antigens.sort(key=lambda x: x.lstrip("-"))

                pheno_numeric = f"{bg.type}:{','.join(unique_antigens)}"
                if "1w," in pheno_numeric and "1," in pheno_numeric:
                    pheno_numeric = pheno_numeric.replace("1w,", "")  # Helgeson

                # Process pheno_alphanumeric
                phenotypes_alt = [allele.phenotype_alt for allele in combo]
                antigens_alt = []
                for phenotype_alt in phenotypes_alt:
                    antigens_alt.extend(phenotype_alt.split(","))

                antigen_expressions_alt: dict[str, str] = {}
                for antigen in antigens_alt:
                    if "-" in antigen:
                        expression = "-"
                        name = antigen.replace("-", "")
                    else:
                        expression = "+"
                        name = antigen.replace("+", "")
                    antigen_expressions_alt = _check_conflict(
                        antigen_expressions_alt, name, expression
                    )

                unique_antigens_alt = []
                for name, expression in antigen_expressions_alt.items():
                    if expression == "-":
                        if ")" in name:
                            unique_antigens_alt.append(name.replace(")", "-)"))
                        else:
                            unique_antigens_alt.append(f"{name}-")
                    elif ")" in name:
                        unique_antigens_alt.append(name.replace(")", "+)"))
                    else:
                        unique_antigens_alt.append(f"{name}+")
                unique_antigens_alt.sort(key=lambda x: x.rstrip("+-"))

                pheno_alphanumeric = ",".join(unique_antigens_alt)
                if "Kn(aw+)," in pheno_alphanumeric and "Kn(a+)," in pheno_alphanumeric:
                    pheno_alphanumeric = pheno_alphanumeric.replace(
                        "Kn(aw+),", ""
                    )  # Helgeson
                mushed_pair.append(
                    make_mushed_allele(
                        combo, geno_numeric, pheno_numeric, pheno_alphanumeric
                    )
                )
            elif len(combo) == 1:
                # Single allele, append as is
                mushed_pair.append(combo[0])

        # Ensure mushed_pair has exactly two alleles
        if len(mushed_pair) == 2:
            result.append(Pair(*mushed_pair))
        else:
            # Handle cases where there's only one allele
            # You can decide how to handle this based on your requirements
            # For now, let's raise an error
            raise ValueError("Each mushed_pair must contain exactly two alleles.")
    if result:
        bg.alleles[AlleleState.CO] = result
    else:
        bg.alleles[AlleleState.CO] = None

    return bg


def geno_str(combo: tuple[Allele, ...]) -> str:
    """Build the genotype string by sorting on allele.genotype and joining with '+'.

    Args:
        combo (tuple[Allele, ...]): A tuple of Allele objects.

    Returns:
        str: A string with sorted allele.genotype values concatenated by '+'.
    """
    return "+".join(sorted([a.genotype for a in combo]))


def make_mushed_allele(
    combo: tuple[Allele, ...], geno: str, pheno_numeric: str, pheno_alpha: str
) -> Allele:
    """Mush co existing alleles into one new allele.

    Args:
        combo (tuple[Allele, ...]): A tuple of Allele objects to be merged.
        geno (str): The genotype string for the new allele.
        pheno_numeric (str): The numeric phenotype string for the new allele.
        pheno_alpha (str): The alphanumeric phenotype string for the new allele.

    Returns:
        Allele: A new Allele object with merged attributes from the input alleles.
    """
    for pheno in [pheno_numeric, pheno_alpha]:
        assert ",," not in pheno
    return Allele(
        genotype=geno,
        phenotype=pheno_numeric,
        genotype_alt="mushed",
        phenotype_alt=pheno_alpha,
        defining_variants=frozenset.union(*[a.defining_variants for a in combo]),
        weight_geno=min([a.weight_geno for a in combo]),
        weight_pheno=min([a.weight_pheno for a in combo]),
        sub_type="/".join({a.sub_type for a in combo}),
    )


@apply_to_dict_values
def list_excluded_co_existing_pairs(
    bg: BloodGroup,
    reference_alleles: dict[str, Allele],
) -> BloodGroup:
    """Check all combinations of alleles against all other combinations to see if
    they can co-exist once the variant pool is reduced by the respective defining
    variants.

    Args:
        bg (BloodGroup): A BloodGroup object containing alleles, type, and misc
            information.
        reference_alleles (dict[str, Allele]): A dictionary mapping blood group types
            to their reference Allele.

    Returns:
        BloodGroup: The updated BloodGroup with filtered out co-existing pairs.
    """
    if bg.type != "KN":
        bg.alleles[AlleleState.CO] = None
        return bg
    tested = []
    ref = reference_alleles[bg.type]
    for combo1 in bg.misc["combos"]:
        for combo2 in bg.misc["combos"]:
            tested.append(
                Pair(
                    make_mushed_allele(combo1, geno_str(combo1), "numeric", "alpha"),
                    make_mushed_allele(combo2, geno_str(combo2), "numeric", "alpha"),
                )
            )
        tested.append(
            Pair(ref, make_mushed_allele(combo1, geno_str(combo1), "numeric", "alpha"))
        )

    bg.filtered_out[AlleleState.CO] = [
        pair for pair in tested if pair not in bg.alleles[AlleleState.CO]
    ]

    return bg
