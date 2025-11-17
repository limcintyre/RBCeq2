from __future__ import annotations

import copy
from collections import Counter, defaultdict
from functools import partial
from itertools import product
from typing import TYPE_CHECKING
from venv import logger

from rbceq2.core_logic.constants import (
    ANTITHETICAL,
    AlleleState,
    BgName,
    PhenoType,
    RHD_ANT_MAP,
)
from rbceq2.core_logic.data_procesing import apply_to_dict_values
from rbceq2.core_logic.utils import (
    BeyondLogicError,
)
from icecream import ic
from typing import Mapping
from rbceq2.db.db import compare_antigen_profiles
from rbceq2.core_logic.alleles import Allele

from . import antigens as an

if TYPE_CHECKING:
    import pandas as pd

    from rbceq2.core_logic.alleles import BloodGroup, Pair


def fut_helper(res: dict[str, BloodGroup]):
    """Generate all possible pairs of alphanumeric phenotypes from FUT1 and FUT2.

    This helper function extracts the alphanumeric phenotype values from the FUT1 and
    FUT2 BloodGroup objects, then returns the Cartesian product of these phenotype sets
    as a list of tuples.

    Args:
        res (dict[str, BloodGroup]): Dictionary of blood group results containing at
            least keys "FUT1" and "FUT2". Each corresponding value is a BloodGroup object
            with a 'phenotypes' attribute, where phenotypes are categorized by PhenoType.

    Returns:
        list[tuple[str, str]]: A list of tuples, where each tuple represents a pair of
            phenotype values, one from FUT1 and one from FUT2.
    """
    fut1 = res["FUT1"]
    fut2 = res["FUT2"]
    pheno_key: PhenoType = PhenoType.alphanumeric
    fut1s = set(fut1.phenotypes[pheno_key].values())
    fut2s = set(fut2.phenotypes[pheno_key].values())
    fut1_and_2 = list(product(fut1s, fut2s))
    return fut1_and_2


def FUT1(res: dict[str, BloodGroup]) -> dict[str, BloodGroup]:
    """Update FUT1 phenotypes based on interactions between FUT1 and FUT2.

    In the context of fucosyltransferase blood group antigens, FUT1 (H), FUT2 (Se),
    and FUT3 (Le) have interrelated phenotypes. This function updates the alphanumeric
    phenotypes of the FUT1 BloodGroup by combining information from FUT1 and FUT2.
    For each phenotype in FUT1, it finds all phenotype pairs from FUT1 and FUT2 that
    include the original FUT1 phenotype, concatenates them, and then stores the result
    as the new phenotype value.

    Args:
        res (dict[str, BloodGroup]): Dictionary of blood group results, where keys
            are 'FUT1', 'FUT2', and 'FUT3', each mapped to their respective BloodGroup
            data.

    Returns:
        dict[str, BloodGroup]: The updated dictionary with modified FUT1 phenotypes.
    """

    new_phenos = {}
    fut1_and_2 = fut_helper(res)
    for pair, pheno in res["FUT1"].phenotypes[PhenoType.alphanumeric].items():
        new = []
        for combo in fut1_and_2:
            if pheno in combo:
                new.append(",".join(combo))
        new_phenos[pair] = "/".join(sorted(new))

    res["FUT1"].phenotypes[PhenoType.alphanumeric] = new_phenos

    return res


def FUT3(res: dict[str, BloodGroup]) -> dict[str, BloodGroup]:
    """
    H = FUT1
    Se = FUT2
    Le = FUT3

    In individuals with an active FUT2 (Secretor or SE) gene, which encodes a fully
    active α(1,2)-fucosyltransferase (see H Blood Group System), predominantly
    Leb (and related Lewis antigens depending on ABO group, e.g. ALeb in group
    A) is made alongside trace amounts of Lea. The trace amounts of Lea produced
    are typically undetectable using serological methods, and usually, a Le(a–b+)
    phenotype is reported.

    functional == active???!!!

    LE	When FUT2 is H+ Se+ and FUT3 is functional = Le(a-b+)
        When FUT2 is H+w Se+ and FUT3 is functional =Le(a+b+)
        When FUT2 is homozygous null (Se-) and FUT3 is functional = Le(a+b-)
        When FUT3 is homozygous null, regardless of FUT2 functionality = Le(a-b-)
    H	FUT1: RBC expresssion
        FUT2: Secretor status - Report phenotype as Se+/-

    Args:
        res (dict): Dictionary of blood group results, where keys are 'FUT1', 'FUT2',
        and 'FUT3', each with their respective `BloodGroup` data.

    Returns:
        dict: Updated dictionary with modified FUT3 phenotypes based on FUT1 and FUT2
        interactions.
    """

    fut1_and_2 = fut_helper(res)

    new_phenos = {}
    for pair, pheno in res["FUT3"].phenotypes[PhenoType.alphanumeric].items():
        new_pheno = []
        if "active" in pheno.lower():
            for combo in fut1_and_2:
                if "Se-" in combo:
                    new_pheno.append("Le(a+b-)")
                elif "H+w" in combo and "Se+" in combo:
                    new_pheno.append("Le(a+b+)")
                elif "H+" in combo and any(se in combo for se in ["Se+", "Se+w"]):
                    new_pheno.append("Le(a-b+)")
                else:
                    raise BeyondLogicError(
                        message="Unexpected FUT3 value.",
                        context=f"Received value: {combo}",
                    )
        elif "Le(a-b-)" in pheno:
            new_pheno.append("Le(a-b-)")
        else:
            raise BeyondLogicError(
                message="Unexpected FUT3 pheno.", context=f"Received value: {pheno}"
            )

        new_phenos[pair] = "/".join(sorted(set(new_pheno)))

    res["FUT3"].phenotypes[PhenoType.alphanumeric] = new_phenos

    return res


def choose_class_type(bg_type, ant_type):
    """Choose the antigen class based on blood group and phenotype type.

    Args:
        bg_type (BgName): The blood group type.
        ant_type (PhenoType): The phenotype type, either alphanumeric or numeric.

    Returns:
        type: The antigen class corresponding to the given types.
    """
    class_types = {
        PhenoType.alphanumeric: {
            BgName.XG: an.AlphaNumericAntigenXG,
            BgName.VEL: an.AlphaNumericAntigenVel,
            BgName.GYPA: an.AlphaNumericAntigenMNS,
            BgName.GYPB: an.AlphaNumericAntigenMNS,
            BgName.ABO: an.AlphaNumericAntigenABO,
            BgName.RHCE: an.AlphaNumericAntigenRHCE,
            BgName.RHD: an.AlphaNumericAntigenRHD,
        },
        PhenoType.numeric: {
            BgName.RHCE: an.NumericAntigenRHCE,
            BgName.VEL: an.NumericAntigenVel,
        },
    }
    return (
        class_types[ant_type].get(bg_type, an.NumericAntigen)
        if ant_type == PhenoType.numeric
        else class_types[ant_type].get(bg_type, an.AlphaNumericAntigen)
    )


def make_values_dict(
    values_strs: list[str], ant_type: PhenoType, bg_type: BgName
) -> dict[str, list[an.Antigen]] | None:
    """Construct a dictionary of Antigen objects from antigen value strings.

    This function processes a list of antigen strings and creates a dictionary where the keys
    are antigen base names and the values are lists of corresponding Antigen objects. The antigen
    objects are instantiated using a class determined by the phenotype type (ant_type) and blood
    group type (bg_type).

    Processing steps:
      1. If any string in values_strs contains a period ('.'), return None. This is because some
         blood groups may not have consistent numeric or alphanumeric allele representations.
      2. Determine the antigen class to use by calling `choose_class_type(bg_type, ant_type)`.
      3. Modify the input strings:
           - For numeric antigens, split each string by ":" and use the second part.
           - For alphanumeric antigens, use the string as is.
      4. Concatenate the modified strings, split by commas, and count the occurrences using a Counter.
      5. For each unique antigen string and its count:
           - Create an antigen object with:
               - given_name set to the antigen string.
               - expressed flag determined by the absence of '-' in the antigen string.
               - homozygous flag set to True if the count equals 2.
               - antithetical_relationships set based on the global ANTITHETICAL dictionary for the given
                 phenotype and blood group types.
           - Assert that the antigen's base_name is not empty.
           - Append the antigen to the dictionary under its base_name.

    Args:
        values_strs (list[str]): A list of strings representing antigen values.
        ant_type (PhenoType): The phenotype type (numeric or alphanumeric) of the antigen.
        bg_type (BgName): The blood group type used to determine the antigen class and relationships.

    Returns:
        dict[str, list[an.Antigen]] | None: A dictionary mapping antigen base names to lists of Antigen objects,
        or None if the input format is invalid.
    """
    antigens = defaultdict(list)
    # all dicts are ordered
    if "." in values_strs:
        return None  # some BGs just don't have numeric for all or alpha for all alleles
    ant_class = choose_class_type(bg_type, ant_type)
    values_strs_mod = []
    for values_str in values_strs:
        values_strs_mod.append(
            values_str.split(":")[1] if ant_type == PhenoType.numeric else values_str
        )

    # NB order only matters for ref, where there's just one str (Counter is ordered)
    for ant_str, count in Counter(",".join(values_strs_mod).split(",")).items():
        ant = ant_class(
            given_name=ant_str,
            expressed="-" not in ant_str,
            homozygous=count == 2,
            antithetical_relationships=ANTITHETICAL[ant_type].get(bg_type, {}),
        )
        assert ant.base_name != ""
        antigens[ant.base_name].append(ant)

    return antigens


@apply_to_dict_values
def add_ref_phenos(bg: BloodGroup, df: pd.DataFrame) -> BloodGroup:
    """Add reference phenotypes to the BloodGroup from a DataFrame.

    This function filters the DataFrame for the reference genotype corresponding
    to the BloodGroup type, then constructs alphanumeric and numeric phenotype
    dictionaries using make_values_dict. These are stored in the bg.misc attribute.

    Args:
        bg (BloodGroup): The BloodGroup object to update.
        df (pd.DataFrame): A DataFrame containing phenotype data. Must include
            'type', 'Reference_genotype', 'Phenotype_alt', and 'Phenotype' columns.

    Returns:
        BloodGroup: The updated BloodGroup with reference phenotypes added.

    Raises:
        AssertionError: If the filtered DataFrame does not have exactly one row.
    """
    df["type"] = df["type"].astype("category")
    df_ref = df.loc[(df["Reference_genotype"] == "Yes") & (df["type"] == bg.type)]

    try:
        assert df_ref.shape == (1, 19)
    except AssertionError:
        ic(df_ref.shape, df_ref)
    bg.misc = {}

    bg.misc["ref_PhenoType.alphanumeric"] = make_values_dict(
        [df_ref.iloc[0]["Phenotype_alt"]],
        PhenoType.alphanumeric,
        BgName.from_string(bg.type),
    )
    bg.misc["ref_PhenoType.numeric"] = make_values_dict(
        [df_ref.iloc[0]["Phenotype"]], PhenoType.numeric, BgName.from_string(bg.type)
    )

    return bg


@apply_to_dict_values
def instantiate_antigens(bg: BloodGroup, ant_type: PhenoType) -> BloodGroup:
    """Instantiate antigen objects for a BloodGroup based on allele phenotype changes.

    This function iterates over allele pairs in the BloodGroup (using either the CO or
    NORMAL state, as determined by the inner helper) and generates a dictionary of antigens.
    For each allele pair, phenotype changes are extracted based on the phenotype type (numeric
    or alphanumeric). If the phenotype changes indicate that the pair is informative (i.e. not
    just missing data represented by '.'), the function uses a pre-filled partial of
    make_values_dict to generate antigen objects. The resulting dictionary is stored in the
    bg.misc attribute under a key formatted as "antigens_{ant_type.name}".

    Args:
        bg (BloodGroup): A BloodGroup object containing allele pairs and associated data.
        ant_type (PhenoType): The phenotype type, which influences the phenotype attribute
            to be used and the antigen class to instantiate.

    Returns:
        BloodGroup: The updated BloodGroup with instantiated antigens added to bg.misc.
    """

    def select_allele_state(bg):
        """Select the allele state to use for antigen instantiation.

        For blood groups of type "KN", if co-existing alleles (CO) are present, then
        AlleleState.CO is returned. Otherwise, AlleleState.NORMAL is used.

        Args:
            bg (BloodGroup): The BloodGroup object.

        Returns:
            AlleleState: The allele state (CO or NORMAL) to be used.
        """
        if AlleleState.CO in bg.alleles and bg.alleles[AlleleState.CO] is not None:
            assert bg.type == "KN"
            return AlleleState.CO
        return AlleleState.NORMAL

    def get_pheno_changes(current_pair: Pair, ant_type: PhenoType) -> list[str]:
        """Extract phenotype changes from an allele pair based on phenotype type.

        For numeric phenotypes, the 'phenotype' attribute is used; for alphanumeric,
        'phenotype_alt' is used. When the BloodGroup type is "ABO", the function uses
        alleles with expressed phenotypes; otherwise, it uses all alleles in the pair.

        Args:
            current_pair (Pair): The allele pair from which to extract phenotype changes.
            ant_type (PhenoType): The phenotype type (numeric or alphanumeric).

        Returns:
            list[str]: A list of phenotype strings derived from the allele pair.
        """

        def filter_abo_types(alleles: list[Allele]) -> frozenset[str]:
            """
            Filters ABO types based on the presence of 'ABO*O.'.

            Args:
                Allele1: The first ABO type string.
                Allele2: The second ABO type string.

            Returns:
                A frozenset containing:
                - Both strs if both contain 'ABO*O.'.
                - Only the str that does NOT contain 'ABO*O.' if one of them does.
                - Both strs if neither contains 'ABO*O.'.
            """
            if len(alleles) == 1:
                return alleles
            elif len(alleles) == 2:
                a1, a2 = alleles
                contains_o1 = "ABO*O." in a1.genotype
                contains_o2 = "ABO*O." in a2.genotype

                if contains_o1 and contains_o2:
                    return frozenset([a1, a2])
                elif contains_o1 and not contains_o2:
                    return frozenset([a2])
                elif not contains_o1 and contains_o2:
                    return frozenset([a1])
                else:  # Neither contains 'ABO*O.'
                    return frozenset([a1, a2])
            else:
                raise ValueError("ABO allele count wrong")

        phenotype_attr = (
            "phenotype" if ant_type == PhenoType.numeric else "phenotype_alt"
        )
        alleles_to_use = (
            filter_abo_types(current_pair.alleles)
            if bg.type == "ABO"
            else current_pair.alleles
        )

        return [getattr(allele, phenotype_attr) for allele in alleles_to_use]

    pair_antigens = {}
    make_values_dict_pre_filled = partial(
        make_values_dict, ant_type=ant_type, bg_type=BgName.from_string(bg.type)
    )
    for pair in bg.alleles[select_allele_state(bg)]:
        pheno_changes = get_pheno_changes(pair, ant_type)
        if set(pheno_changes) != {"."}:  # 1 or more defininitions missing
            if pair.allele1 == pair.allele2:
                pheno_changes *= 2
            pair_antigens[pair] = make_values_dict_pre_filled(values_strs=pheno_changes)
            # TODO need to do this specifically for the PhenoType and only if the other one isn't {"."}
            # TODO replace bg.type globally

    bg.misc[f"antigens_{ant_type.name}"] = pair_antigens

    return bg


@apply_to_dict_values
def get_phenotypes1(bg: BloodGroup, ant_type: PhenoType) -> BloodGroup:
    """Integrate antigen information with reference phenotypes for a BloodGroup.

    This function retrieves the reference antigen data from bg.misc using the key
    corresponding to the provided phenotype type. It then iterates over the antigen
    pairs stored in bg.misc for that phenotype type and merges each pair with the
    reference data when necessary. The merging behavior depends on whether the allele
    is homozygous, heterozygous, or has antithetical antigens. The combined antigen
    information is stored back in bg.misc under a key formed by concatenating
    "antigens_and_ref_" with the antigen type name.

    Args:
        bg (BloodGroup): The BloodGroup object containing antigen and reference data.
        ant_type (PhenoType): The phenotype type (numeric or alphanumeric) for which
            antigen data is being processed.

    Returns:
        BloodGroup: The updated BloodGroup with combined antigen and reference
            information in its misc attribute.

    Raises:
        BeyondLogicError: If an unexpected configuration of allele antigens is encountered.
    """
    reference = bg.misc[f"ref_{str(ant_type)}"]
    if reference is None:
        return bg
    else:
        ref = copy.deepcopy(reference)

    d = {}
    for pair, antigens in bg.misc[f"antigens_{ant_type.name}"].items():
        if antigens is None:
            continue
        antigens_with_ref_if_needed = {}
        for ant_pos, allele_antigens in antigens.items():
            if (
                bg.type in ["ABO"]
                or allele_antigens[0].homozygous
                and len(allele_antigens) == 1
            ):
                antigens_with_ref_if_needed[ant_pos] = allele_antigens
            elif len(allele_antigens) == 2:
                # ic(allele_antigens)
                assert all(not allele.homozygous for allele in allele_antigens), (
                    "Expected both alleles to be heterozygous"
                )
                antigens_with_ref_if_needed[ant_pos] = allele_antigens
            elif len(allele_antigens) == 1:
                ant = allele_antigens[0]
                assert not ant.homozygous
                add_ref = True
                if ant.antithetical_antigen is not None:
                    for anti_ant in ant.antithetical_antigen:
                        if anti_ant.base_name in antigens.keys():
                            add_ref = False
                if add_ref:
                    ant_pair = allele_antigens + ref[ant_pos]
                    antigens_with_ref_if_needed[ant_pos] = ant_pair
                else:
                    antigens_with_ref_if_needed[ant_pos] = allele_antigens
            else:
                raise BeyondLogicError(
                    message="Unexpected allele_antigens.",
                    context=f"Received value: {allele_antigens}",
                )
            d[pair] = antigens_with_ref_if_needed

    bg.misc[f"antigens_and_ref_{ant_type.name}"] = d

    return bg


@apply_to_dict_values
def get_phenotypes2(bg: BloodGroup, ant_type: PhenoType) -> BloodGroup:
    """Merge antigen phenotype options and update BloodGroup phenotype mapping.

    This function integrates antigen data with reference information by merging
    phenotype options for each allele pair. It iterates over the antigen pairs stored
    in bg.misc under the key "antigens_and_ref_{ant_type.name}" and applies merging rules
    based on antigen equality, given names, and weight comparisons. The resulting merged
    antigen data is then stored in bg.phenotypes for the provided phenotype type.

    Args:
        bg (BloodGroup): A BloodGroup object containing antigen, reference, and phenotype data.
        ant_type (PhenoType): The phenotype type (numeric or alphanumeric) to process.

    Returns:
        BloodGroup: The updated BloodGroup with merged phenotype information.

    Raises:
        BeyondLogicError: If an unexpected antigen pattern is encountered during merging.
    """
    reference = bg.misc[f"ref_{str(ant_type)}"]
    if reference is None:
        return bg

    for pair, antigens_with_ref_if_needed in bg.misc[
        f"antigens_and_ref_{ant_type.name}"
    ].items():
        if antigens_with_ref_if_needed is None:
            continue

        merged_pheno2 = []
        for _, options in antigens_with_ref_if_needed.items():
            if len(options) == 2:
                ant1, ant2 = options
                assert ant1.base_name == ant2.base_name
                if ant1 == ant2 and ant1.given_name != ant2.given_name:
                    merged_pheno2.append(ant1)
                    merged_pheno2.append(ant2)
                elif ant1 == ant2 and ant1.given_name == ant2.given_name:
                    if ant1.expressed and ant2.expressed:
                        ant1.homozygous = (
                            True  # TODO fix - can freeze class withotut this
                        )
                    merged_pheno2.append(ant1)
                elif ant1 > ant2:
                    if ant1.expressed and ant2.expressed:
                        ant2.homozygous = True
                    merged_pheno2.append(ant2)
                elif ant2 > ant1:
                    if ant1.expressed and ant2.expressed:
                        ant1.homozygous = True
                    merged_pheno2.append(ant1)
                else:
                    raise BeyondLogicError(
                        message="Unexpected antigne pattern.",
                        context=f"Received value: {options}",
                    )
            elif len(options) == 1:
                merged_pheno2.append(options[0])
            else:
                raise BeyondLogicError(
                    message="Unexpected antigne pattern.",
                    context=f"Received value: {options}",
                )

        bg.phenotypes[ant_type][pair] = merged_pheno2

    return bg


def count_expressed_ants(current_ant, current_base_names):
    """Count the number of expressed antigens including antithetical ones.

    The count is determined by:
      - Adding 1 if the current antigen is expressed.
      - Multiplying by 2 if the current antigen is homozygous.
      - For each antithetical antigen in current_ant.antithetical_antigen,
        if it exists in current_base_names and is expressed, add 1 if heterozygous,
        or 2 if homozygous.

    Args:
        current_ant (an.Antigen): The antigen to evaluate.
        current_base_names (dict[str, an.Antigen]): A mapping from antigen base names to
            corresponding antigen objects.

    Returns:
        int: The total number of expressed antigens.
    """
    number_of_expressed_ants = 1 if current_ant.expressed else 0
    if current_ant.homozygous:
        number_of_expressed_ants *= 2

    for antithetical_ant in current_ant.antithetical_antigen:
        ant_to_count = current_base_names.get(antithetical_ant.base_name)
        if ant_to_count is None:
            continue
        if ant_to_count.expressed:
            if ant_to_count.homozygous:
                number_of_expressed_ants += 2
            else:
                number_of_expressed_ants += 1

    return number_of_expressed_ants


@apply_to_dict_values
def internal_anithetical_consistency_HET(
    bg: BloodGroup, ant_type: PhenoType
) -> BloodGroup:
    """
    Ensure correct antigen expression for heterozygous antithetical pairs.

    In certain cases, when an antithetical pair is present and only one is
    marked with '-', the code sets it to positive if heterozygous. This often
    arises in tricky alleles like KEL*02.03 (which flips several antigens).

    Steps:
        1. count expressed antithetical antigens (from antigens not
        ant.antithetical_antigen)
        2. If two are already expressed, do nothing.
        3. If not, retrieve reference states to correct expression for the
           missing antigens.

    *NB some antigens have 3 antitheticals described ie (3,4, and 21 are
    anithetical in KEL, this is well tested). Others seem to have 3 but only
    2 are known - listed in sudo nulls (less well tested)

    Args:
        bg (BloodGroup): The blood group data, containing phenotypes and references.
        ant_type (PhenoType): The phenotype group to process.

    Returns:
        A potentially updated BloodGroup with corrected antithetical expressions.

    """

    def check_expression(antigen):
        ref_ant = reference[antigen.base_name][0]
        if ref_ant.expressed:
            new_antigens.append(ref_ant)
            already_checked.add(ref_ant.base_name)

    new_phenos = []
    ref = bg.misc[f"ref_{str(ant_type)}"]
    if ref is not None:
        reference = copy.deepcopy(ref)
    else:
        return bg

    for pair, antigens in bg.phenotypes[ant_type].items():
        null = pair.allele1.null or pair.allele2.null
        already_checked = set()
        new_antigens = []
        if pair.allele1.phenotype == ".":
            new_phenos.append((pair, antigens))
            continue
        if (
            "N." in pair.allele1.genotype.upper()
            or "N." in pair.allele2.genotype.upper()
        ):
            # TODO think about how to handle nulls better
            new_phenos.append((pair, antigens))
            continue
        base_names = {ant2.base_name: ant2 for ant2 in antigens}

        for ant in antigens:
            if ant.base_name in already_checked:
                continue
            if ant.antithetical_antigen:
                no_expressed = count_expressed_ants(ant, base_names)
                # no change
                if no_expressed == 2:
                    for ant2 in antigens:
                        if ant2.base_name == ant.base_name:
                            new_antigens.append(ant2)
                            already_checked.add(ant2.base_name)
                    continue
                # add expressed ants
                for ant2 in antigens:
                    if ant2.base_name == ant.base_name and ant2.expressed:
                        new_antigens.append(ant2)
                        already_checked.add(ant2.base_name)
                # add expressed ref ants
                if ant.base_name not in already_checked:
                    check_expression(ant)
                for antithetical_ant in ant.antithetical_antigen:
                    if antithetical_ant.base_name in already_checked:
                        continue
                    check_expression(antithetical_ant)
            else:
                new_antigens.append(ant)
            already_checked.add(ant.base_name)

        base_names_new = {ant2.base_name: ant2 for ant2 in new_antigens}
        for ant_base_name, ant in base_names.items():
            if ant_base_name not in base_names_new:
                assert not ant.expressed
                new_antigens.append(ant)
        new_phenos.append((pair, new_antigens))
        base_names_new = {ant2.base_name: ant2 for ant2 in new_antigens}

        for ant in base_names:
            assert ant in base_names_new

        sudo_nuls = [
            "CO*01.-04",
            "ER*01.-03",
            "CO*M.01",
            "KEL*02.-14.1",
            "DI*02.15",
            "DI*02.16",
            "DI*02.11",
            "DI*02.12",
            "DI*02.17",
            "DI*02.18",
            "DI*02.22",
            "DI*02.09",
            "KEL*02M.05",
            "GYPB*21",  # these GYPBs can be S or s
            "GYPB*23",
            "GYPB*24",
            "RHCE*01.34",
            "RHCE*01.35",
            "GYPA*11",
        ]
        # Di11/12 and 15/16 17/18 antithetical and same in ref so, yes, sudo null
        # (or more accurately, ref is third [unamed] ant)
        # DI*02.22/DI*02.09 (ant 9/22) seem to have a third antigen as well,
        #  ~300 examples in UKB - not a sudo null, but here till I find a better place
        # KEL*02M.05', is just KO, no other info....
        for ant in new_antigens:
            if (
                ant.antithetical_antigen
                and pair.allele1.genotype not in sudo_nuls
                and pair.allele2.genotype not in sudo_nuls
            ):  # TODO move to dif func
                final_no_expressed = count_expressed_ants(ant, base_names_new)
                to_neg = "to_neg" in str(pair) and "RHCE" in str(pair)
                if not to_neg:
                    try:
                        assert final_no_expressed == 2
                    except AssertionError:
                        ic(
                            "Expressed antigens != 2! plz report to devs",
                            bg.sample,
                            new_antigens,
                            pair.allele1,
                            pair.allele2,
                            str(pair),
                            ant,
                            final_no_expressed,
                            null,
                        )

    for pair, merged_pheno in new_phenos:
        bg.phenotypes[ant_type][pair] = merged_pheno

    return bg


def make_ant(antithetical_ant, ant_type, bg_type, expressed):
    """Create an antigen instance from an antithetical antigen.

    The antigen class is chosen based on the blood group type and phenotype type.
    This function uses the base name of the provided antithetical antigen and sets
    the antigen as homozygous with the specified expressed state. Antithetical
    relationships are assigned from the global ANTITHETICAL dictionary.

    Args:
        antithetical_ant (an.Antigen): An antigen whose base name is used.
        ant_type (PhenoType): The phenotype type (numeric or alphanumeric) of the antigen.
        bg_type (BgName): The blood group type to determine the antigen class.
        expressed (bool): Flag indicating whether the antigen is expressed.

    Returns:
        an.Antigen: A newly instantiated antigen object.
    """
    ant_class = choose_class_type(bg_type, ant_type)
    return ant_class(
        given_name=antithetical_ant.base_name,
        expressed=expressed,  # False if ant.expressed else True,
        homozygous=True,
        antithetical_relationships=ANTITHETICAL[ant_type].get(bg_type, {}),
    )


@apply_to_dict_values
def internal_anithetical_consistency_HOM(
    bg: BloodGroup, ant_type: PhenoType
) -> BloodGroup:
    """Ensure internal antithetical consistency for homozygous antigen pairs.

    This function enforces consistency among antigens for a BloodGroup by
    supplementing the existing phenotype mappings with antithetical antigen objects
    when necessary. For each allele pair in bg.phenotypes corresponding to the given
    phenotype type, if the allele pair is informative (i.e. its first allele's phenotype
    is not missing) and an antigen is homozygous with defined antithetical relationships,
    the function checks whether each antithetical antigen (based on its base name) is
    already present. If not, a new antigen is created using the `make_ant` helper function,
    with its expressed state set inversely to the original antigen's expressed flag.

    Example scenario:
        - Given a sample where KN*01.07/KN*01.07 is observed and KN*01.07 has a defined
          set of antigens such as KN:7, -9, 10.
        - The expected output may include a full set like:
          KN:1, -2, 3, -6, -4, 7, 5, 8, -9, 10, 11, -12, 13.
        - Explanation:
            * If KN7 is homozygous, then its antithetical antigen KN4 should not be present
              as a positive value, enforcing the relationship that 4 is antithetical to 7.
            * Similarly, if there is no supporting evidence for KN6 on the positive side,
              it remains negative.

    Args:
        bg (BloodGroup): The BloodGroup object containing phenotype data.
        ant_type (PhenoType): The phenotype type (e.g., numeric or alphanumeric) for which
            the antigen consistency should be enforced.

    Returns:
        BloodGroup: The updated BloodGroup with merged antigen lists that respect antithetical
            consistency for homozygous alleles.
    """

    new_phenos = []
    for pair, antigens in bg.phenotypes[ant_type].items():
        new_antigens = []
        if pair.allele1.phenotype == ".":
            continue
        base_names = [ant2.base_name for ant2 in antigens]
        base_names_dict = {ant2.base_name: ant2 for ant2 in antigens}

        for ant in antigens:
            if ant.antithetical_antigen and ant.homozygous:
                if len(ant.antithetical_antigen) > 1:
                    no_expressed = count_expressed_ants(ant, base_names_dict)
                    if pair.allele1.null and pair.allele2.null:
                        assert no_expressed == 0
                    elif pair.allele1.null or pair.allele2.null:
                        assert no_expressed == 1
                    else:
                        assert no_expressed == 2
                for antithetical_ant in ant.antithetical_antigen:
                    if antithetical_ant.base_name not in base_names:
                        bg_type = BgName.from_string(bg.type)
                        # TODO freeze ants again??? - should never set to expressed?
                        new_antigens.append(
                            make_ant(
                                antithetical_ant,
                                ant_type,
                                bg_type,
                                False if ant.expressed else True,
                            )
                        )

        new_phenos.append((pair, new_antigens + antigens))
    for pair, merged_pheno in new_phenos:
        bg.phenotypes[ant_type][pair] = merged_pheno

    return bg


@apply_to_dict_values
def include_first_antithetical_pair(bg: BloodGroup, ant_type: PhenoType) -> BloodGroup:
    """Include the first antithetical pair from a reference into the blood group's
    phenotype.

    If the blood group is one of the FUT or ABO types, the function returns immediately
    with no changes. It then checks how many primary antitheticals are required for
    this blood group type. If zero, it returns immediately again. If no reference is
    found, it also returns immediately. Otherwise, it appends the required number of
    reference antigens to each pair in the phenotype if they are missing.

    Args:
        bg (BloodGroup): The blood group to update.
        ant_type (PhenoType): The domain of phenotypes to modify.

    Returns:
        BloodGroup: The updated (or unmodified) blood group.
    """
    if bg.type in ["FUT1", "FUT2", "FUT3", "ABO", "RHCE", "RHD"]:
        return bg
    no_of_positions_required = number_of_primary_antitheticals.get(bg.type, 2)
    if no_of_positions_required == 0:
        return bg
    new_phenos = []
    ref = bg.misc[f"ref_{str(ant_type)}"]
    if ref is not None:
        # reference = ref.copy()
        reference = copy.deepcopy(ref)
    else:
        return bg
    for pair, merged_pheno in bg.phenotypes[ant_type].items():
        for i, name_ant in enumerate(reference.items(), start=1):
            _, ref_ant_list = name_ant
            # ic(pair, ref_ant_list)
            assert len(ref_ant_list) == 1
            ref_ant = ref_ant_list[0]

            if ref_ant.base_name not in [
                allele_ant.base_name for allele_ant in merged_pheno
            ]:
                merged_pheno.append(ref_ant)
            if i == no_of_positions_required:
                break
        new_phenos.append((pair, merged_pheno))
    for pair, merged_pheno in new_phenos:
        bg.phenotypes[ant_type][pair] = merged_pheno

    return bg


@apply_to_dict_values
def sort_antigens(bg: BloodGroup, ant_type: PhenoType) -> BloodGroup:
    """Sort antigens for each allele pair in a BloodGroup based on a reference order.

    This function retrieves reference antigen data from bg.misc corresponding to the
    provided phenotype type. It then creates an ordering map from the reference keys
    (antigen base names) and sorts the merged antigen lists (bg.phenotypes) for each allele
    pair according to this order. Antigens with base names not found in the reference
    are placed at the end.

    Args:
        bg (BloodGroup): A BloodGroup object containing phenotype and miscellaneous data.
        ant_type (PhenoType): The phenotype type (numeric or alphanumeric) for which antigens
            are being sorted.

    Returns:
        BloodGroup: The updated BloodGroup with sorted antigen lists in bg.phenotypes.
    """
    new_phenos = []
    reference = bg.misc[f"ref_{str(ant_type)}"]
    if reference is None:
        return bg
    else:
        ref = copy.deepcopy(reference)
    order_map = {name: position for position, name in list(enumerate(ref.keys()))}
    for pair, merged_pheno in bg.phenotypes[ant_type].items():
        sorted_merged_pheno = sorted(
            merged_pheno,
            key=lambda antigen: order_map.get(antigen.base_name, float("inf")),
        )
        new_phenos.append((pair, sorted_merged_pheno))
    for pair, sorted_merged_pheno in new_phenos:
        bg.phenotypes[ant_type][pair] = sorted_merged_pheno
    # if bg.type == "RHCE":
    #     ic(2222222, bg.phenotypes[ant_type])
    return bg


@apply_to_dict_values
def phenos_to_str(bg: BloodGroup, ant_type: PhenoType) -> BloodGroup:
    """Convert a list of antigen objects into a consolidated phenotype string.

    If the blood group has a RAW allele, its phenotype (before the first colon)
    is used as a prefix. If there's no RAW allele, the POS allele's prefix is used.
    Then, for each phenotype pair in ``bg.phenotypes[ant_type]``, a comma-separated
    string of antigen names is created. If ``bg.type`` is 'ABO', the list of antigens
    is sorted before joining. Otherwise, the order remains as is.

    Args:
        bg (BloodGroup): The blood group object containing alleles and phenotypes.
        ant_type (PhenoType): The phenotype category (e.g., alphanumeric or numeric).

    Returns:
        BloodGroup: The updated blood group with phenotype strings in
        ``bg.phenotypes[ant_type]`` instead of lists of antigen objects.

    Raises:
        IndexError: If neither RAW nor POS allele lists contain any alleles.
    """

    allele_name = bg.alleles[AlleleState.RAW][0].phenotype.split(":")[0]

    for pair, merged_pheno in bg.phenotypes[ant_type].items():
        ants = [ant.name for ant in merged_pheno]
        as_str = ",".join(sorted(ants)) if bg.type == "ABO" else ",".join(ants)
        # if bg.type == "RHCE":
        #     ic(111111222222, as_str)
        pheno = (
            as_str if ant_type == PhenoType.alphanumeric else f"{allele_name}:{as_str}"
        )
        # if bg.type == "RHCE":
        #     ic(11111133333, pheno)

        bg.phenotypes[ant_type][pair] = pheno
    # if bg.type == "RHCE":
    #     ic(111111, bg.phenotypes[ant_type])
    return bg


@apply_to_dict_values
def combine_anitheticals(bg: BloodGroup) -> BloodGroup:
    def split_list(input_list: list[str]) -> tuple[list[str], list[str]]:
        """
        Splits the input list into two lists based on specific criteria.

        Args:
            input_list (List[str]): The list of strings to be split.

        Returns:
            Tuple[List[str], List[str]]: A tuple containing two lists:
                - The first list contains strings that have a '(' and there are more
                than one string starting with the same prefix (substring before '(').
                - The second list contains all other strings.
        Example:

        Ie
        Example1 = ['Kn(a+)','Kn(b-)','McC(a+)','Sl1+', 'Yk(a+)','McC(b-)', 'Vil-',
        'Sl3+','KCAM-','KDAS+','DACY-','YCAD+','KNMB+']
        Or
        Example2 = ['AUG1+', 'At(a+)', 'ATML-', 'ATAM+']

        And splits the list into 2 based on if the strs have a '(' and if there are
        more than 1 of the same type ('Kn(a+)' and 'Kn(b-)' are the same type because
        they both start with 'kn')

        The output from Example1 would be
        List_with_parens = ['Kn(a+)','Kn(b-)','McC(a+)','McC(b-)']
        Rest = ['Sl1+', 'Yk(a+)', 'Vil-','Sl3+','KCAM-','KDAS+','DACY-','YCAD+','KNMB+']

        Example2
        List_with_parens = []
        Rest = ['AUG1+', 'At(a+)', 'ATML-', 'ATAM+']
        """
        parens: list[str] = []
        rest: list[str] = []
        strings_with_parens: list[str] = [s for s in input_list if "(" in s]
        prefixes: list[str] = [s.split("(")[0] for s in strings_with_parens]
        prefix_counts: Counter = Counter(prefixes)
        for s in input_list:
            if "(" in s:
                prefix = s.split("(")[0]
                if prefix_counts[prefix] > 1:
                    parens.append(s)
                else:
                    rest.append(s)
            else:
                rest.append(s)

        return parens, rest

    to_update = []
    for pair, pheno in bg.phenotypes[PhenoType.alphanumeric].items():
        antigens = pheno.split(",")
        assert len(antigens) > 0

        has_paren, no_paren = split_list(antigens)
        if has_paren and no_paren:
            combined = ",".join([combine_expressions(has_paren)] + no_paren)
        elif has_paren:
            combined = ",".join([combine_expressions(has_paren)])
        elif no_paren:
            combined = ",".join(no_paren)

        to_update.append((pair, combined))

    for pair, combined in to_update:
        bg.phenotypes[PhenoType.alphanumeric][pair] = combined
    return bg


def combine_expressions(components: list[str]) -> str:
    """
    Combines multiple expressions with the same prefix into a single expression.

    This function takes a string containing multiple expressions separated by commas,
    where each expression has the same prefix (e.g., 'In', 'LW', 'Tc') followed by
    content in parentheses. It combines these expressions into a single expression
    with the same prefix and concatenated content.

    Args:
        expression (str): A string containing comma-separated expressions to be
        combined.

    Returns:
        str: A single combined expression.

    Examples:
        ['In(a-)','In(b+)'] > 'In(a-b+)'
        ['LW(a+)','LW(b-)'] > 'LW(a+b-)'
        ['Tc(a+)','Tc(b-)','Tc(c-)'] > 'Tc(a+b-c-)'
        ['Au(a-)', 'Lu(a-)', 'Au(b+)', 'Lu(b+)'] > 'Au(a-b+),Lu(a-b+)'

    """

    def extract_content(comp: str) -> str:
        """
        Extracts the content between parentheses from a component.

        Args:
            comp (str): A component string containing content in parentheses.

        Returns:
            str: The content between the parentheses.
        """
        start: int = comp.index("(") + 1
        end: int = comp.rindex(")")
        return comp[start:end]

    clustered = defaultdict(list)
    for phenotype in components:
        prefix_end: int = phenotype.index("(")
        prefix: str = phenotype[:prefix_end]
        clustered[prefix].append(phenotype)
    pheno = []
    for prefix, sub_components in clustered.items():
        combined_content: str = "".join(
            sorted(extract_content(comp) for comp in sub_components)
        )
        pheno.append(f"{prefix}({combined_content})")
    return ",".join(pheno)


def null_or_mod(pair: Pair, check: str) -> bool:
    """Check if the given pattern is present in both allele genotypes of a pair.

    This function constructs a pattern by appending a period ('.') to the provided
    check string, converts both allele genotypes to uppercase, and verifies whether the
    pattern exists in both.

    Args:
        pair (Pair): A Pair object containing two allele objects.
        check (str): A substring used to build the pattern to search for.

    Returns:
        bool: True if the pattern is present in both allele genotypes, otherwise False.
    """

    pattern = f"{check}."  # TODO AUG*01N not covered!
    return (
        pattern in pair.allele1.genotype.upper()
        and pattern in pair.allele2.genotype.upper()
    )


@apply_to_dict_values
def modify_FY(bg: BloodGroup, ant_type: PhenoType) -> BloodGroup:
    """Modify FY BloodGroup phenotypes to indicate erythroid cell restriction.

    For BloodGroup of type "FY", this function appends a suffix to the phenotype string
    if either allele's genotype is one of the specified erythroid null alleles and if both
    alleles contain the 'N.' pattern. This indicates that the phenotype is restricted to
    erythroid cells only.

    Args:
        bg (BloodGroup): The BloodGroup object to be modified.
        ant_type (PhenoType): The phenotype type (numeric or alphanumeric) determining
            which phenotype mapping to modify.

    Returns:
        BloodGroup: The updated BloodGroup with modified FY phenotypes.
    """
    if bg.type != "FY":
        return bg
    for pair, pheno in bg.phenotypes[ant_type].items():
        erythroid = ["FY*02N.01", "FY*01N.08", "FY*01N.01"]
        if (
            pair.allele1.genotype in erythroid or pair.allele2.genotype in erythroid
        ) and null_or_mod(pair, "N"):
            # only if 2 nulls
            # TODO think about how to handle nulls better
            pheno += "_erythroid_cells_only"
        bg.phenotypes[ant_type][pair] = pheno

    return bg


@apply_to_dict_values
def modify_FY2(bg: BloodGroup, ant_type: PhenoType) -> BloodGroup:
    """
    add FYX to FY*02W.01 and FY*02W.02

    Args:
        bg (BloodGroup): The BloodGroup object to be modified.
        ant_type (PhenoType): The phenotype type (alphanumeric) whose
            corresponding phenotype mapping is to be updated.

    Returns:
        BloodGroup: The updated BloodGroup with modified FY phenotypes.
    """
    if bg.type != "FY":
        return bg
    for pair, pheno in bg.phenotypes[ant_type].items():
        FYX = ["FY*02W.01", "FY*02W.02"]

        if (
            pair.allele1.genotype in FYX or pair.allele2.genotype in FYX
        ) and "b+w)" in pheno:
            pheno += ",Fyx"

        bg.phenotypes[ant_type][pair] = pheno

    return bg


@apply_to_dict_values
def modify_KEL(bg: BloodGroup, ant_type: PhenoType) -> BloodGroup:
    """Modify KEL BloodGroup phenotypes based on allele null/modification patterns.

    For BloodGroup of type "KEL", this function updates the phenotype mapping in
    bg.phenotypes for the provided phenotype type. If an allele pair matches the
    'M.' pattern (as determined by null_or_mod), the phenotype is set to "Kmod".
    Similarly, if an allele pair matches the 'N.' pattern, the phenotype is set to "KO".

    Args:
        bg (BloodGroup): The BloodGroup object to be modified.
        ant_type (PhenoType): The phenotype type (numeric or alphanumeric) whose
            corresponding phenotype mapping is to be updated.

    Returns:
        BloodGroup: The updated BloodGroup with modified KEL phenotypes.
    """
    if bg.type != "KEL":
        return bg
    for pair in bg.phenotypes[ant_type]:
        if null_or_mod(pair, "M"):
            bg.phenotypes[ant_type][pair] = "Kmod"
    for pair in bg.phenotypes[ant_type]:
        if null_or_mod(pair, "N"):
            bg.phenotypes[ant_type][pair] = "KO"

    return bg


@apply_to_dict_values
def modify_CROM(bg: BloodGroup, ant_type: PhenoType) -> BloodGroup:
    """Modify CROM BloodGroup phenotypes based on allele null/modification patterns.

    For BloodGroup of type "CROM", this function updates the phenotype mapping in
    bg.phenotypes for the provided phenotype type.
    if an allele pair matches the 'N.' pattern, the phenotype is changed from
    IFC- to Inab

    Args:
        bg (BloodGroup): The BloodGroup object to be modified.
        ant_type (PhenoType): The phenotype type (numeric or alphanumeric) whose
            corresponding phenotype mapping is to be updated.

    Returns:
        BloodGroup: The updated BloodGroup with modified CROM phenotypes.
    """
    if bg.type != "CROM":
        return bg

    for pair in bg.phenotypes[ant_type]:
        if null_or_mod(pair, "N"):
            bg.phenotypes[ant_type][pair] = "Inab"

    return bg


@apply_to_dict_values
def modify_MNS(bg: BloodGroup, ant_type: PhenoType) -> BloodGroup:
    """
    U+w -> U+var

    Args:
        bg (BloodGroup): The BloodGroup object to be modified.
        ant_type (PhenoType): The phenotype type (alphanumeric) whose
            corresponding phenotype mapping is to be updated.

    Returns:
        BloodGroup: The updated BloodGroup with modified MNS (GYPB) phenotypes.
    """
    if bg.type != "GYPB":
        return bg
    for pair, pheno in bg.phenotypes[ant_type].items():
        weak_mod = "GYPB*04.03"
        var_mod = [
            "GYPB*04N.03",
            "GYPB*04N.01",
            "GYPB*03N.01",
            "GYPB*03N.02",
            "GYPB*03N.03",
            "GYPB*03N.04",
            "GYPB*03N.05",
            "GYPB*03N.06",
            "GYPB*03N.07",
        ]
        if (
            pair.allele1.genotype == weak_mod or pair.allele2.genotype == weak_mod
        ) and (pair.allele1.genotype in var_mod or pair.allele2.genotype in var_mod):
            pheno += "var"
        elif pair.allele1.genotype in var_mod or pair.allele2.genotype in var_mod:
            pheno = pheno.replace("U+w", "U+var")

        bg.phenotypes[ant_type][pair] = pheno

    return bg


@apply_to_dict_values
def modify_RHD(bg: BloodGroup, ant_type: PhenoType) -> BloodGroup:
    """Annotate RHD partial/weak phenotypes with specific variant types.

    Modifies RHD phenotype annotations by adding variant-specific labels to
    partial and weak D antigens. For example, 'D+partial' becomes
    'D+partial(DAU0.01)' based on the underlying genotype.

    Args:
        bg: BloodGroup object containing phenotype and genotype information.
        ant_type: Phenotype representation type (e.g., verbose, short).

    Returns:
        The modified BloodGroup object with annotated RHD phenotypes. Returns
        the original object unchanged if not RHD type or if ant_type is numeric.

    Example:
    'D+partial' -> 'D+partial(DAU0.01)'
    """

    if bg.type != "RHD":
        return bg
    if ant_type == PhenoType.numeric:
        return bg

    for pair, pheno in bg.phenotypes[ant_type].items():
        # Start with all antigens from the original phenotype
        pheno_parts = pheno.split(",")
        new_pheno = []

        for ant in pheno_parts:
            modified = False
            # Check each allele to see if it should annotate this antigen
            for allele in pair:
                if (
                    ant in allele.phenotype_alt.split(",")
                    and allele.genotype in RHD_ANT_MAP
                ):
                    genotype_label = RHD_ANT_MAP[allele.genotype]
                    if "partial" in ant:
                        new_pheno.append(
                            ant.replace("partial", f"partial({genotype_label})")
                        )
                        modified = True
                        break
                    elif "weak" in ant:
                        new_pheno.append(ant.replace("weak", f"weak({genotype_label})"))
                        modified = True
                        break

            # If not modified, keep original
            if not modified:
                new_pheno.append(ant)
        if new_pheno:
            updated_pheno = ",".join(new_pheno)
            # ic(
            #     bg.sample,
            #     pair,
            #     pheno,
            #     updated_pheno,
            #     allele.phenotype_alt,
            #     allele.phenotype_alt.split(","),
            # )
            bg.phenotypes[ant_type][pair] = updated_pheno

    return bg


@apply_to_dict_values
def re_order_KEL(bg: BloodGroup, ant_type: PhenoType) -> BloodGroup:
    """Reorder KEL phenotypes so that K+ and K- appear first.

    This function processes phenotype strings for BloodGroup objects of type "KEL". It
    reorders the antigens in the phenotype string such that antigens starting with "K+"
    or "K-" are placed at the beginning, followed by the remaining antigens. For example,
    a phenotype string "Js(a+b+),K-,k+" will be rearranged to "K-,k+,Js(a+b+)".

    Args:
        bg (BloodGroup): The BloodGroup object containing phenotype mappings.
        ant_type (PhenoType): The phenotype type (e.g., numeric or alphanumeric) to be
            processed.

    Returns:
        BloodGroup: The updated BloodGroup with re-ordered KEL phenotypes.
    """
    if bg.type != "KEL":
        return bg

    big_k = ("K+", "K-")

    def sorter(unsorted_pheno):
        sorted_antigens = []
        for i in range(2):
            for ant in unsorted_pheno.strip().split(","):
                if (
                    i == 0
                    and ant.upper().startswith(big_k)
                    or i == 1
                    and not ant.upper().startswith(big_k)
                ):
                    sorted_antigens.append(ant)

        return ",".join(sorted_antigens)

    for pair, pheno in bg.phenotypes[ant_type].items():
        if not pheno.upper().startswith(big_k):
            pheno = sorter(pheno)
        bg.phenotypes[ant_type][pair] = pheno

    return bg


@apply_to_dict_values
def compare_numeric_ants_to_alphanumeric(
    bg: BloodGroup, mapping: Mapping[str, Mapping[str, str]]
) -> BloodGroup:
    """sanity checker
    Args:
        bg (BloodGroup): The BloodGroup object containing phenotype mappings.
        ant_type (PhenoType): The phenotype type (e.g., numeric or alphanumeric) to be
            processed.

    Returns:
        BloodGroup: The updated BloodGroup with re-ordered KEL phenotypes.
    """
    if (
        bg.phenotypes.get(PhenoType.alphanumeric) == {}
        or bg.phenotypes.get(PhenoType.numeric) == {}
    ):
        return bg
    # if bg.type == 'GBGT':
    #     return bg
    bg_name_map = {
        "GBGT1": "FORS",
        # "GYPA": "MNS",
        # "GYPB": "MNS",
        "ABCC4": "PEL",
        "PIGG": "EMM",
        "GCNT2": "I",
    }
    for pair_alpha, pheno_alpha in bg.phenotypes[PhenoType.alphanumeric].items():
        compare = True
        # if bg.type == "RHCE":
        #     ic(
        #         bg.type,
        #         pair_alpha,
        #         pheno_alpha,
        #         bg.phenotypes[PhenoType.numeric][pair_alpha],
        #     )
        pheno_numeric = bg.phenotypes[PhenoType.numeric][pair_alpha]
        for skip in ["Vel+strong", "erythroid"]:
            if skip in pheno_alpha or skip in pheno_numeric:
                compare = False
        if not compare:
            continue
        if not compare_antigen_profiles(
            pheno_numeric,
            pheno_alpha,
            mapping,
            bg_name_map.get(bg.type, bg.type),
        ):
            logger.warning(
                f"WARNING: Modifiers dont match for sample {bg.sample}:\n{pheno_numeric}\n{pheno_alpha}\n plz report"
            )

    return bg


numeric_calculators = {}

number_of_primary_antitheticals = {
    "RHAG": 1,
    "JMH": 1,
    "OK": 1,
    "KN": 13,
    "CROM": 1,
    "LW": 1,
    "Vel": 1,
    "LU": 0,
}
