from __future__ import annotations

from collections import defaultdict
from functools import partial, reduce
from itertools import zip_longest
from multiprocessing import Pool
from typing import TYPE_CHECKING, Any, Callable

if TYPE_CHECKING:
    from src.core_logic.alleles import Allele, BloodGroup


class BeyondLogicError(Exception):
    """Custom exception for scenarios beyond logical comprehension.

    Attributes:
        message (str): Explanation of the error.
        context (str | None): Additional context or metadata to describe the issue.
    """

    def __init__(self, message: str, context: str | None = None):
        """Initialize BeyondLogicError with a message and optional context.

        Args:
            message (str): The error message describing the beyond logic situation.
            context (str | None): Optional context to provide more information.
        """
        self.message = message
        self.context = context
        super().__init__(self._format_error())

    def _format_error(self) -> str:
        """Format the error message.

        Returns:
            str: Formatted error message including context if available.
        """
        if self.context:
            return f"{self.message} | Context: {self.context}"
        return self.message

    def __str__(self) -> str:
        """Return string representation of the error.

        Returns:
            str: The formatted error message.
        """
        return self._format_error()


class Zygosity:
    HOM = "Homozygous"  # hom alt
    HET = "Heterozygous"
    REF = "Reference"  # hom ref


Preprocessor = Callable[[dict[str, "BloodGroup"]], dict[str, "BloodGroup"]]


def compose(*functions: Preprocessor) -> Preprocessor:
    """Helper function to call all BloodGroup functions sequentially.

    Args:
        *functions (Preprocessor): A variable number of functions that take a
            dictionary of BloodGroup objects and return a modified dictionary of
            BloodGroup objects.

    Returns:
        Preprocessor: A composed function that applies all given functions in
            sequence.
    """
    return reduce(lambda func1, func2: lambda x: func2(func1(x)), functions)


def chunk_list_by_rank(input_list: list[Allele]) -> list[list[Allele]]:
    """Split a list into chunks of Allele objects with the same phenotype weight.

    Args:
        input_list (list[Allele]): A list of Allele objects to be grouped by their
            weight_pheno attribute.

    Returns:
        list[list[Allele]]: A list of lists, where each inner list contains Allele
            objects sharing the same weight_pheno, ordered by increasing weight_pheno.
    """
    result = defaultdict(list)
    for allele in input_list:
        result[allele.weight_pheno].append(allele)

    return [result[i] for i in sorted(result.keys())]


def apply_to_dict_values(func: Callable[..., Any]) -> Callable[..., dict]:
    """
    Decorator to apply a function to all values of a dictionary.

    Args:
        func (Callable[..., Any]): The function to apply to each value.

    Returns:
        Callable[..., Dict]: A new function that takes a dictionary and additional
                             arguments, and returns a dictionary with the function
                             applied to its values.
    """

    def decorator(
        input_dict: dict[str, list[BloodGroup]], *args, **kwargs
    ) -> dict[str, list[BloodGroup]]:
        return {key: func(value, *args, **kwargs) for key, value in input_dict.items()}

    return decorator


def one_HET_variant(allele: Allele, pool: dict[str, str]) -> bool:
    """Check if the allele has one HET variant, all HOM variants, or just one variant.

    Args:
        allele (Allele): An Allele object containing defining variants and
            related attributes.
        pool (dict[str, str]): A dictionary mapping variant identifiers to their
            zygosity. Defaults to "HOM" if a variant is not found.

    Returns:
        bool: True if the allele has exactly one HET variant, if it has only one
            defining variant, or if all defining variants are HOM; otherwise,
            False.
    """
    proper_vars = [
        variant for variant in allele.defining_variants if not variant.endswith(".")
    ]

    return (
        sum([1 for variant in proper_vars if pool.get(variant, "HOM") == Zygosity.HET])
        == 1
        or allele.number_of_defining_variants == 1
        or all(pool.get(variant, "HOM") == Zygosity.HOM for variant in proper_vars)
    )


def one_HET_all_HOM_ref_or_1variant(allele: Allele, pool: dict[str, str]) -> bool:
    """Check if an allele meets one of the following conditions:
    1. It has exactly one HET variant.
    2. It has exactly one defining variant.
    3. It is a reference allele.
    4. All defining variants are HOM.

    Args:
        allele (Allele): An Allele object with its defining variants.
        pool (dict[str, str]): A dictionary mapping variant identifiers to
            their zygosity. Defaults to "HOM" if a variant is not found.

    Returns:
        bool: True if any of the above conditions is met, False otherwise.
    """

    return (
        sum(
            [
                1
                for variant in allele.defining_variants
                if pool.get(variant, "HOM") == Zygosity.HET
            ]
        )
        == 1
        or allele.number_of_defining_variants == 1
        or allele.reference
        or all(
            pool.get(variant, "HOM") == Zygosity.HOM
            for variant in allele.defining_variants
        )
    )


def check_available_variants(
    count: int,
    variant_pool: dict[str, int],
    operator_func: Callable[[int, int], bool],
    allele: Allele,
) -> bool:
    """Check variant counts for a given allele.

    Unaltered counts: 2 indicates HOM (homozygous) and 1 indicates HET 
    (heterozygous). An altered count of 0 means none are left.

    Args:
        count (int): The count value to compare against.
        variant_pool (dict[str, int]): A dictionary mapping variant identifiers
            to their available counts.
        operator_func (Callable[[int, int], bool]): A function that compares the
            available count from the variant pool with the provided count.
        allele (Allele): An Allele object whose defining variants are to be checked.

    Returns:
        bool: A list of boolean results for each defining variant check.
    """
    return [
        operator_func(variant_pool.get(variant, 0), count)
        for variant in allele.defining_variants
    ]


def get_non_refs(options: list[Allele]) -> list[Allele]:
    """Get non-reference alleles.

    Args:
        options (list[Allele]): A list of Allele objects to filter.

    Returns:
        list[Allele]: A list containing only those Allele objects that are not
            reference alleles.
    """

    return [allele for allele in options if not allele.reference]


def chunk_geno_list_by_rank(input_list: list[Allele]) -> list[list[Allele]]:
    """Split a list into chunks of alleles with the same weight.

    Args:
        input_list (List[Allele]): List of Allele objects.

    Returns:
        List[List[Allele]]: List of lists, where each sublist contains alleles of
        the same weight.
    """
    result = defaultdict(lambda: defaultdict(list))
    subs = sorted({a.sub_type for a in input_list})
    weights = sorted({a.weight_geno for a in input_list})

    for allele in input_list:
        result[allele.sub_type][allele.weight_geno].append(allele)

    # Initialize an empty list to hold the final chunks
    final_chunks = []

    # Iterate over each Sub_type and accumulate the weighted alleles
    for sub in subs:
        sub_chunks = [result[sub][weight] for weight in weights if result[sub][weight]]
        if not final_chunks:
            # For the first Sub_type, just assign the chunks directly
            final_chunks = sub_chunks
        else:
            # For subsequent Sub_types, merge the chunks with the existing ones
            final_chunks = [
                x + y
                for x, y in zip_longest(final_chunks, sub_chunks, fillvalue=[])
                if x or y
            ]

    return final_chunks


def sub_alleles_relationships(
    all_alleles: dict[str, list[Allele]], key: str
) -> tuple[dict[str, bool], str]:
    """
    Calculate relationships between alleles under a given key and return the
    relationships along with the key.
    This function iterates over a list of Allele objects, compares each
    pair, and determines if one allele is considered to be 'in' another.
    The results are stored in a dictionary with a custom description as
    the key.

    Args:
        all_alleles (dict[str, list[Allele]]): Dictionary mapping keys to
            lists of Allele objects.
        key (str): The key to access the list of Alleles for comparison.

    Returns:
        tuple[dict[str, bool], str]: A tuple containing a dictionary of
            relationship statuses between alleles and the key.
    """
    relationship = {}
    lst = all_alleles[key]
    for allele1 in lst:
        for allele2 in lst:
            relationship[f"{allele1.genotype}_isin_{allele2.genotype}"] = (
                allele1 in allele2
            )

    return relationship, key


def get_allele_relationships(
    all_alleles: dict[str, list[Allele]], processes: int
) -> dict[str, dict[str, bool]]:
    """
    Compute relationships between alleles across multiple processes. 
    (Was needed when doing & on ABO, left it here incase KN grows)

    Args:
        all_alleles (dict[str, list[Allele]]): Dictionary where keys
            are allele identifiers and values are lists of Allele objects.
        processes (int): Number of processes to use for parallel computation.

    Returns:
        dict[str, dict[str, bool]]: A dictionary where keys are allele
            identifiers and values are dictionaries detailing relationships
            between alleles.

    This function uses multiprocessing to compute relationships
    between alleles in a parallel manner to improve performance.
    Each process handles a subset of the alleles as specified by the
    keys in `all_alleles`.
    """

    relationships = {}
    with Pool(processes=processes) as pool:
        sub_all = partial(sub_alleles_relationships, all_alleles)
        for result, key in pool.imap_unordered(
            sub_all, ["KN"]
        ):  # all_alleles.keys(): (no &)
            relationships[key] = result
    return relationships
