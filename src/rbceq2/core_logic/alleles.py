"""This module contains the Allele class and related functions."""

from __future__ import annotations

import operator
from collections import defaultdict
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any, Callable, Iterator

from loguru import logger

from rbceq2.core_logic.utils import Zygosity, collapse_variant
from rbceq2.core_logic.constants import AlleleState
from frozendict import frozendict

if TYPE_CHECKING:
    from core_logic.constants import PhenoType
    from phenotype.antigens import Antigen


@dataclass(slots=False, frozen=True)
class Allele:
    """A data class representing an Allele.

    Data about genotype, phenotype, and defining_variants of an allele
    ----------

    Attributes
    ----------
    genotype: str
        Genotype of the allele, e.g., "JK*02.01".
    phenotype: str
        Phenotype of the allele.
    genotype_alt: str
        Alternative genotype representation (if any).
    phenotype_alt: str
        Alternative phenotype representation (if any).
    defining_variants: frozenset[str]
        Variants (positions) that define this allele.
    weight_geno: int
        Weight used for ordering or ranking genotype.
    reference: bool
        Whether this allele is the reference allele.
    sub_type: str
        Additional subtype string if needed.
    # phases: tuple[str] | None
    #     Phase IDs associated with this allele (if phased).
    number_of_defining_variants: int
        Automatically set based on the size of defining_variants.
    """

    genotype: str
    phenotype: str
    genotype_alt: str
    phenotype_alt: str
    defining_variants: frozenset[str]
    null: bool
    weight_geno: int = 1
    reference: bool = False
    sub_type: str = ""
    big_variants: frozendict[str, str] = field(default_factory=frozendict)
    number_of_defining_variants: int = field(init=False)

    def __post_init__(self: Allele) -> None:
        """Perform post-initialization tasks.

        This method is automatically called after the object is initialized.
        """
        object.__setattr__(
            self, "number_of_defining_variants", len(self.defining_variants)
        )

    def with_big_variants(self, new: dict[str, str]) -> "Allele":
        """Return a new Allele with updated big_variants."""
        return Allele(
            genotype=self.genotype,
            phenotype=self.phenotype,
            genotype_alt=self.genotype_alt,
            phenotype_alt=self.phenotype_alt,
            defining_variants=self.defining_variants,
            null=self.null,
            weight_geno=self.weight_geno,
            reference=self.reference,
            sub_type=self.sub_type,
            big_variants=frozendict(new),
        )

    def __contains__(self, other: Allele) -> bool:
        """Check if other.defining_variants are subset of self.defining_variants

        Args:
            other Allele: Another Allele object.

        Returns:
            bool: True if all variants in other defining_variants, False otherwise.
        """
        if self.__eq__(other):
            return False  # if compared to self
        if other.reference and self.sub_type == other.sub_type:
            return True  # ref is always 'in' in any child allele
            # especially necesary for KN and any BGs with no SNV for ref
        if self.number_of_defining_variants > other.number_of_defining_variants:
            return other.defining_variants.issubset(self.defining_variants)
        return False

    def __eq__(self, other: Allele) -> bool:
        """Check if another pair object has exactly the same alleles as this one.

        Args:
            other (Allele): Another Allele object to compare against.

        Returns:
            bool: True if both Alleles contain the same geno, pheno and vars,
            False otherwise.
        """
        return (
            self.genotype == other.genotype
            and self.genotype_alt == other.genotype_alt
            and self.phenotype == other.phenotype
            and self.phenotype_alt == other.phenotype_alt
            and self.defining_variants == other.defining_variants
        )

    def __gt__(self, other: Allele) -> bool:
        """Greater than comparison, inverted.

        Args:
            other (Allele) : Another Allele object.

        Returns:
            bool: True if self is greater than other, False otherwise.
        """
        return self._rank(other, operator.gt)

    def __lt__(self, other: Allele) -> bool:
        """Less than comparison, inverted.

        Args:
            other (Allele) : Another Allele object.

        Returns:
            bool: True if self is less than other, False otherwise.
        """
        return self._rank(other, operator.lt)

    def _format_allele(self) -> str:
        """Generate a string representation of the allele."""
        sep_var = "\n\t\t"

        return (
            f"Allele \n "
            f"genotype: {self.genotype} \n "
            f"defining_variants: {sep_var}{sep_var.join([collapse_variant(variant) for variant in self.defining_variants])} \n "
            f"weight_geno: {self.weight_geno} \n "
            f"phenotype: {self.phenotype} or {self.phenotype_alt} \n "
            f"reference: {self.reference} \n"
        )

    def __str__(self) -> str:
        """Return a user-friendly string representation of this Allele."""
        return self._format_allele()

    def __repr__(self) -> str:
        """Return a detailed string representation of this Allele."""
        return self._format_allele()

    def _rank(self, other: Allele, operator_func: Callable[[int, int], bool]) -> bool:
        """Rank one allele vs. another using an operator function.

        Args:
            other: The other allele to compare.
            operator_func: The comparison function, e.g., operator.gt or operator.lt.

        Returns:
            True if the comparison passes, else False.
        """
        if self.sub_type != other.sub_type:
            logger.warning(f"Sub type miss match! {self.sub_type, other.sub_type}")
        if self.weight_geno == other.weight_geno:
            return operator_func(
                self.number_of_defining_variants, other.number_of_defining_variants
            )
        return not operator_func(
            self.weight_geno, other.weight_geno
        )  # inverted because rank 1 is higher weight than 2 etc

    @property
    def blood_group(self) -> str:
        """Get blood group type from genotype.

        Returns:
            str: The blood group of the allele.

        Example:
            KN*01 -> KN
        """
        bg = (
            self.genotype.split("*")[0]
            if self.genotype != "."
            else self.genotype_alt.split("*")[0]
        )
        if "KLF" in bg.upper():
            bg = "KLF"  # most are KLF1, but need to all be the same
        return bg


@dataclass(slots=True, frozen=False)
class Line:
    """A data class used to build allele objects.

    Attributes:
        allele_defining_variants: str
            New variants defining the allele.
        geno: str
            Genotype information.
        pheno: str
            Phenotype information.
        geno_alt: str
            Alternative genotype information.
        pheno_alt: str
            Alternative phenotype information.
        chrom_str: str
            Chromosome information.
        weight_geno: int
            Weight of the genotype.
        weight_pheno: int
            Weight of the phenotype.
        ref: str
            Reference genotype information.
        Sub_type: str
            Subtype of the allele.
        chrom: str
            Chromosome number (only).
    """

    allele_defining_variants: str
    geno: str
    pheno: str
    geno_alt: str
    pheno_alt: str
    chrom_str: str
    weight_geno: int
    ref: str
    sub_type: str
    chrom: str = field(init=False)

    def __post_init__(self):
        """Post-initialization processing to refine chromosome field."""
        self.chrom = self.chrom_str.replace("chr", "")


@dataclass(slots=True, frozen=False)
class BloodGroup:
    type: str
    alleles: dict[
        AlleleState, list[Allele | Pair]
    ]  # fix all ANy as per Arjan's comment
    sample: str
    variant_pool: dict[str, Zygosity] = field(default_factory=dict)
    variant_pool_phase: dict[str, str] = field(default_factory=dict)
    variant_pool_phase_set: dict[str, str] = field(default_factory=dict)
    genotypes: list[str] = field(default_factory=list)
    phenotypes: dict[PhenoType, dict[Pair, list[Antigen]]] = field(
        default_factory=lambda: defaultdict(dict)
    )

    filtered_out: dict[str, list[Allele | Pair]] = field(
        default_factory=lambda: defaultdict(list)
    )
    len_dict: dict[str, int] = field(
        default_factory=lambda: {
            Zygosity.HOM: 2,
            Zygosity.HET: 1,
            Zygosity.REF: 2,
            Zygosity.HEM: 1,
        }
    )
    misc: dict[Any, Any] = None  # TODO - pheno separately??

    """A data class representing a Blood Group
    ----------

    Data about blood group type, genotype, phenotype, and alleles
    to be built up and modified as it runs through the pipeline
    ----------

   Attributes:
        type (str): 
            Blood group type.
        alleles (Dict[str, List[Allele]]): 
            Dictionary mapping allele types to lists of Allele objects.
        sample (str): 
            Sample identifier.
        variant_pool (Dict[str, str]): 
            Mapping of variants to zygosity states.
        variant_pool_phase (Dict[str, str]): 
            Mapping of variants to phase states, ie 1|0.
        variant_pool_phase_set (Dict[str, str]): 
            Mapping of variants to phase sets ie, 126354.
        genotypes (List[str]): 
            List of genotypes associated with the blood group.
        phenotypes (List[str]): 
            List of phenotypes associated with the blood group.
        filtered_out (Dict[str, List[Allele | Pair]]): 
            Alleles filtered out during processing, categorized by reason.
        len_dict (Dict[str, int]): 
            Dictionary mapping zygosity states to their associated numerical values.
        misc (Dict[Any, Any]): 
            Dictionary for miscellaneous stuff. A little promiscuous, probably too much so
    """

    @property
    def variant_pool_numeric(self) -> dict[str, int]:
        """Convert variant pool zygosity states to their numerical values.

        Returns:
            Dict[str, int]: Dictionary of variants mapped to their numerical zygosity values.
        """
        return {k: self.len_dict[v] for k, v in self.variant_pool.items()}

    @property
    def number_of_putative_alleles(self) -> int:
        """Count the number of putative alleles classified as 'raw'.

        Returns:
            int: The number of raw alleles.
        """
        return len(self.alleles[AlleleState.RAW])

    def remove_pairs(
        self, to_remove: list[Pair], filter_name: str, allele_type: str = "pairs"
    ) -> None:
        """Remove pairs of alleles based on specific criteria.

        Args:
            to_remove (list[Pair]]): List of allele pairs to be removed.
            filter_name (str): Category name for the filtering reason.
            allele_type (str): Type of allele group from which pairs are removed.
        """
        already_removed = set()
        if to_remove:
            for pair in to_remove:
                pair_id = ".".join(pair.genotypes)
                if pair_id in already_removed:
                    continue
                self.alleles[allele_type].remove(pair)
                self.filtered_out[filter_name].append(pair)
                already_removed.add(pair_id)
        if not self.alleles[allele_type]:
            logger.warning(
                f"all pairs removed (reverting to reference allele, if possible): {self.sample} {self.type} {filter_name}"
            )

    def remove_alleles(
        self, to_remove: list[str], filter_name: str, allele_type: str = "raw"
    ) -> None:
        """Remove specific alleles from the raw set.

        Args:
            to_remove (List[Allele]): List of alleles to be removed.
            filter_name (str): Category name for the filtering reason.
        """
        for allele in to_remove:
            self.alleles[allele_type].remove(allele)
            self.filtered_out[filter_name].append(allele)
        if not self.alleles[allele_type]:
            logger.warning(
                f"all alleles removed (will revert to reference): {self.sample} {self.type} {filter_name}"
            )


@dataclass(slots=True, frozen=True)
class Pair:
    """Class for managing sets of alleles as immutable pairs.

    Attributes:
        alleles (FrozenSet[Allele]): A set of `Allele` objects that defines this pair.
    """

    allele1: Allele
    allele2: Allele
    alleles: frozenset[Allele] = field(init=False)

    def __post_init__(self) -> None:
        """Perform post-initialization tasks. #TODO sort on ISBT tables once available

        This method is automatically called after the object is initialized.
        """

        object.__setattr__(self, "alleles", frozenset([self.allele1, self.allele2]))

    def __iter__(self) -> Iterator[Allele]:
        """Return an iterator over the alleles in the pair.

        Returns:
            Iterator[Allele]: An iterator over the alleles in the pair.
        """
        return iter([self.allele1, self.allele2])

    def __eq__(self, other: Pair) -> bool:
        """Check if another pair object has exactly the same alleles as this one.

        Args:
            other (Pair): Another pair object to compare against.

        Returns:
            bool: True if both pairs contain the same alleles, False otherwise.
        """
        return self.alleles == other.alleles

    def __contains__(self, other: Allele) -> bool:
        """Check if an allele is in the pair.

        Args:
            other (Allele): Another allele object to check.

        Returns:
            bool: True if the allele is in the pair, False otherwise.
        """
        return other in self.alleles

    def _format_pair(self) -> str:
        """Generate a string representation of the pair.

        Args:
            include_alphanumeric (bool): Whether to include the alphanumeric phenotype.

        Returns:
            str: A formatted string representation of the pair.
        """
        geno = "/".join(self.genotypes)
        pheno = "/".join(self.phenotypes)
        pheno_alpha = "/".join(self.phenotypes_alphanumeric)
        return f"Pair(Genotype: {geno} Phenotype: numeric {pheno} alphanumeric: {pheno_alpha}) #"

    def __repr__(self) -> str:
        """Return a detailed string representation of the pair."""
        return self._format_pair()

    def __str__(self) -> str:
        """Return a user-friendly string representation of the pair."""
        return self._format_pair()

    def _ordered(self):
        """Sort alleles in the pair"""

        return sorted([self.allele1, self.allele2], key=lambda x: x.genotype)

    @property
    def genotypes(self) -> list[str]:
        """Get a list of the genotypes of the alleles in the pair.

        Returns:
            List[str]: A list of the genotypes of the alleles in the pair.
        """
        return [allele.genotype for allele in self._ordered()]

    @property
    def phenotypes(self) -> list[str]:
        """Get a list of the phenotypes of the alleles in the pair.

        Returns:
            List[str]: A list of the genotypes of the alleles in the pair.
        """
        return [allele.phenotype for allele in self._ordered()]

    @property
    def phenotypes_alphanumeric(self) -> list[str]:
        """Get a list of the phenotypes of the alleles in the pair.

        Returns:
            List[str]: A list of the genotypes of the alleles in the pair.
        """
        return [allele.phenotype_alt for allele in self._ordered()]

    @property
    def contains_reference(self) -> bool:
        """Check if the pair contains a reference allele.

        Returns:
            bool: True if the pair contains a reference allele, False otherwise.
        """
        return any(allele.reference for allele in self.alleles)

    @property
    def all_reference(self) -> bool:
        """Check if both alleles in the pair are the reference allele.

        Returns:
            bool: True if the pair contains 2 reference alleles, False otherwise.
        """
        return all(allele.reference for allele in self.alleles)

    @property
    def comparable(self) -> list[frozenset[str]]:  # was parse_bio_info1
        """Parse biological information string into a set of frozensets representing
        each side of '/'.

        Args:
            info (str): Biological information string.

        Returns:
            A set of frozensets, each representing unique substrings
        """
        # Split the info by '/' and then by '+' or '&' to get individual substrings

        return [
            frozenset(allele.genotype.split("+"))
            for allele in [self.allele1, self.allele2]
        ]

    @property
    def same_subtype(self) -> bool:
        """checks if both sub types are equal

        Args:
            info (str):

        Returns:
           bool
        """

        return self.allele1.sub_type == self.allele2.sub_type
