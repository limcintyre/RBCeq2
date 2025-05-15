from __future__ import annotations

from dataclasses import dataclass, field
import importlib.resources
from io import StringIO
from typing import Any, Iterable

import pandas as pd
from rbceq2.core_logic.alleles import Allele, Line
from rbceq2.core_logic.constants import LOW_WEIGHT
from loguru import logger
from icecream import ic
from collections import defaultdict


import re
from abc import abstractmethod
from typing import Mapping, Protocol


class VariantCountMismatchError(ValueError):
    """Exception raised when the number of GRCh37 variants does not match the number of
    GRCh38 variants.

    Attributes
    ----------
    grch37 : str
        The string representation of GRCh37 variants.
    grch38 : str
        The string representation of GRCh38 variants.
    message : str
        Explanation of the error.
    """

    def __init__(self, grch37: str, grch38: str):
        self.grch37 = grch37
        self.grch38 = grch38
        self.message = (
            f"Number of GRCh37 variants must equal the number of GRCh38 variants: "
            f"{grch37} vs {grch38}"
        )
        super().__init__(self.message)


def load_db() -> str:
    """Load the db.tsv file from package resources."""
    # Use importlib.resources.files() which is preferred for Python >= 3.9
    # It needs the package name ('rbceq2') as the anchor.
    try:
        # Correctly reference the resource within the 'rbceq2' package
        resource_path = importlib.resources.files("rbceq2").joinpath(
            "resources", "db.tsv"
        )
        logger.debug(f"Attempting to load db from resource path: {resource_path}")
        return resource_path.read_text(
            encoding="utf-8"
        )  # Always good practice to specify encoding
    except Exception as e:
        logger.error(
            f"Failed to load resource 'resources/db.tsv' from package 'rbceq2': {e}"
        )
        raise  # Re-raise the exception after logging


@dataclass(slots=True, frozen=True)
class Db:
    """A data class representing a genomic database configuration.

    Attributes

        ref (str):
            The reference column name used for querying data within the database.
        df (DataFrame):
            DataFrame loaded from the database file, initialized post-construction.
        lane_variants (dict[str, Any]):
            Dictionary mapping chromosome to its lane variants, initialized
            post-construction.
        antitheticals (dict[str, list[str]]):
            Dictionary mapping blood groups to antithetical alleles, initialized
              post-construction.
        reference_alleles (dict[str, Allele]):
            Dictionary mapping genotype identifiers to reference Allele objects,
            initialized post-construction.
    """

    ref: str
    df: pd.DataFrame #= field(init=False)
    lane_variants: dict[str, Any] = field(init=False)
    antitheticals: dict[str, list[str]] = field(init=False)
    reference_alleles: dict[str, Any] = field(init=False)

    def __post_init__(self):
        #object.__setattr__(self, "df", prepare_db())
        object.__setattr__(self, "antitheticals", self.get_antitheticals())
        object.__setattr__(self, "lane_variants", self.get_lane_variants())
        object.__setattr__(self, "reference_alleles", self.get_reference_allele())
        # self.grch37_38_def_var_count_equal()
        # #self.build_antigen_map()
        # self.check_antigens_are_same()

    # def check_antigens_are_same(self):
    #     """Ensure that the number of antigens in the numeric and alphanumeric descriptions
    #     match and that they haev the same expression and modifiers (weak etc)"""
    #     mapping = self.build_antigen_map()
    #     for num, alpha, sub in zip(
    #         list(self.df.Phenotype_change),
    #         list(self.df.Phenotype_alt_change),
    #         list(self.df.Sub_type),
    #     ):
    #         if num == '.' or alpha == '.':
    #             continue
    #         system = num.strip().split(':')[0]
    #         ic(num, alpha, sub, system)
    #         if system in ['RHD', 'CH', 'RG', 'Ch+Rg+WH+']: #TODO rm C4A
    #             continue
    #         if '?' in num or '?' in alpha:
    #             continue
    #         assert compare_antigen_profiles(
    #             num,
    #             alpha,
    #             mapping,
    #             system,
    #         )


    # def build_antigen_map(self) -> dict[str, dict[str, str]]:
    #     """Build ``{SYSTEM: {numeric_id: canonical_alpha}}`` mapping.

    #     Args:
    #         df: DataFrame with *Phenotype* (numeric) and *Phenotype_alt* (α) columns.

    #     Returns:
    #         Nested mapping suitable for :pyfunc:`compare_antigen_profiles`.

    #     Raises:
    #         ValueError: Any row where token counts differ or parsing fails.
    #     """
    #     mapping: dict[str, dict[str, str]] = defaultdict(dict)
    #     ic(self.df.columns)
    #     for num_raw, α_raw in zip(
    #         self.df.Phenotype, self.df.Phenotype_alt, strict=True
    #     ):
    #         if num_raw == "." or α_raw == ".":
    #             continue  # no cross‑walk for this allele

    #         system, _, num_payload = num_raw.partition(":")
    #         system = system.upper()

    #         num_tokens = [
    #             _NUM_ID_RE.match(tok.strip()).group(1)
    #             for tok in num_payload.split(",")
    #             if tok.strip()
    #         ]
    #         α_tokens = [
    #             _canonical_alpha(tok) for tok in α_raw.split(",") if tok.strip()
    #         ]
    #         ic(num_tokens, α_tokens)
    #         if len(num_tokens) != len(α_tokens):
    #             raise ValueError(
    #                 f"Token mismatch in {system}: "
    #                 f"{len(num_tokens)} numeric vs {len(α_tokens)} alpha"
    #             )

    #         for n, a in zip(num_tokens, α_tokens, strict=True):
    #             mapping[system][n] = a

    #     return mapping

    # def grch37_38_def_var_count_equal(self):
    #     """Ensure that the number of GRCh37 variants == the number of GRCh38 variants
    #     for each allele, and the no of transcript changes"""
    #     for grch37, grch38 in zip(list(self.df.GRCh37), list(self.df.GRCh38)):
    #         if len(grch37.strip().split(",")) != len(grch38.strip().split(",")):
    #             raise VariantCountMismatchError(grch37, grch38)

    # def prepare_db(self) -> pd.DataFrame:
    #     """Read and prepare the database from a TSV file, applying necessary transformations.

    #     Returns:
    #         DataFrame: The prepared DataFrame with necessary data transformations applied.
    #     """
    #     logger.info("Attempting to load database content...")
    #     try:
    #         db_content_str = load_db()
    #         db_content = StringIO(db_content_str)
    #         logger.info("Database content loaded successfully.")
    #     except FileNotFoundError:
    #         logger.error("CRITICAL: db.tsv not found within the package resources!")
    #         # You might want to provide a more informative error or exit here
    #         raise  # Re-raise the specific error
    #     except Exception as e:
    #         logger.error(f"An unexpected error occurred during db loading: {e}")
    #         raise

    #     logger.info("Preparing database from loaded content...")
    #     df: pd.DataFrame = pd.read_csv(db_content, sep="\t")
    #     logger.debug(f"Initial DataFrame shape: {df.shape}")

    #     df["type"] = df.Genotype.apply(lambda x: str(x).split("*")[0])
    #     update_dict = (
    #         df.groupby("Sub_type").agg({"Weight_of_genotype": "max"}).to_dict()
    #     )
    #     mapped_values = df["Sub_type"].map(update_dict)

    #     df["Weight_of_genotype"] = df["Weight_of_genotype"].where(
    #         df["Weight_of_genotype"].notna(), mapped_values
    #     )

    #     pd.set_option("future.no_silent_downcasting", True)
    #     df.Weight_of_genotype = df.Weight_of_genotype.fillna(LOW_WEIGHT)
    #     df.Weight_of_phenotype = df.Weight_of_phenotype.fillna(LOW_WEIGHT)
    #     df = df.fillna(".")
    #     df = df.infer_objects(copy=False)

    #     logger.debug(f"Final DataFrame shape after processing: {df.shape}")
    #     logger.info("Database preparation completed.")
    #     return df

    def get_antitheticals(self) -> dict[str, list[str]]:
        """
        Retrieve antithetical relationships defined in the database.

        Returns:
            dict[str, list[str]]: A dictionary mapping blood groups to lists of
            antithetical alleles.
        """
        antithetical = self.df.query("Antithetical == 'Yes'")
        logger.info(
            f"Antithetical relationships generated: {len(antithetical)} entries."
        )
        return {
            blood_group: list(df[self.ref])
            for blood_group, df in antithetical.groupby("type")
        }

    def get_lane_variants(self) -> dict[str, set[str]]:
        """
        Extract lane variants grouping by chromosome.

        Returns:
            dict[str, set[str]]: A dictionary mapping chromosomes to sets of lane
            variants.
        """

        lane: dict[str, Any] = {}
        for chrom, df in self.df.query("Lane == True").groupby("Chrom"):
            options = {
                sub_variant
                for variant in df[self.ref].unique()
                for sub_variant in variant.split(",")
            }

            lane[chrom] = {
                variant.strip().split("_")[0]
                for variant in options
                if variant.endswith("_ref")
            }
        logger.info(f"Lane positions generated: {len(lane)} entries.")
        return lane

    def line_generator(self, df: pd.DataFrame) -> Iterable[Line]:
        """Yields AlleleData objects from DataFrame columns.

        Args:
            df: DataFrame containing allele data.

        Yields:
            Line objects populated with data from the DataFrame.
        """
        for cols in zip(
            df[self.ref],
            df.Genotype,
            df.Phenotype_change,
            df.Genotype_alt,
            df.Phenotype_alt_change,
            df.Chrom,
            df.Weight_of_genotype,
            #df.Weight_of_phenotype,
            df.Reference_genotype,
            df.Sub_type,
        ):
            yield Line(*cols)

    def get_reference_allele(self) -> dict[str, Allele]:
        """
        Generate reference alleles based on specified criteria.

        Returns:
            Dict[str, Allele]: A dictionary mapping genotype identifiers to reference
            Allele objects.
        """
        refs = self.df.query('Reference_genotype == "Yes"')
        res = {}

        for line in self.line_generator(refs):
            key = line.geno.split("*")[0]
            res[key] = Allele(
                genotype=line.geno,
                phenotype=line.pheno,
                genotype_alt=line.geno_alt,
                phenotype_alt=line.pheno_alt,
                defining_variants=frozenset(
                    [
                        f"{line.chrom}:{a}"
                        for a in line.allele_defining_variants.split(",")
                    ]
                ),
                null=False,
                weight_geno=int(line.weight_geno),
                reference=True,
                sub_type=line.sub_type,
            )
        logger.info(f"Reference alleles generated: {len(res)} entries.")

        return res

    def make_alleles(self) -> Iterable[Allele]:
        """
        Generate Allele objects from the database rows.

        Yields:
            Allele: Allele objects constructed from data rows.
        """
        # df2 = self.df.copy(deep=True)
        # df2["type"] = df2["type"].astype("category")
        # df_ref = df2.loc[(df2["Reference_genotype"] == "Yes") & (df2["type"] == bg.type)]

        # assert df_ref.shape == (1, 20)

        for line in self.line_generator(self.df):
            if line.allele_defining_variants == ".":
                continue
            allele_defining_variants = [
                f"{line.chrom}:{var}"
                for var in map(str.strip, line.allele_defining_variants.split(","))
            ]
            yield Allele(
                line.geno,
                line.pheno,
                line.geno_alt,
                line.pheno_alt,
                frozenset(allele_defining_variants),
                'N.' in line.geno.upper(),
                int(line.weight_geno),
                line.ref == "Yes",
                sub_type=line.sub_type,
            )

    @property
    def unique_variants(self) -> list[str]:
        """
        Compute unique variants from the alleles.

        Returns:
            List[str]: A list of unique variant positions extracted from alleles.
        """
        unique_vars = {
            variant
            for allele in self.make_alleles()
            for variant in allele.defining_variants
        }
        return [f"{pos.split('_')[0]}" for pos in unique_vars]


def prepare_db() -> pd.DataFrame:
    """Read and prepare the database from a TSV file, applying necessary transformations.

    Returns:
        DataFrame: The prepared DataFrame with necessary data transformations applied.
    """
    logger.info("Attempting to load database content...")
    try:
        db_content_str = load_db()
        db_content = StringIO(db_content_str)
        logger.info("Database content loaded successfully.")
    except FileNotFoundError:
        logger.error("CRITICAL: db.tsv not found within the package resources!")
        # You might want to provide a more informative error or exit here
        raise  # Re-raise the specific error
    except Exception as e:
        logger.error(f"An unexpected error occurred during db loading: {e}")
        raise

    logger.info("Preparing database from loaded content...")
    df: pd.DataFrame = pd.read_csv(db_content, sep="\t")
    logger.debug(f"Initial DataFrame shape: {df.shape}")

    df["type"] = df.Genotype.apply(lambda x: str(x).split("*")[0])
    update_dict = (
        df.groupby("Sub_type").agg({"Weight_of_genotype": "max"}).to_dict()
    )
    mapped_values = df["Sub_type"].map(update_dict)

    df["Weight_of_genotype"] = df["Weight_of_genotype"].where(
        df["Weight_of_genotype"].notna(), mapped_values
    )

    pd.set_option("future.no_silent_downcasting", True)
    df.Weight_of_genotype = df.Weight_of_genotype.fillna(LOW_WEIGHT)
    #df.Weight_of_phenotype = df.Weight_of_phenotype.fillna(LOW_WEIGHT)
    df = df.fillna(".")
    df = df.infer_objects(copy=False)

    logger.debug(f"Final DataFrame shape after processing: {df.shape}")
    logger.info("Database preparation completed.")
    return df
######new ####



@dataclass(slots=True, frozen=True)
class DbInternalConsistencyCheck:
    """A data class representing a genomic database configuration.

    Attributes

        ref (str):
            The reference column name used for querying data within the database.
        df (DataFrame):
            DataFrame loaded from the database file, initialized post-construction.
        lane_variants (dict[str, Any]):
            Dictionary mapping chromosome to its lane variants, initialized
            post-construction.
        antitheticals (dict[str, list[str]]):
            Dictionary mapping blood groups to antithetical alleles, initialized
              post-construction.
        reference_alleles (dict[str, Allele]):
            Dictionary mapping genotype identifiers to reference Allele objects,
            initialized post-construction.
    """

    ref: str
    df: pd.DataFrame #= field(init=False)
    lane_variants: dict[str, Any] = field(init=False)
    antitheticals: dict[str, list[str]] = field(init=False)
    reference_alleles: dict[str, Any] = field(init=False)

    def __post_init__(self):
        #object.__setattr__(self, "df", prepare_db())
        # object.__setattr__(self, "antitheticals", self.get_antitheticals())
        # object.__setattr__(self, "lane_variants", self.get_lane_variants())
        # object.__setattr__(self, "reference_alleles", self.get_reference_allele())
        self.grch37_38_def_var_count_equal()
        #self.build_antigen_map()
        self.check_antigens_are_same()

    def check_antigens_are_same(self):
        """Ensure that the number of antigens in the numeric and alphanumeric descriptions
        match and that they haev the same expression and modifiers (weak etc)"""
        mapping = self.build_antigen_map()
        for num, alpha, sub in zip(
            list(self.df.Phenotype_change),
            list(self.df.Phenotype_alt_change),
            list(self.df.Sub_type),
        ):
            if num == '.' or alpha == '.':
                continue
            system = num.strip().split(':')[0]
            ic(num, alpha, sub, system)
            if system in ['RHD', 'CH', 'RG', 'Ch+Rg+WH+']: #TODO rm C4A
                continue
            if '?' in num or '?' in alpha:
                continue
            assert compare_antigen_profiles(
                num,
                alpha,
                mapping,
                system,
            )


    def build_antigen_map(self) -> dict[str, dict[str, str]]:
        """Build ``{SYSTEM: {numeric_id: canonical_alpha}}`` mapping.

        Args:
            df: DataFrame with *Phenotype* (numeric) and *Phenotype_alt* (α) columns.

        Returns:
            Nested mapping suitable for :pyfunc:`compare_antigen_profiles`.

        Raises:
            ValueError: Any row where token counts differ or parsing fails.
        """
        mapping: dict[str, dict[str, str]] = defaultdict(dict)
        for num_raw, α_raw in zip(
            self.df.Phenotype, self.df.Phenotype_alt, strict=True
        ):
            if num_raw == "." or α_raw == ".":
                continue  # no cross‑walk for this allele

            system, _, num_payload = num_raw.partition(":")
            system = system.upper()

            num_tokens = [
                _NUM_ID_RE.match(tok.strip()).group(1)
                for tok in num_payload.split(",")
                if tok.strip()
            ]
            α_tokens = [
                _canonical_alpha(tok) for tok in α_raw.split(",") if tok.strip()
            ]
            if len(num_tokens) != len(α_tokens):
                raise ValueError(
                    f"Token mismatch in {system}: "
                    f"{len(num_tokens)} numeric vs {len(α_tokens)} alpha"
                )

            for n, a in zip(num_tokens, α_tokens, strict=True):
                mapping[system][n] = a

        return mapping

    def grch37_38_def_var_count_equal(self):
        """Ensure that the number of GRCh37 variants == the number of GRCh38 variants
        for each allele, and the no of transcript changes"""
        for grch37, grch38 in zip(list(self.df.GRCh37), list(self.df.GRCh38)):
            if len(grch37.strip().split(",")) != len(grch38.strip().split(",")):
                raise VariantCountMismatchError(grch37, grch38)



# ────────────────────── helper regexes ──────────────────────
_NUM_ID_RE = re.compile(r"-?(\d+)")  # leading '-' allowed
_ALPHA_CANON_RE = re.compile(r"^(.*?)\s|[+-]", re.S)  # up‑to first space/+/‑


# ────────────────────── internal helpers ───────────────────
def _canonical_alpha(token: str) -> str:
    """Return antigen name stripped of sign/modifiers."""
    # stop at first space or sign; then strip trailing sign if still present
    cut = _ALPHA_CANON_RE.split(token.strip(), maxsplit=1)[0]
    return cut.rstrip("+-")


@dataclass(frozen=True, slots=True)
class Antigen:
    """Canonical antigen description.

    Attributes
    ----------
    system : str
        Blood‑group system code (RH, LU …).
    name : str
        Antigen name in its canonical (α) form – e.g. ``"C"`` or ``"Lu4"``.
    expressed : bool
        ``True`` → positive ( ‘+’ ); ``False`` → negative ( ‘‑’ ).
    modifiers : frozenset[str]
        One‑letter modifier codes {'w', 'p', 'n', 'm', 'i', 'r'}.
    """

    system: str
    name: str
    expressed: bool
    modifiers: frozenset[str]


# ──────────────────────────────── parsing ────────────────────────────────


class AntigenParser(Protocol):
    """A parser returns a sequence of canonical :class:`Antigen` objects."""

    @abstractmethod
    def parse(self, text: str) -> list[Antigen]: ...


_NUMERIC_RE = re.compile(r"(?P<sign>-)?(?P<num>\d+)(?P<mods>[a-z]+)?", re.IGNORECASE)
_ALPHA_MOD = {
    "weak": "w",
    "very_weak": "v",
    "partial": "p",
    "neg": "n",
    "negative": "n",
    "monoclonal": "m",
    "inferred": "i",
    "robust": "r",
    "strong": "s",
    "positive_to_neg": "n",
    "weak_to_neg": "n",
    "very_weak_to_neg": "n",
}
# # Assume _ALPHA_MOD and _MOD_LETTERS are defined as in your provided code:
# # _ALPHA_MOD = {
# #     "weak": "w", "very_weak": "v", "partial": "p", "neg": "n", "negative": "n",
# #     "monoclonal": "m", "inferred": "i", "robust": "r", "strong": "s",
# #     "positive_to_neg": "n", "weak_to_neg": "n", "very_weak_to_neg": "n",
# # }
# # _MOD_LETTERS = set(_ALPHA_MOD.values())


# NEW – valid single‑letter codes
_MOD_LETTERS = set(_ALPHA_MOD.values())

class NumericParser:
    """Convert *numeric* antigen strings into :class:`Antigen` objects."""

    def __init__(self, system: str):
        self._system = system.upper()

    def parse(self, text: str) -> list[Antigen]:
        if text.strip() in {"", "."}:
            return []

        _, _, tokens = text.partition(":")  # ignore leading "RH:" / "LU:"
        antigens: list[Antigen] = []

        for raw in tokens.split(","):
            tok = raw.strip()
            if not tok:
                continue

            m = _NUMERIC_RE.fullmatch(tok)
            if not m:
                raise ValueError(f"Bad numeric token: {tok}")

            sign, num, mods = m.group("sign", "num", "mods")
            antigens.append(
                Antigen(
                    system=self._system,
                    name=num,  # mapped to α later
                    expressed=sign != "-",
                    modifiers=frozenset(mods or ""),
                )
            )
        return antigens

# import re # Ensure re is imported

# # Assume _ALPHA_MOD and _MOD_LETTERS are defined as in your provided code:
# # _ALPHA_MOD = {
# #     "weak": "w", "very_weak": "v", "partial": "p", "neg": "n", "negative": "n",
# #     "monoclonal": "m", "inferred": "i", "robust": "r", "strong": "s",
# #     "positive_to_neg": "n", "weak_to_neg": "n", "very_weak_to_neg": "n",
# # }
# # _MOD_LETTERS = set(_ALPHA_MOD.values())


class AlphaParser:
    """Convert *alphanumeric* antigen strings into :class:`Antigen` objects."""

    def __init__(self, system: str):
        self._system = system.upper()

    def parse(self, text: str) -> list[Antigen]:
        if text.strip() in {"", "."}:
            return []

        antigens: list[Antigen] = []

        for raw in text.split(","):
            tok = raw.strip()
            if not tok:
                continue
                
            try:
                idx = next(i for i, ch in enumerate(tok) if ch in "+-")
            except StopIteration:
                raise ValueError(f"Missing +/- in token: {tok}")

            # name_part_before_rstrip is the segment of the token before the sign.
            # E.g., for "Fy(a+w)", this is "Fy(a".
            # E.g., for "K+", this is "K".
            name_part_before_rstrip = tok[:idx]
            
            name = name_part_before_rstrip.rstrip(" (") # Final antigen name, e.g., "Fy(a)" or "K"
            expr = tok[idx] == "+"
            
            # Initial tail, may contain structural parentheses.
            # E.g., for "Fy(a+w)", tail is "w)". For "Fy(b-)", tail is ")". For "K+", tail is "".
            tail = tok[idx + 1 :].lower()
            
            # --- MINIMAL FIX FOR TRAILING PARENTHESES IN MODIFIERS ---
            # If the tail ends with ')' and the part of the token that formed the 'name'
            # had an unmatched opening parenthesis, remove the ')' from the tail.
            # This handles cases like "Fy(a+w)" -> tail "w)" becoming "w",
            # and "Fy(b-)" -> tail ")" becoming "".
            if tail.endswith(')') and \
               name_part_before_rstrip.count('(') > name_part_before_rstrip.count(')'):
                tail = tail[:-1]
            # --- END OF MINIMAL FIX ---
                
            mods: set[str] = set()
            
            # Existing modifier parsing logic (should now work with cleaned tail)
            components = [comp for comp in re.split(r"\s+", tail.strip()) if comp]

            for comp in components:
                comp_code = _ALPHA_MOD.get(comp)
                if comp_code:
                    mods.add(comp_code)
                    continue

                if "_" in comp:
                    found_mods_in_underscore_split = False
                    for part in comp.split("_"):
                        part_code = _ALPHA_MOD.get(part)
                        if part_code:
                            mods.add(part_code)
                            found_mods_in_underscore_split = True
                    if found_mods_in_underscore_split:
                        continue
                
                # Check _MOD_LETTERS is not empty before using it in set operations
                # to prevent error if it was somehow misconfigured.
                if _MOD_LETTERS and set(comp) <= _MOD_LETTERS:
                    mods.update(list(comp))
                    continue
            
            antigens.append(
                Antigen(
                    system=self._system,
                    name=name,
                    expressed=expr,
                    modifiers=frozenset(mods),
                )
            )
        return antigens


    # def parse(self, text: str) -> list[Antigen]:
    #     if text.strip() in {"", "."}:
    #         return []

    #     antigens: list[Antigen] = []

    #     for raw in text.split(","):
    #         tok = raw.strip()
    #         # locate first + or -
    #         try:
    #             idx = next(i for i, ch in enumerate(tok) if ch in "+-")
    #         except StopIteration:
    #             raise ValueError(f"Missing +/- in token: {tok}")

    #         name = tok[:idx].rstrip(" (")
    #         expr = tok[idx] == "+"
    #         tail = tok[idx + 1 :].lower()

    #         mods: set[str] = set()
    #         for word in re.split(r"\W+", tail):
    #             #for word in re.split(r"[_\W]+", tail):
    #             code = _ALPHA_MOD.get(word)
    #             if code:
    #                 mods.add(code)

       
    #         first_seq = re.match(r"[a-z]+", tail.lstrip())
    #         if first_seq and set(first_seq.group(0)) <= _MOD_LETTERS:
    #             mods.update(first_seq.group(0))          # e.g. "+w", "+pw", "+pwn"

    #         antigens.append(
    #             Antigen(
    #                 system=self._system,
    #                 name=name,
    #                 expressed=expr,
    #                 modifiers=frozenset(mods),
    #             )
    #         )
    #     return antigens


# ──────────────────────────────── comparison ────────────────────────────────


def compare_antigen_profiles(
    numeric: str,
    alpha: str,
    mapping: Mapping[str, Mapping[str, str]],
    system: str,
    *,
    strict: bool = True,
) -> bool:
    """Return *True* when the profiles are equivalent.

    Args:
        numeric: Numeric string, e.g. ``"RH:2w,-3"``.
        alpha:   Alphanumeric string, e.g. ``"C+ weak,E-"``.
        mapping: ``{"RH": {"2": "C", ...}, "LU": ...}``.
        system:  Blood‑group code (RH, LU …).
        strict:  Require one‑to‑one match; ``False`` allows extras.

    Raises:
        ValueError: Unknown antigen or malformed token.
    """
    num_ants = NumericParser(system).parse(numeric)
    α_ants = AlphaParser(system).parse(alpha)
    # translate numeric → canonical α‑name
    num_by_name: dict[str, Antigen] = {}
    sys_map = mapping.get(system.upper(), {})
    for n in num_ants:
        try:
            α = sys_map[n.name]
        except KeyError as exc:
            raise ValueError(f"Missing map for {system}:{n.name}") from exc
        num_by_name[α] = Antigen(
            system=n.system,
            name=α,
            expressed=n.expressed,
            modifiers=n.modifiers,
        )

    seen = set()
    for a in α_ants:
        counterpart = num_by_name.get(a.name)
        if not counterpart:
            if strict:
                return False
            continue
        seen.add(a.name)
        if (a.expressed != counterpart.expressed) or (
            a.modifiers != counterpart.modifiers
        ):
            return False

    return not (strict and (seen != set(num_by_name)))

