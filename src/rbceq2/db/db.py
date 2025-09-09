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
        resource_path = importlib.resources.files("rbceq2").joinpath(
            "resources", "db.tsv"
        )
        logger.debug(f"Attempting to load db from resource path: {resource_path}")
        return resource_path.read_text(
            encoding="utf-8"
        )
    except Exception as e:
        logger.error(
            f"Failed to load resource 'resources/db.tsv' from package 'rbceq2': {e}"
        )
        raise


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
    df: pd.DataFrame 
    lane_variants: dict[str, Any] = field(init=False)
    antitheticals: dict[str, list[str]] = field(init=False)
    reference_alleles: dict[str, Any] = field(init=False)

    def __post_init__(self):
        object.__setattr__(self, "antitheticals", self.get_antitheticals())
        object.__setattr__(self, "lane_variants", self.get_lane_variants())
        object.__setattr__(self, "reference_alleles", self.get_reference_allele())
        

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
        antithetical = self.df.query("Antithetical == 'Yes'")
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
        
        ic(1111,lane, self.antitheticals, antithetical)
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
            bg_key = line.geno.split("*")[0]
            if 'KLF' in bg_key:
                bg_key = 'KLF'
            res[bg_key] = Allele(
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
                _is_null_genotype(line.geno),
                int(line.weight_geno),
                line.ref == "Yes",
                sub_type=line.sub_type,
            )

    @property
    def unique_variants(self) -> set[str]:
        """
        Compute unique variants from the alleles.

        Returns:
            set[str]: A set of unique variant positions extracted from alleles.
        """
        unique_vars = {
            variant
            for allele in self.make_alleles()
            for variant in allele.defining_variants
        }
        lanes = []
        for chrom, poses in self.lane_variants.items():
            for pos in poses:
                lanes.append(f"{chrom.replace('chr','')}:{pos}")
        return set([f"{pos.split('_')[0]}" for pos in unique_vars] + lanes)

def _is_null_genotype(genotype: str) -> bool:
    """
    Determines if a blood group genotype is a null allele.

    According to blood group nomenclature, a null allele is typically indicated
    by the presence of "N." in the allele designation or by ending in 'N'.
    This function checks for these two patterns.

    Args:
        genotype: A string representing the blood group genotype
                (e.g., 'JK*02N.22', 'FY*01', 'AUG*01N').

    Returns:
        True if the genotype is identified as null, False otherwise.
    """
    geno_upper = genotype.upper()
    return "N." in geno_upper or geno_upper.endswith('N')

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
    
    #defaults weights; null = LOW_WEIGHT/2 and normal = LOW_WEIGHT
    is_null_mask = df['Genotype'].apply(_is_null_genotype)
    df.loc[is_null_mask & df['Weight_of_genotype'].isnull(), 'Weight_of_genotype'] = LOW_WEIGHT / 2
    df['Weight_of_genotype'] = df['Weight_of_genotype'].fillna(LOW_WEIGHT)

    df = df.fillna(".")
    df = df.infer_objects(copy=False)

    logger.debug(f"Final DataFrame shape after processing: {df.shape}")
    logger.info("Database preparation completed.")
    df.loc[df['type'] == 'KLF1', 'type'] = 'KLF'
    return df


class DbDataConsistencyChecker:
    """
    Provides static methods to check the internal consistency of the database DataFrame
    before it's used to create a Db object.
    """

    @staticmethod
    def _build_antigen_map_for_checks(df: pd.DataFrame) -> dict[str, dict[str, str]]:
        """
        Build ``{SYSTEM: {numeric_id: canonical_alpha}}`` mapping from DataFrame.
        Helper for antigen consistency checks.
        (This is the logic from your original build_antigen_map method,
         now static and taking df as a parameter.)
        """
        mapping: dict[str, dict[str, str]] = defaultdict(dict)
        for num_raw, α_raw in zip(
            df.Phenotype, df.Phenotype_alt, strict=True
        ):
            if num_raw == "." or α_raw == ".":
                continue

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
                    f"Token mismatch in {system} for Phenotype/Phenotype_alt mapping: "
                    f"{len(num_tokens)} numeric vs {len(α_tokens)} alpha ({num_raw} / {α_raw})"
                )

            for n, a in zip(num_tokens, α_tokens, strict=True):
                mapping[system][n] = a
        return mapping

    @staticmethod
    def _perform_antigen_consistency_check(
        df: pd.DataFrame,
        phenotype_col_name: str,
        phenotype_alt_col_name: str
    ):
        """
        Helper to check antigen consistency for a given pair of phenotype columns.
        Uses the mapping built from primary Phenotype/Phenotype_alt columns.
        """
        # The mapping should ideally be built once from the canonical Phenotype columns
        # and then used for all checks.
        # Assuming _build_antigen_map_for_checks uses df.Phenotype and df.Phenotype_alt
        # as the source of truth for the mapping.
        mapping = DbDataConsistencyChecker._build_antigen_map_for_checks(df)
        
        phenotype_series = df[phenotype_col_name]
        phenotype_alt_series = df[phenotype_alt_col_name]

        for i, (num, alpha) in enumerate(zip(list(phenotype_series), list(phenotype_alt_series))):
            if num == '.' or alpha == '.':
                continue
            
            current_system_from_num = num.strip().split(':')[0]
            
            if current_system_from_num in ['RHD', 'CH', 'RG', 'Ch+Rg+WH+']: #TODO rm C4A
                continue
            if '?' in num or '?' in alpha:
                continue

            # Ensure consistent system context for compare_antigen_profiles
            # Usually, the system is derived from the numeric part.
            system_context = current_system_from_num

            if not compare_antigen_profiles(
                numeric=num, # Use the full string e.g. "RH:1+"
                alpha=alpha,   # Use the full string e.g. "D+"
                mapping=mapping,       # The global mapping
                system=system_context, # System derived from the numeric string
            ):
                allele_info = df.loc[i, "Genotype"] if "Genotype" in df.columns else f"Row index {i}"
                error_msg = (
                    f"Antigen profile mismatch for allele '{allele_info}' (System: {system_context}) "
                    f"between '{phenotype_col_name}' ('{num}') and "
                    f"'{phenotype_alt_col_name}' ('{alpha}')"
                )
                # Consider raising a more specific error if desired
                raise AssertionError(error_msg)

    @staticmethod
    def check_phenotype_change_antigens(df: pd.DataFrame):
        """Ensure consistency between Phenotype_change and Phenotype_alt_change."""
        logger.debug("Checking antigen consistency for Phenotype_change / Phenotype_alt_change...")
        DbDataConsistencyChecker._perform_antigen_consistency_check(
            df, "Phenotype_change", "Phenotype_alt_change"
        )
        logger.debug("Phenotype_change / Phenotype_alt_change antigen consistency check passed.")

    @staticmethod
    def check_phenotype_antigens(df: pd.DataFrame):
        """Ensure consistency between Phenotype and Phenotype_alt."""
        logger.debug("Checking antigen consistency for Phenotype / Phenotype_alt...")
        DbDataConsistencyChecker._perform_antigen_consistency_check(
            df, "Phenotype", "Phenotype_alt"
        )
        logger.debug("Phenotype / Phenotype_alt antigen consistency check passed.")

    @staticmethod
    def check_grch37_38_variant_counts(df: pd.DataFrame):
        """Ensure GRCh37 and GRCh38 variant counts match for each allele."""
        logger.debug("Checking GRCh37/38 defining variant counts...")
        for index, row in df.iterrows(): # Iterate for potentially better error context
            grch37_vars_str = str(row.GRCh37)
            grch38_vars_str = str(row.GRCh38)
            
            # Handle potential empty strings or "." consistently before splitting
            grch37_list = [v for v in grch37_vars_str.strip().split(",") if v and v != "."]
            grch38_list = [v for v in grch38_vars_str.strip().split(",") if v and v != "."]

            if len(grch37_list) != len(grch38_list):
                # The VariantCountMismatchError takes the raw strings
                raise VariantCountMismatchError(grch37_vars_str, grch38_vars_str)
        logger.debug("GRCh37/38 defining variant count check passed.")

    @staticmethod
    def run_all_checks(df: pd.DataFrame, ref_genome_name: str | None = None):
        """
        Runs all internal consistency checks on the DataFrame.
        
        Args:
            df (pd.DataFrame): The DataFrame to check.
            ref_genome_name (str | None, optional): The reference genome name (e.g., "GRCh37").
                Currently unused by these specific checks but included for future extensibility
                if some checks become dependent on the reference context.
        """
        logger.info("Running database data consistency checks...")
        DbDataConsistencyChecker.check_grch37_38_variant_counts(df)
        DbDataConsistencyChecker.check_phenotype_change_antigens(df)
        DbDataConsistencyChecker.check_phenotype_antigens(df)
        
        # Example of how you might use ref_genome_name if a check needed it:
        # if ref_genome_name:
        #     DbDataConsistencyChecker.some_check_dependent_on_ref(df, ref_genome_name)
            
        logger.info("All database data consistency checks passed successfully.")

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

OVERRIDING_INTENSITY_PHRASES = {
    "very_weak": "v",}

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

            name_part_before_rstrip = tok[:idx]
            name = name_part_before_rstrip.rstrip(" (")
            expr = tok[idx] == "+"
            tail = tok[idx + 1 :].lower()

            if tail.endswith(')') and \
               name_part_before_rstrip.count('(') > name_part_before_rstrip.count(')'):
                tail = tail[:-1]
                
            current_antigen_mods: set[str] = set() # Mods for the current antigen token (e.g. e+partial_weak_to_neg)
            
            # Split tail into space-separated components (e.g., "very_weak", "partial", "wp")
            components = [comp for comp in re.split(r"\s+", tail.strip()) if comp]

            for comp in components: # Process each component (e.g., "partial_weak_to_neg")
                # 1. Check for overriding intensity phrases first
                overriding_code = OVERRIDING_INTENSITY_PHRASES.get(comp)
                if overriding_code:
                    current_antigen_mods.add(overriding_code)
                    continue # This component is fully handled by the overriding phrase

                # 2. If not an overriding phrase, accumulate modifiers from:
                #    a) The direct match of the component in _ALPHA_MOD
                #    b) Underscore-separated parts of the component found in _ALPHA_MOD
                #    (This allows "weak_to_neg" (if mapped to "n") to also pick up "w" from its "weak" part)
                
                component_processed_by_phrase_or_parts = False

                # 2a. Direct match of the whole component
                direct_comp_code = _ALPHA_MOD.get(comp)
                if direct_comp_code:
                    current_antigen_mods.add(direct_comp_code)
                    component_processed_by_phrase_or_parts = True
                
                # 2b. Underscore-separated parts
                if "_" in comp:
                    for part in comp.split("_"):
                        part_code = _ALPHA_MOD.get(part)
                        if part_code:
                            current_antigen_mods.add(part_code)
                            component_processed_by_phrase_or_parts = True # Mark as processed if any part matches
                
                # 3. If the component was NOT processed by direct phrase match (2a)
                #    NOR by underscore parts (2b effectively, because if it had parts, flag would be true),
                #    AND it does not contain underscores itself (ensuring it's a candidate for single letters),
                #    THEN try to parse as concatenated single modifier letters.
                if not component_processed_by_phrase_or_parts and not ("_" in comp):
                    if _MOD_LETTERS and set(comp) <= _MOD_LETTERS:
                        current_antigen_mods.update(list(comp))
            
            antigens.append(
                Antigen(
                    system=self._system,
                    name=name,
                    expressed=expr,
                    modifiers=frozenset(current_antigen_mods),
                )
            )
        return antigens


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

