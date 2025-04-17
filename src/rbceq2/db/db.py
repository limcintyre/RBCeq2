from __future__ import annotations

from dataclasses import dataclass, field
import importlib.resources
from io import StringIO
from typing import Any, Iterable

import pandas as pd
from rbceq2.core_logic.alleles import Allele, Line
from rbceq2.core_logic.constants import LOW_WEIGHT
from loguru import logger


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
        resource_path = importlib.resources.files("rbceq2").joinpath("resources", "db.tsv")
        logger.debug(f"Attempting to load db from resource path: {resource_path}")
        return resource_path.read_text(encoding="utf-8") # Always good practice to specify encoding
    except Exception as e:
        logger.error(f"Failed to load resource 'resources/db.tsv' from package 'rbceq2': {e}")
        raise # Re-raise the exception after logging


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
    df: pd.DataFrame = field(init=False)
    lane_variants: dict[str, Any] = field(init=False)
    antitheticals: dict[str, list[str]] = field(init=False)
    reference_alleles: dict[str, Any] = field(init=False)

    def __post_init__(self):
        object.__setattr__(self, "df", self.prepare_db())
        object.__setattr__(self, "antitheticals", self.get_antitheticals())
        object.__setattr__(self, "lane_variants", self.get_lane_variants())
        object.__setattr__(self, "reference_alleles", self.get_reference_allele())
        self.grch37_38_def_var_count_equal()
        # TODO same for antigen counts

    def grch37_38_def_var_count_equal(self):
        """Ensure that the number of GRCh37 variants == the number of GRCh38 variants
        for each allele, and the no of transcript changes"""
        for grch37, grch38 in zip(list(self.df.GRCh37), list(self.df.GRCh38)):
            if len(grch37.strip().split(",")) != len(grch38.strip().split(",")):
                raise VariantCountMismatchError(grch37, grch38)

    def prepare_db(self) -> pd.DataFrame:
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
            raise # Re-raise the specific error
        except Exception as e:
            logger.error(f"An unexpected error occurred during db loading: {e}")
            raise

        logger.info("Preparing database from loaded content...")
        df: pd.DataFrame = pd.read_csv(db_content, sep="\t")
        logger.debug(f"Initial DataFrame shape: {df.shape}")

        # db_content = StringIO(load_db())
        # # current_dir = Path(__file__).resolve().parent.parent.parent
        # # db_path = current_dir / "resources" / "db.tsv"
        # # logger.info(f"Preparing database from path: {db_path}")
        # df: pd.DataFrame = pd.read_csv(db_content, sep="\t")
        # logger.debug(f"Initial DataFrame shape: {df.shape}")

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
        df.Weight_of_phenotype = df.Weight_of_phenotype.fillna(LOW_WEIGHT)
        df = df.fillna(".")
        df = df.infer_objects(copy=False)

        logger.debug(f"Final DataFrame shape after processing: {df.shape}")
        logger.info("Database preparation completed.")
        return df

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
            df.Weight_of_phenotype,
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
                weight_geno=int(line.weight_geno),
                weight_pheno=int(line.weight_pheno),
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
                int(line.weight_geno),
                int(line.weight_pheno),
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
