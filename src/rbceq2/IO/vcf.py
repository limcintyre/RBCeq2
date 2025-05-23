import gzip
import io
import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import pandas as pd

os.environ["POLARS_MAX_THREADS"] = "7"  # Must be set before polars import
import polars as pl
from loguru import logger

from rbceq2.core_logic.constants import COMMON_COLS, HOM_REF_DUMMY_QUAL


@dataclass(slots=True, frozen=False)
class VCF:
    """A data class to process VCF files for variant calling and analysis.

    Attributes
        input_vcf (Path | pd.DataFrame):
            The input VCF file path or DataFrame.
        lane_variants (dict[str, Any]):
            Mapping of chromosome to variants specific to lanes.
        unique_variants (set[str]):
            A set of unique variant identifiers.
    """

    input_vcf: Path | pd.DataFrame | None
    lane_variants: dict[str, Any]
    unique_variants: set[str]
    sample: str = field(init=False)
    df: pd.DataFrame = field(init=False)
    loci: set[str] = field(init=False)
    variants: dict[str, str] = field(init=False)

    def __post_init__(self):
        """Handle initialization after data class creation."""
        object.__setattr__(self, "df", self.handle_single_or_multi())
        object.__setattr__(self, "sample", self.get_sample())
        self.rename_chrom()
        self.remove_home_ref()
        self.encode_variants()
        self.add_loci()
        object.__setattr__(self, "loci", self.set_loci())
        self.add_lane_variants()
        object.__setattr__(self, "variants", self.get_variants())

    def handle_single_or_multi(self) -> pd.DataFrame:
        """Handle single or multiple entries in the VCF, returning a DataFrame.

        Returns:
            pd.DataFrame: The DataFrame representation of the VCF data.
        """
        if isinstance(self.input_vcf, Path):
            vcf = read_vcf(self.input_vcf)
            return filter_VCF_to_BG_variants(vcf, self.unique_variants)
        else:
            return self.input_vcf[0]

    def rename_chrom(self) -> None:
        """Rename chromosome identifiers by removing the 'chr' prefix."""
        self.df["CHROM"] = self.df["CHROM"].apply(lambda x: x.replace("chr", ""))

    def remove_home_ref(self) -> None:
        """Remove homozygous reference calls from the DataFrame."""
        self.df = self.df[~self.df["SAMPLE"].str.startswith("0/0")].copy(deep=True)

    def encode_variants(self) -> None:
        """Encode variants into a unified format in the DataFrame."""

        def join_vars(chrom: str, pos: str, ref: str, alts: str) -> bool:
            return ",".join([f"{chrom}:{pos}_{ref}_{alt}" for alt in alts.split(",")])

        self.df["variant"] = self.df.apply(
            lambda x: join_vars(
                x["CHROM"],
                x["POS"],
                x["REF"],
                x["ALT"],
            ),
            axis=1,
        )

    def add_loci(self) -> None:
        """Add loci identifiers to the DataFrame."""
        self.df["loci"] = self.df.CHROM + ":" + self.df.POS

    def get_sample(self) -> str:
        """Determine the sample name from the VCF path or DataFrame.

        Returns:
            str: The sample name.
        """
        if isinstance(self.input_vcf, Path):
            return self.input_vcf.stem
        else:
            return self.input_vcf[-1]

    def set_loci(self) -> set[str]:
        """Create a set of loci identifiers from the DataFrame.

        Returns:
            set[str]: The set of loci identifiers.
        """
        return set(self.df.loci)

    def add_lane_variants(self) -> None:
        """Add lane-specific variants to the DataFrame based on existing loci,
        where Lane variants are those of the type first brought to my attention
        by a paper by Dr. Lane. they are those where a variant in the context
        of a given transcript is just wildtype in a genomic reference. ie
        GRCh37/8

        Example:

        1 - new lanes - HOM loci:
        Generic middle cols = ... =
        ID  REF  ALT  QUAL  FILTER  INFO  GT:GQ:DP:AD:AF:PS  ./.:3,89:92:99:0.967:.

        1   25643553  ...   1:25643553_ref  loci

        2 - HET at loci:
        1  159175354   G   A  ... 1:159175354_G_A  1:159175354
        Becomes:
        1  159175354   G   A  ... 1:159175354_G_A,1:159175354_ref  1:159175354

        """

        new_lanes = {}

        for chrom, loci in self.lane_variants.items():
            # assert 'CHR' not in chrom.upper() #del after a few months of use TODO
            chrom = chrom.replace("chr", "")
            for pos in loci:
                # TODO blindly adding is problematic,  what if there's just no read
                # depth - vcf should be forced to report these
                lane_loci = f"{chrom}:{pos}"
                if lane_loci in self.loci:
                    GT = (
                        self.df.loc[self.df.loci == lane_loci, "SAMPLE"]
                        .values[0]
                        .split(":")[0]
                    )
                    assert GT.count("/") == 1 or GT.count("|") == 1
                    assert (
                        "2" not in GT
                    )  # these 2 asserts are designed to find examples so I can sort
                    # complex/multi var
                    if GT.startswith(("0/1", "0|1", "1/0", "1|0")):
                        self.df.loc[self.df.loci == lane_loci, "variant"] = (
                            self.df.loc[self.df.loci == lane_loci, "variant"].values[0]
                            + f",{lane_loci}_ref"
                        )
                else:
                    new_lanes[lane_loci] = (
                        [chrom, pos]
                        + COMMON_COLS[2:-1]
                        + ["GT:AD:GQ:DP:PS"]
                        + [
                            HOM_REF_DUMMY_QUAL,
                            f"{lane_loci}_ref",
                            "loci",
                        ]
                    )
        if new_lanes:
            new_lanes_df = pd.DataFrame.from_dict(new_lanes, orient="index")
            new_lanes_df.columns = self.df.columns
            self.df = pd.concat([self.df, new_lanes_df])

    def get_variants(self) -> dict[str, str]:
        """Retrieve variant information from the DataFrame.

        Returns:
            dict[str, str]: A dictionary of variants and their associated metrics.
        """
        vcf_variants = {}
        for varaint, metrics, format in zip(
            list(self.df.variant), list(self.df["SAMPLE"]), list(self.df["FORMAT"])
        ):
            if isinstance(metrics, float):
                continue
            assert format.startswith("GT")  # needed for add_lane_variants
            mapped_metrics = dict(
                zip(format.strip().split(":"), metrics.strip().split(":"))
            )
            mapped_metrics["GT"] == mapped_metrics["GT"].replace("|", "/")
            if mapped_metrics["GT"] == "0/0":
                continue
            if "," in varaint:
                for variant in varaint.split(","):
                    vcf_variants[variant] = mapped_metrics
            else:
                vcf_variants[varaint] = mapped_metrics

        return vcf_variants


def split_vcf_to_dfs(vcf_df: pd.DataFrame) -> pd.DataFrame:
    """Split multi-sample VCF DataFrame into individual sample DataFrames.

    Args:
        vcf_df (pd.DataFrame): Multi-sample VCF loaded into a DataFrame.

    Returns:
        Dict[str, pd.DataFrame]: Dictionary of sample-specific DataFrames.
    """
    # Extract column names related to samples
    sample_cols = [col for col in vcf_df.columns if col not in COMMON_COLS]

    for sample in sample_cols:
        try:
            assert all(row[1] in ("|", "/") for row in vcf_df[sample])
        except TypeError:
            logger.info(f"Sample {sample} is not diploid")
        cols: list[str] = COMMON_COLS + [sample]
        sample_vcf_df = vcf_df[cols].copy(deep=True)
        sample_vcf_df.columns = COMMON_COLS + ["SAMPLE"]
        yield sample_vcf_df, sample


def filter_VCF_to_BG_variants(df: pl.DataFrame, unique_variants) -> pd.DataFrame:
    """Filter a VCF represented as a Polars DataFrame to only include specified variants.

    This function creates a temporary column 'LOCI' by concatenating the 'CHROM' and
    'POS' columns, filters the DataFrame to retain only rows where 'LOCI' is in the
    provided unique_variants list, converts the result to a Pandas DataFrame, and
    removes the temporary 'LOCI' column.

    Args:
        df (pl.DataFrame): A Polars DataFrame containing VCF file data with columns
            such as "CHROM" and "POS".
        unique_variants (list[str]): A list of unique variant identifiers (e.g.,
            "chr:pos") to filter the DataFrame.

    Returns:
        pd.DataFrame: A Pandas DataFrame containing only the filtered variants from
            the original DataFrame, with the temporary 'LOCI' column removed.
    """
    df = df.with_columns(
        pl.concat_str(pl.col("CHROM"), pl.lit(":"), pl.col("POS")).alias("LOCI")
    )
    filtered_df = df.filter(pl.col("LOCI").is_in(unique_variants))
    if filtered_df.height == 0:  # empty
        pandas_df = df.to_pandas(use_pyarrow_extension_array=False)
    else:
        pandas_df = filtered_df.to_pandas(use_pyarrow_extension_array=False)

    del pandas_df["LOCI"]

    return pandas_df


class VcfMissingHeaderError(Exception):
    """
    Custom exception raised when a VCF file's header is missing,
    empty, or critically invalid (e.g., missing ##fileformat or #CHROM line).
    """

    def __init__(
        self, filename=None, message="VCF header is missing or invalid", reason=None
    ):
        """
        Initializes the VcfMissingHeaderError exception.

        Args:
            filename (str, optional): The path or name of the VCF file. Defaults to None.
            message (str, optional): The base error message.
                                     Defaults to "VCF header is missing or invalid".
            reason (str, optional): A specific reason for the header failure
                                    (e.g., "File is empty",
                                     "Missing '##fileformat' line",
                                     "Missing '#CHROM' line"). Defaults to None.
        """
        self.filename = filename
        self.base_message = message
        self.reason = reason

        full_message = self.base_message
        if filename:
            base_filename = os.path.basename(filename)
            full_message += f" in file: '{base_filename}'"
        if reason:
            full_message += f". Reason: {reason}"

        super().__init__(full_message)


class VcfNoDataError(Exception):
    """
    Custom exception raised when a VCF file is found to contain no
    variant data records (potentially only a header).
    """

    def __init__(self, filename=None, message="VCF contains no data records"):
        """
        Initializes the VcfNoDataError exception.

        Args:
            filename (str, optional): The path or name of the VCF file. Defaults to None.
            message (str, optional): The base error message.
                                     Defaults to "VCF contains no data records".
        """
        self.filename = filename
        self.message = message

        if filename:
            base_filename = os.path.basename(filename)
            full_message = f"{message} in file: '{base_filename}'"
        else:
            full_message = message

        super().__init__(full_message)

    def __str__(self):
        return super().__str__()


def read_vcf(file_path: str) -> pl.DataFrame:
    """Read a VCF file using polars while preserving the header and sample names.

    This function manually extracts the header (line starting with "#CHROM")
    and skips meta-information lines (starting with "##"). It then constructs a
    CSV-formatted string and parses it with polars.

    Args:
        file_path (str): Path to the VCF file (can be gzipped).

    Returns:
        pl.DataFrame: DataFrame containing the VCF data.
    """
    header = None
    # Use gzip.open if file is gzipped, else standard open.
    open_func = gzip.open if str(file_path).endswith(".gz") else open
    with open_func(file_path, "rt") as f:
        # Find header line starting with "#CHROM"
        for line in f:
            if line.startswith("##"):
                continue
            if line.startswith("#"):
                header = line.lstrip("#").strip().split("\t")
                if len(header) == 10:
                    header[-1] = "SAMPLE"  # for single sample
                break

        if header is None:
            raise VcfMissingHeaderError(filename=file_path)

        header_line = "\t".join(header) + "\n"
        data = f.read()

    csv_content = header_line + data
    df = pl.read_csv(
        io.StringIO(csv_content),
        separator="\t",
        schema_overrides=dict.fromkeys(["CHROM", "POS", "QUAL"], str),
    )
    if df.is_empty():
        raise VcfNoDataError(filename=file_path)
    df = df.with_columns(df["CHROM"].str.replace("chr", "", literal=True))
    return df


def check_if_multi_sample_vcf(file_path: str) -> bool:
    """Read a VCF file header.

    This function manually extracts the header (line starting with "#CHROM")
    to check if multi sample

    Args:
        file_path (str): Path to the VCF file (can be gzipped).

    Returns:
        bool: True if there's multiple samples

    """
    header = None
    # Use gzip.open if file is gzipped, else standard open.
    open_func = gzip.open if str(file_path).endswith(".gz") else open
    with open_func(file_path, "rt") as f:
        # Find header line starting with "#CHROM"
        for line in f:
            if line.startswith("##"):
                continue
            if line.startswith("#"):
                header = line.lstrip("#").strip().split("\t")
                if len(header) == 10:
                    return False
                elif len(header) < 10:
                    raise VcfMissingHeaderError(filename=file_path)
                else:
                    assert len(header) == len(set(header))
                    # unique sample names
            break

    return True
