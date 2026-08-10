import gzip
import io
import os
import re
from dataclasses import dataclass, field
from typing import Any
import pandas as pd
from icecream import ic
os.environ["POLARS_MAX_THREADS"] = "1"  # Must be set before polars import
import polars as pl
from loguru import logger
from collections import defaultdict
from rbceq2.core_logic.constants import (
    COMMON_COLS,
    HOM_REF_DUMMY_QUAL,
    HOM_REF_GTS,
    LANE,
    PAR,
)
from rbceq2.IO.encoders import VariantEncoderFactory


def gt_of(sample_field: str) -> str:
    """Pull the GT out of a SAMPLE column value.

    GT is always the first field - get_variants asserts FORMAT starts with 'GT'
    (vcf.py) and add_lane_variants relies on the same thing.

    Args:
        sample_field (str): A raw SAMPLE column value, ie '0/1:41,47:88:99' or, on an
            array export whose FORMAT is just 'GT', '0/1'.

    Returns:
        str: The GT, or '' if there is nothing parseable there.
    """
    if not isinstance(sample_field, str):
        return ""
    return sample_field.split(":")[0].strip()


def is_haploid_gt(GT: str) -> bool:
    """True if a GT names one allele rather than two, ie '1', '0' or '.'.

    Deliberately a string test and nothing more. Whether a haploid GT is *legitimate*
    depends on the coordinate, and that is is_single_copy's job - a haploid GT inside
    PAR is a caller error, not a ploidy statement.
    """
    return bool(GT) and "/" not in GT and "|" not in GT


def recode_gt_for_alt_index(GT: str, alt_index: int) -> str:
    """Rewrite a genotype so it describes one alternate allele of a multi-allelic row.

    A row with two alternates carries one genotype describing both, ie '1/2' for a sample
    heterozygous for each. The encoder already splits such a row into one token per
    alternate, but a token is a statement about *its* alternate only, so it needs a
    genotype that says how many copies of that one the sample has. This is that rewrite:
    the named index becomes 1, every other index becomes 0, '.' is left alone, and the
    separator is preserved so a phased genotype stays phased.

        '1/2' index 1 -> '1/0'      '1/2' index 2 -> '0/1'
        '2/2' index 1 -> '0/0'      '2/2' index 2 -> '1/1'
        '1|2' index 1 -> '1|0'      '1|2' index 2 -> '0|1'
        './2' index 2 -> './1'      '2'   index 2 -> '1'

    This is what makes the rest of the pipeline work unchanged rather than needing to
    learn about multi-allelic sites. get_ref sees only 0 and 1 again, so its existing
    zygosity logic is correct as written; and the phased filters compare phase strings
    literally - `phase == "1|0"` - so recoding is also what stops two different
    alternates at one position from both carrying '1|2' and comparing as though they
    were in phase with each other.

    Recoding to hom ref is meaningful and is the caller's to act on: '2/2' read for
    index 1 gives '0/0', which says the sample has zero copies of that alternate. See
    get_variants, which drops the token rather than storing it, the same rule
    remove_home_ref applies to a row.

    Deliberately not applied to a row with a single alternate. There the only valid
    index is 1, so the rewrite would be a no-op for every well formed genotype and would
    silently turn a malformed one - '2/2' where there is no second alternate - into
    '0/0' and drop it. That is input worth raising about, so it is left to reach get_ref.

    Args:
        GT (str): The genotype as written, ie '1/2', '1|2', '2'.
        alt_index (int): Which alternate this token is, 1 based, matching the order the
            alternates appear in the ALT column.

    Returns:
        str: The genotype rewritten for that alternate.
    """
    if not GT:
        return GT
    separator = "|" if "|" in GT else "/"
    recoded = [
        call if call == "." else ("1" if call == str(alt_index) else "0")
        for call in re.split(r"[/|]", GT)
    ]
    return separator.join(recoded)


def is_single_copy(chrom: str, pos: int, reference_genome: str) -> bool:
    """True if this coordinate sits on a chromosome that a male carries once.

    Non-PAR X and non-PAR Y. Everything else - every autosome, and X/Y inside PAR - is
    two copies in every sample, so this returns False for them and a haploid GT there is
    an input problem rather than constitutional haploidy.

    Note this asks about the *coordinate*, not the sample. A female is diploid at these
    coordinates too; whether a given sample is single copy is decided by
    VCF._infer_haploid_chroms, which needs both this and positive evidence from the GTs.

    Args:
        chrom (str): Chromosome with the 'chr' prefix already stripped, ie 'X'.
        pos (int): 1 based position.
        reference_genome (str): 'GRCh37' or 'GRCh38'.

    Returns:
        bool: True if the coordinate is outside PAR on X or Y.

    Example:
        >>> is_single_copy("X", 37686068, "GRCh38")   # XK
        True
        >>> is_single_copy("X", 2748343, "GRCh38")    # XG, inside PAR1
        False
        >>> is_single_copy("1", 25272548, "GRCh38")   # RHD, autosomal
        False
    """
    par_for_build = PAR.get(reference_genome)
    if par_for_build is None:
        return False
    intervals = par_for_build.get(chrom)
    if intervals is None:
        return False
    return not any(start <= pos <= end for start, end in intervals)


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
        reference_genome (str | None):
            'GRCh37' or 'GRCh38'. Optional, and None means ploidy inference is off and
            every region is treated as two copies - the behaviour before v2.4.4. Kept
            optional so the many test constructions of VCF keep working unchanged.
        haploid_chroms (frozenset[str]):
            Chromosomes this sample's caller reported as single copy outside PAR.
        haploid_loci (dict[str, frozenset[int]]):
            Chromosome -> positions this sample's caller wrote a haploid GT at.
        diploid_loci (dict[str, frozenset[int]]):
            Chromosome -> positions this sample's caller wrote a diploid GT at.
    """

    input_vcf: pl.DataFrame | pd.DataFrame
    lane_variants: dict[str, Any]
    unique_variants: set[str]
    sample: str  # field(init=False)
    reference_genome: str | None = None
    df: pd.DataFrame = field(init=False)
    loci: set[str] = field(init=False)
    variants: dict[str, str] = field(init=False)
    phase_sets: dict[str, dict[int, tuple[int, int]]] = field(init=False)
    haploid_chroms: frozenset[str] = field(init=False)
    haploid_loci: dict[str, frozenset[int]] = field(init=False)
    diploid_loci: dict[str, frozenset[int]] = field(init=False)

    def __post_init__(self):
        """Handle initialization after data class creation."""
        object.__setattr__(self, "df", self.handle_single_or_multi())
        # object.__setattr__(self, "sample", self.get_sample())
        self.rename_chrom()
        # Has to run before remove_home_ref: whether a haploid '0' is a hom ref call
        # depends on the answer, and remove_home_ref is where that row gets dropped.
        object.__setattr__(self, "haploid_chroms", self._infer_haploid_chroms())
        haploid_loci, diploid_loci = self._infer_locus_ploidy()
        object.__setattr__(self, "haploid_loci", haploid_loci)
        object.__setattr__(self, "diploid_loci", diploid_loci)
        self.remove_home_ref()
        self.encode_variants()
        self.add_loci()
        object.__setattr__(self, "loci", self.set_loci())
        self.add_lane_variants()
        object.__setattr__(self, "variants", self.get_variants())
        self._create_phase_sets()

    def _create_phase_sets(self) -> None:
        """Parses the VCF DataFrame to find phased variants and stores their
        chromosomal ranges (min and max position) by phase set ID.

        This method populates the `self.phase_sets` attribute with a nested
        dictionary: {chromosome: {phase_set_id: (min_position, max_position)}}.
        """
        # Temporary dict to aggregate all positions for each phase set
        # Format: {chrom: {ps_id: [pos1, pos2, ...]}}
        temp_phase_data = {}

        # Filter for rows that might contain phasing info to reduce work
        df_phased = self.df[self.df["FORMAT"].str.contains("PS", na=False)].copy()

        # Convert POS to integer once for performance
        df_phased["POS"] = pd.to_numeric(df_phased["POS"])

        for _, row in df_phased.iterrows():
            format_keys = row["FORMAT"].split(":")
            try:
                ps_index = format_keys.index("PS")
            except ValueError:
                continue  # 'PS' not in this row's FORMAT string

            sample_values = row["SAMPLE"].split(":")
            genotype = sample_values[0]

            # Ensure genotype is phased ('|') and the PS value is valid
            if "|" in genotype and len(sample_values) > ps_index:
                ps_value = sample_values[ps_index]
                if ps_value != ".":
                    chrom = row["CHROM"]
                    pos = row["POS"]
                    ps_id = int(ps_value)

                    # Initialize nested dicts if they don't exist and append pos
                    temp_phase_data.setdefault(chrom, {}).setdefault(ps_id, []).append(
                        pos
                    )

        # Convert the lists of positions to (min, max) tuples
        final_phase_sets = {}
        for chrom, ps_groups in temp_phase_data.items():
            final_phase_sets[chrom] = {
                ps_id: (min(positions), max(positions))
                for ps_id, positions in ps_groups.items()
            }

        object.__setattr__(self, "phase_sets", final_phase_sets)

    def handle_single_or_multi(self) -> pd.DataFrame:
        """Handle single or multiple entries in the VCF, returning a DataFrame.

        Returns:
            pd.DataFrame: The DataFrame representation of the VCF data.
        """

        if isinstance(self.input_vcf, pl.DataFrame):
            return filter_VCF_to_BG_variants(self.input_vcf, self.unique_variants)
        else:
            return self.input_vcf[0]

    def rename_chrom(self) -> None:
        """Rename chromosome identifiers by removing the 'chr' prefix."""
        self.df["CHROM"] = self.df["CHROM"].apply(lambda x: x.replace("chr", ""))

    def _infer_haploid_chroms(self) -> frozenset[str]:
        """Chromosomes this sample's caller reported as single copy outside PAR.

        Positive evidence only: at least one haploid GT at a non-PAR X/Y coordinate. That
        is the one signal a VCF actually carries. Sex is not in a VCF, so a caller that
        codes a male as '1/1' across non-PAR X leaves nothing to find and the sample stays
        diploid - that is state table row B3, and it is a known, permanent limitation
        rather than an oversight. B3 gets the phenotype right and only the genotype string
        wrong.

        Deliberately per sample, not per file. A multi-sample array export mixes both
        codings in one file - males haploid across non-PAR X, females diploid - and each
        VCF object here already holds one sample.

        Deliberately restricted to X and Y. A haploid GT on an autosome is a copy number
        statement about a *locus*, not about how many chromosomes the sample has, and the
        two have different correct answers - see get_ref. Returning early when
        reference_genome is None keeps every existing caller on the old behaviour.

        Returns:
            frozenset[str]: ie frozenset({'X'}), or an empty set for a female sample, an
            autosome-only VCF, or when reference_genome was not supplied.
        """
        if self.reference_genome is None:
            return frozenset()

        haploid: set[str] = set()
        for chrom, pos, sample_field in zip(
            self.df["CHROM"], self.df["POS"], self.df["SAMPLE"]
        ):
            if chrom in haploid or chrom not in PAR.get(self.reference_genome, {}):
                continue
            if not is_haploid_gt(gt_of(sample_field)):
                continue
            try:
                position = int(pos)
            except (TypeError, ValueError):
                continue
            if is_single_copy(chrom, position, self.reference_genome):
                haploid.add(chrom)

        if haploid:
            logger.info(
                f"{self.sample}: haploid GTs found outside PAR on "
                f"{','.join(sorted(haploid))}, so those regions are reported with one "
                f"allele slot"
            )
        return frozenset(haploid)

    def _infer_locus_ploidy(
        self,
    ) -> tuple[dict[str, frozenset[int]], dict[str, frozenset[int]]]:
        """Which loci this sample's caller wrote haploid GTs at, and which diploid.

        Raw observation, no interpretation. What a haploid GT *means* depends on where it
        is - constitutional haploidy outside PAR, gene copy number on an autosome, a caller
        error inside PAR - and answering that needs the database's gene boundaries, which
        this layer does not have. So the reading is done here and the meaning in
        locus_copies_for_bg, which has both halves.

        Has to run before remove_home_ref, and this is the whole reason it is a separate
        pass rather than something read off the variant pool later. A gene called entirely
        wildtype has every row dropped, so by pool time 'one copy, and it is reference' and
        'nobody typed this gene' are the same empty set. Here they are still distinct.

        Deliberately records both halves rather than just the haploid one. The question
        downstream is whether a gene is *consistently* single copy, and a gene with some
        haploid and some diploid loci is a contradiction that must not be read as copy
        number - knowing only where the haploid GTs were could not tell the two apart.

        Returns:
            tuple[dict[str, frozenset[int]], dict[str, frozenset[int]]]: haploid loci and
            diploid loci, each chromosome -> positions.
        """
        haploid: dict[str, set[int]] = defaultdict(set)
        diploid: dict[str, set[int]] = defaultdict(set)

        for chrom, pos, sample_field in zip(
            self.df["CHROM"], self.df["POS"], self.df["SAMPLE"]
        ):
            GT = gt_of(sample_field)
            if not GT:
                continue
            try:
                position = int(pos)
            except (TypeError, ValueError):
                continue
            target = haploid if is_haploid_gt(GT) else diploid
            target[str(chrom)].add(position)

        return (
            {chrom: frozenset(pos) for chrom, pos in haploid.items()},
            {chrom: frozenset(pos) for chrom, pos in diploid.items()},
        )

    def _is_haploid_hom_ref(self, sample_field: Any) -> bool:
        """True for a haploid '0', anywhere.

        The single-copy equivalent of '0/0'. Split out of remove_home_ref so the row-wise
        test stays readable next to the vectorised startswith it sits beside.

        Deliberately asks nothing about the coordinate, unlike everything else here that
        touches ploidy. Whether a haploid GT means one chromosome, one copy of a gene, or a
        caller having a bad day is undecided at this point and needs the database to
        settle. It does not need settling: under every one of those readings the caller has
        said the copies it can see are reference, so the ALT token has zero copies and
        absence is the right encoding. This is the one ploidy question that can be answered
        without knowing which question it is.

        Dropping the row does not lose the ploidy evidence. _infer_locus_ploidy has already
        recorded that this locus was haploid, which is what a mixed coding is detected from
        later.
        """
        return gt_of(sample_field) == "0"

    def remove_home_ref(self) -> None:
        """Remove homozygous reference calls from the DataFrame.

        A hom ref call carries no allele, so the row is dropped and the token's absence
        from the variant pool is what records 'zero copies of this ALT'. Kept as a row it
        would instead become an ALT token asserting the sample carries the variant on both
        chromosomes - the exact inverse of what the GT says.

        Both separators are matched. Only '0/0' was matched until v2.4.3, so a phased
        '0|0' survived and produced that inverted HOM ALT token. No caller in any test
        dataset emits '0|0' (they all write hom ref unphased even in phased VCFs), so the
        bug was latent rather than observed.

        Haploid '0' is matched too, and since v2.4.5 anywhere rather than only outside PAR
        on X/Y. It means exactly what '0/0' means wherever it appears: the copies the
        caller can see are reference, so the token has zero copies. Until v2.4.4 an
        autosomal one was left in place to be rejected by name in get_ref, which was a
        raise about a row carrying no allele - state table row D3. The reading of what the
        haploidy *means* still happens later and elsewhere; see _is_haploid_hom_ref.
        """
        hom_ref = self.df["SAMPLE"].str.startswith(HOM_REF_GTS, na=False)

        haploid_hom_ref = pd.Series(
            [self._is_haploid_hom_ref(sample_field) for sample_field in self.df["SAMPLE"]],
            index=self.df.index,
            dtype=bool,
        )

        self.df = self.df[~(hom_ref | haploid_hom_ref)].copy(deep=True)

    def encode_variants(self) -> None:
        """Encode variants into a unified format using the encoder factory."""
        factory = VariantEncoderFactory()

        self.df["variant"] = self.df.apply(
            lambda row: factory.encode_variant(row), axis=1
        )

    def add_loci(self) -> None:
        """Add loci identifiers to the DataFrame."""
        self.df["loci"] = self.df.CHROM + ":" + self.df.POS


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

        Added variant details (LANE constant) in v2.3.4 after seeing instances
        where 133257521 had a SNP, not ins,
        so lane was created erroneously

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
        new_rows = []
        for chrom, loci in self.lane_variants.items():
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
                    # try:
                    #     assert GT.count("/") == 1 or GT.count("|") == 1
                    #     assert ("2" not in GT) TODO, not sure if this is too defensive...
                    # except AssertionError:
                    #     print('multi allele loci are not supported, please use bcftools norm -m -both ...')
                    #     raise ValueError('multi allele loci are not supported, please use bcftools norm -m -both on your VCF/s')
                    if (
                        GT.startswith(("0/1", "0|1", "1/0", "1|0"))
                        and len(self.df.loc[self.df.loci == lane_loci, "SAMPLE"]) == 1
                    ):
                        # HET and not multi allelic = ref
                        ref = self.df.loc[self.df.loci == lane_loci, "REF"].values[0]
                        alt = self.df.loc[self.df.loci == lane_loci, "ALT"].values[0]
                        lane = LANE[f"chr{chrom}"][pos]
                        if f"{ref}_{alt}" == lane or lane == "no_ALT":
                            original_row = (
                                self.df.loc[self.df.loci == lane_loci].iloc[0].copy()
                            )

                            # Create the reference variant row
                            ref_row = original_row.copy()
                            ref_row["variant"] = f"{lane_loci}_ref"
                            ref_row["ALT"] = original_row[
                                "REF"
                            ]  # ALT becomes the original REF

                            # Flip the genotype in FORMAT column
                            gt_field = ref_row["SAMPLE"].split(":")[0]
                            if "|" in gt_field:
                                # Phased: flip 0|1 to 1|0 or 1|0 to 0|1
                                flipped_gt = "|".join(reversed(gt_field.split("|")))
                            elif "/" in gt_field:
                                # Unphased: flip 0/1 to 1/0 or 1/0 to 0/1
                                flipped_gt = "/".join(reversed(gt_field.split("/")))
                            else:
                                raise ValueError("GT formated wrong")

                            # Replace the GT field in SAMPLE
                            sample_fields = ref_row["SAMPLE"].split(":")
                            sample_fields[0] = flipped_gt
                            ref_row["SAMPLE"] = ":".join(sample_fields)
                            self.df.loc[self.df.loci == lane_loci, "variant"] = (
                                f"{lane_loci}_{lane}"
                            )
                            new_rows.append(ref_row)
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
        if new_rows:
            self.df = pd.concat([self.df, pd.DataFrame(new_rows)], ignore_index=True)
        if new_lanes:
            new_lanes_df = pd.DataFrame.from_dict(new_lanes, orient="index")
            new_lanes_df.columns = self.df.columns
            self.df = pd.concat([self.df, new_lanes_df])

    def get_variants(self) -> dict[str, str]:
        """Retrieve variant information from the DataFrame.

        The hom ref skip here is the second of two, and is a safety net rather than the
        primary filter: remove_home_ref has already dropped these rows before the loci set
        is built. Both are needed and they are not interchangeable - dropping the row in
        remove_home_ref also removes the locus from self.loci, which is what lets
        add_lane_variants synthesise the '_ref' partner for a lane locus. Skipping only
        here would leave the locus looking called and no '_ref' token would be made.

        This is also where a multi-allelic row becomes usable. The encoder has already
        emitted one token per alternate, in ALT column order, comma joined - so the split
        below is not a string convenience, it is the alternates. What each token was
        missing is a genotype of its own: until v2.4.6 they all took the row's genotype
        unchanged, so both tokens of a '1/2' were handed '1/2' and get_ref refused it as
        multi-allelic. Each now gets its own recoded genotype and its own metrics dict.

        The dicts were shared before, which nothing depended on but which meant a write
        through one token's metrics would have reached the other's.

        A token recoding to hom ref is dropped rather than stored, which is the row level
        rule of remove_home_ref applied per alternate: '2/2' read for the first alternate
        is '0/0', the sample has none of it, and absence from the pool is how zero copies
        is encoded.

        Rows with one alternate keep the old path exactly, sharing no code with the split
        one - see recode_gt_for_alt_index for why rewriting them would be worse than
        leaving them alone.

        Returns:
            dict[str, str]: A dictionary of variants and their associated metrics.
        """
        vcf_variants = {}
        for variant, metrics, format in zip(
            list(self.df.variant), list(self.df["SAMPLE"]), list(self.df["FORMAT"])
        ):
            if isinstance(metrics, float):
                continue
            assert format.startswith("GT")  # needed for add_lane_variants
            mapped_metrics = dict(
                zip(format.strip().split(":"), metrics.strip().split(":"))
            )
            if mapped_metrics["GT"] in HOM_REF_GTS:
                continue
            if "," in variant:
                for alt_index, alt_variant in enumerate(variant.split(","), start=1):
                    per_alt_metrics = dict(mapped_metrics)
                    per_alt_metrics["GT"] = recode_gt_for_alt_index(
                        mapped_metrics["GT"], alt_index
                    )
                    if per_alt_metrics["GT"] in HOM_REF_GTS or (
                        per_alt_metrics["GT"] == "0"
                    ):
                        continue
                    vcf_variants[alt_variant] = per_alt_metrics
            else:
                vcf_variants[variant] = mapped_metrics

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
        # Informational, not validation - a non-diploid GT is reported and kept, and it is
        # get_ref that decides whether it is a legitimate ploidy statement. This was a bare
        # assert guarded by 'except TypeError', which meant it did the intended thing for a
        # NaN cell but raised IndexError on the case it names: a haploid GT like '0' has no
        # index 1
        gts = [gt_of(row) for row in vcf_df[sample]]
        if any(not is_haploid_gt(GT) and len(GT) < 2 for GT in gts):
            logger.info(f"Sample {sample} has GTs that are neither haploid nor diploid")
        elif any(is_haploid_gt(GT) for GT in gts):
            logger.info(f"Sample {sample} is not diploid at every locus")
        cols: list[str] = COMMON_COLS + [sample]
        sample_vcf_df = vcf_df[cols].copy(deep=True)
        sample_vcf_df.columns = COMMON_COLS + ["SAMPLE"]
        yield sample_vcf_df, sample


def find_phased_neighbors(df: pd.DataFrame) -> set[str]:
    """
    rescues ABOs c.261delG - indels don't always get phased properly"""
    central_loci_to_find = {
        "9:133257521": 133257521,
        "9:136132908": 136132908,
    }

    # 1. Create the sorted list of all PHASED loci on Chromosome 9
    phased_on_chrom9 = (
        df.filter(pl.col("CHROM") == "9")
        .with_columns(pl.col("POS").cast(pl.Int64))
        .filter(pl.col("FORMAT").str.contains("PS"))
        .sort("POS")
    )

    # If there are no phased loci at all, exit
    if phased_on_chrom9.height == 0:
        return set()

    # Extract the positions and loci as separate series for quick lookups
    phased_positions = phased_on_chrom9.get_column("POS")
    phased_loci_series = phased_on_chrom9.get_column("LOCI")

    results = []
    # 2. For each central locus, find its place in the sorted list
    for locus_id, locus_pos in central_loci_to_find.items():
        # search_sorted finds the index where `locus_pos` would be inserted
        # to maintain the sort order. This is the index of the first variant
        # at or after our central locus.
        idx = phased_positions.search_sorted(locus_pos)

        # 3. Use this index to find the neighbors from our list of phased loci
        # We need to be careful about edges (e.g., asking for index -2 when idx is 0 or 1)
        results.append(
            {
                "LOCI": locus_id,
                "prev_2": phased_loci_series[idx - 2] if idx > 1 else None,
                "prev_1": phased_loci_series[idx - 1] if idx > 0 else None,
                "next_1": phased_loci_series[idx]
                if idx < len(phased_loci_series)
                else None,
                "next_2": phased_loci_series[idx + 1]
                if idx < len(phased_loci_series) - 1
                else None,
            }
        )
    neighbor_cols = ["prev_2", "prev_1", "next_1", "next_2"]
    neighbours_df = pl.from_dicts(results)
    # 4. Convert the list of dictionaries to a final DataFrame
    unique_loci_set = {
        locus
        for row in neighbours_df.select(neighbor_cols).rows()
        for locus in row
        if locus is not None  # Filter out the nulls
    }
    return unique_loci_set


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
    # TODO maybe best to switch to tabix?
    # although fuzzy mtaching won't work with tabix...
    df = df.with_columns(
        df["CHROM"].str.replace("chr", "", literal=True).alias("CHROM")
    )
    df = df.with_columns(
        pl.concat_str(pl.col("CHROM"), pl.lit(":"), pl.col("POS")).alias("LOCI")
    )
    large_vars = set(
        df.filter((df["REF"].str.len_chars() > 50) | (df["ALT"].str.len_chars() > 50))[
            "LOCI"
        ]
    )
    massive_vars = set(
        df.filter((df["ALT"].str.contains("<")) | (df["ALT"].str.contains(">")))["LOCI"]
    )
    neighbours = find_phased_neighbors(df)
    merged_set = neighbours | unique_variants | large_vars | massive_vars
    filtered_df = df.filter(pl.col("LOCI").is_in(merged_set))
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


@dataclass(frozen=True, slots=True)
class Interval:
    start: int
    end: int


def parse_positions(db_col: str) -> list[int]:
    """Extract numeric positions from a database column entry string."""
    if pd.isna(db_col):
        return []
    positions = []
    for tok in str(db_col).split(","):
        if "_" in tok:
            pos = tok.split("_", 1)[0]
            if pos.isdigit():
                positions.append(int(pos))
    return positions


def build_intervals(
    db: pd.DataFrame, genome: str, flank: int = 500_000
) -> dict[str, list[Interval]]:
    """Construct per-chrom intervals ±flank from database DataFrame.

    Args:
        db (pd.DataFrame): DataFrame with at least 'Chrom' and genome columns.
        genome (str): Column to use ('GRCh37' or 'GRCh38').
        flank (int): Window size on either side of each variant.

    Returns:
        dict[str, list[Interval]]: Per-chromosome merged intervals.
    """
    intervals = defaultdict(list)

    for row in db.itertuples(index=False):
        chrom = getattr(row, "Chrom").removeprefix("chr")
        genome_col = getattr(row, genome)
        for pos in parse_positions(genome_col):
            intervals[chrom].append(Interval(max(0, pos - flank), pos + flank))

    # merge overlapping intervals per chromosome
    merged = {}
    for chrom, ivals in intervals.items():
        ivals.sort(key=lambda x: x.start)
        merged_list: list[Interval] = []
        for iv in ivals:
            if not merged_list or iv.start > merged_list[-1].end:
                merged_list.append(iv)
            else:
                merged_list[-1] = Interval(
                    merged_list[-1].start, max(merged_list[-1].end, iv.end)
                )
        merged[chrom] = merged_list
    return merged


def variant_in_intervals(
    chrom: str, pos: int, intervals: dict[str, list[Interval]]
) -> bool:
    """Check if a variant lies in any interval for that chrom."""
    if chrom not in intervals:
        return False
    for iv in intervals[chrom]:
        if iv.start <= pos <= iv.end:
            return True
    return False


def row_is_wanted(
    chrom: str,
    pos: int,
    ref: str,
    alt: str,
    fmt: str,
    unique_variants: set[str],
) -> bool:
    """Will filter_VCF_to_BG_variants keep this row? Answered from the raw line.

    Everything that function decides *after* parsing is decidable *before* it, and doing it
    here rather than there is the difference between parsing thousands of rows and parsing
    hundreds of thousands. The set kept is a superset of what it goes on to select, so the
    end result is identical - this only stops rows being parsed to be thrown away.

    Why it matters: read_vcf's interval test keeps anything within 500kb of a database
    position, which is sparse for a genotyping array but dense for a jointly called cohort.
    One 3,209 sample cohort parsed 379,547 rows x 3,218 columns to keep 459, which is tens
    of GB in one process before the Pool is even created - so --processes could not affect
    it. With this it parses 1,868.

    The four clauses mirror filter_VCF_to_BG_variants exactly:

    - at a database locus                      -> its unique_variants set
    - REF or ALT longer than 50 characters     -> its large_vars
    - ALT symbolic, ie '<DEL>'                 -> its massive_vars
    - chr9 and phased                          -> its find_phased_neighbors

    The last is the awkward one and the reason this is not simply 'is it a database locus'.
    find_phased_neighbors rescues ABO c.261delG by taking the two phased loci either side of
    two fixed positions, and how far away those are depends on how sparse phasing is there.
    Windowing chr9 by position was measured against ten phased files and even 100kb picked a
    different set on one of them, so the test is on the PS field rather than on distance:
    every row that function could possibly select is kept, and no guess is involved.

    Args:
        chrom (str): Chromosome, 'chr' prefix already stripped.
        pos (int): 1 based position.
        ref (str): REF field.
        alt (str): ALT field.
        fmt (str): FORMAT field.
        unique_variants (set[str]): Database loci as 'chrom:pos'.

    Returns:
        bool: True if the row can still be needed once the frame is parsed.
    """
    return (
        f"{chrom}:{pos}" in unique_variants
        or len(ref) > 50
        or len(alt) > 50
        or "<" in alt
        or ">" in alt
        or (chrom == "9" and "PS" in fmt)
    )


def read_vcf(
    vcf_path: str,
    intervals: dict[str, list[Interval]],
    unique_variants: set[str] | None = None,
) -> pl.DataFrame:
    """Stream a VCF, keep only relevant lines, return as Polars DataFrame.
    read a VCF file using polars while preserving the header and sample names.

    This function manually extracts the header (line starting with "#CHROM")
    and skips meta-information lines (starting with "##"). It then constructs a
    CSV-formatted string and parses it with polars.

    Two filters, not one, and the second is optional. The interval test is a coarse bound -
    anything within 500kb of a database position - and on its own it is only cheap when the
    caller wrote sparsely. Supplying unique_variants adds the fine test, row_is_wanted,
    which is the same decision filter_VCF_to_BG_variants makes once the frame exists. It
    changes nothing about the result and a great deal about how much is parsed to get there;
    see row_is_wanted for the measurements.

    Optional so that every existing caller and test keeps the old behaviour untouched. Omit
    it and the interval test is all that applies, exactly as before.

    Args:
        vcf_path (str): Path to the VCF file (can be gzipped).
        intervals (dict[str, list[Interval]]): Per chromosome windows around database
            positions, from build_intervals.
        unique_variants (set[str] | None): Database loci as 'chrom:pos'. None skips the
            fine filter and keeps every row inside the intervals.

    Returns:
        pl.DataFrame: DataFrame containing the VCF data."""

    open_func = gzip.open if vcf_path.endswith(".gz") else open
    header = None
    rows: list[str] = []
    with open_func(vcf_path, "rt") as f:
        for line in f:  # TODO Pool
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                header = line.lstrip("#").strip().split("\t")
                if len(header) == 10:
                    header[-1] = "SAMPLE"  # for single sample
                continue

            # parse variant
            fields = line.split("\t")
            try:
                chrom, pos = fields[0].removeprefix("chr"), int(fields[1])
            except:
                raise
            if not variant_in_intervals(chrom, pos, intervals):
                continue
            if unique_variants is not None and not row_is_wanted(
                chrom,
                pos,
                fields[3],
                fields[4],
                # A sites-only VCF has no FORMAT column. Nothing is phased in one either,
                # so an empty string is the right answer rather than an IndexError.
                fields[8] if len(fields) > 8 else "",
                unique_variants,
            ):
                continue
            rows.append(line)

    if header is None:
        raise VcfMissingHeaderError(filename=vcf_path)
    header_line = "\t".join(header) + "\n"
    csv_content = header_line + "".join(rows)
    try:
        df = pl.read_csv(
            io.StringIO(csv_content),
            separator="\t",
            # Every VCF field is text, and inferring types from the first rows gets it
            # wrong in a way that is silent rather than loud. A file whose leading rows
            # are haploid genotypes - which is what a caller encoding gene copy number
            # produces - infers the sample column as an integer, and then every ordinary
            # '0/0' further down the file is unparseable and lands as null. Reading
            # everything as text removes the whole class.
            infer_schema_length=0,
            schema_overrides={"CHROM": str, "POS": str, "QUAL": str},
        )
    except pl.exceptions.ComputeError:
        df = pl.read_csv(
            io.StringIO(csv_content),
            separator="\t",
            # Every VCF field is text, and inferring types from the first rows gets it
            # wrong in a way that is silent rather than loud. A file whose leading rows
            # are haploid genotypes - which is what a caller encoding gene copy number
            # produces - infers the sample column as an integer, and then every ordinary
            # '0/0' further down the file is unparseable and lands as null. Reading
            # everything as text removes the whole class.
            infer_schema_length=0,
            schema_overrides={"CHROM": str, "POS": str, "QUAL": str},
            truncate_ragged_lines=True,
        )
    except MemoryError:
        message = "VCF is too big, plz trim ie bcftools view -R regions.bed ..."
        print(message)
        logger.error(message)
        raise

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
