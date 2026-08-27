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
from rbceq2.core_logic.filter_semantics import FILTER_VALUE_SEPARATOR
from rbceq2.IO.encoders import VariantEncoderFactory

# How many contradicted tokens get named in the warning before it says "and N more".
# One line per sample either way - a jointly called cohort can contradict a dozen at
# once and the point of the line is that the user goes and looks, not that it lists
# them all.
MAX_CONTRADICTED_TOKENS_LOGGED = 4


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
    PAR is a caller error, not a ploidy statement. Whether it carries any information at
    all is gt_names_an_allele's, and '.' passes this test while failing that one.
    """
    return bool(GT) and "/" not in GT and "|" not in GT


def gt_names_an_allele(GT: str) -> bool:
    """True if a GT states at least one allele rather than declining to call.

    A no-call is not evidence of ploidy in either direction, and it has two spellings
    that this is the only thing separating from real calls. '.' is haploid by shape and
    './.' is diploid by shape, so without this test one of them argues for a single copy
    and the other against it - two ways of writing 'I do not know' landing on opposite
    sides of the question.

    That is not hypothetical. Six '.' no-calls at non-PAR X coordinates were enough to
    make _infer_haploid_chroms declare a sample single copy across X, which reported
    XK, ATP11C and GATA1 with one allele slot and lost a chromosome the sample has - the
    exact failure the positive-evidence rule was written to avoid, arriving through the
    one GT that looks like evidence and is not.

    Args:
        GT (str): A genotype string, ie '0/1', '1', '.', './.'.

    Returns:
        bool: True if any position in the genotype names an allele.

    Example:
        >>> gt_names_an_allele("0")
        True
        >>> gt_names_an_allele(".")
        False
        >>> gt_names_an_allele("./.")
        False
        >>> gt_names_an_allele("./1")
        True
    """
    return any(call not in ("", ".") for call in re.split(r"[/|]", GT))


def is_single_base_substitution(ref: str, alt: str, info: str) -> bool:
    """Is this row a plain substitution of one base for another?

    The test that decides which rows reconcile_duplicate_rows may touch, and it is
    deliberately narrower than "not structural". The structural reader takes any row
    whose INFO names a type, whose ALT is symbolic, or whose sequences differ by more
    than --min_size, and that threshold is a run time argument this layer never sees.
    One base on both sides sits below every threshold that argument could hold, so a
    row matching this can never be one the structural reader wanted.

    Conservative on purpose. Every duplicate small variant observed so far is a
    substitution, so nothing is lost by it, and widening it later needs a real example
    rather than an argument.

    Args:
        ref (str): The REF column.
        alt (str): The ALT column.
        info (str): The INFO column, checked for a structural type.

    Returns:
        bool: True where the row is one base substituted for another, carrying no
        structural claim.
    """
    return (
        len(str(ref)) == 1
        and len(str(alt)) == 1
        and "SVTYPE" not in str(info)
    )


def rows_make_the_same_call(gts: list[str]) -> bool:
    """Do these genotypes say the same thing about the alternate?

    Not string equality. Two callers can write the same call three ways and only one of
    the differences is a disagreement:

    - '0|1' and '0/1' - the same call, one of them phased. Agree.
    - '1/1' and '1|1' - likewise, though phase says nothing about a homozygote.
    - '1/1/1/1' and '1/1' - different ploidy, same claim: every copy carries it. Agree.
      A gene conversion caller reporting a paralogue pair writes the first, the general
      caller the second, and neither is wrong.
    - '0|1' and '1/1' - a real disagreement, and nothing here resolves it.

    So the test is: identical once the separator is normalised, or every copy reference
    on both sides, or every copy alternate on both sides. Ploidy is per genotype in a
    VCF, so differing ploidy is not by itself a disagreement.

    Args:
        gts (list[str]): The genotypes of the rows carrying one token, in file order.

    Returns:
        bool: True where every row makes the same claim about the alternate.
    """
    normalised = {gt.replace("|", "/") for gt in gts}
    if len(normalised) == 1:
        return True
    calls = [[c for c in gt.replace("|", "/").split("/")] for gt in gts]
    if any("." in call for call in calls):
        return False
    if all(set(call) == {"0"} for call in calls):
        return True

    return all(set(call) == {"1"} for call in calls)


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


@dataclass
class PloidyScan:
    """Collects evidence of constitutional haploidy while a VCF is being read.

    The evidence is there and then it is gone. read_vcf drops rows in two passes - a
    coarse interval test and then row_is_wanted, which keeps only what the frame will
    actually need - and a sample's haploid genotypes are usually not at a database locus,
    so row_is_wanted discards them. Measured on one per-sample file: 480 haploid called
    genotypes outside PAR, 338 still there after the interval test, **none** after
    row_is_wanted. The sample is male and was reported diploid on XK, ATP11C and GATA1.

    So this reads the evidence where it still exists rather than moving the inference.
    The alternative was to make row_is_wanted keep every non-PAR X/Y row, which on one
    cohort took the retained rows from 5,031 to 25,682 and parsed twenty million extra
    cells - re-spending most of what the read filter was added to save, to compute one
    boolean per sample.

    Deliberately a separate object rather than two more arguments to read_vcf. It carries
    the one input the scan needs and the evidence it produces, so read_vcf keeps four
    parameters and the whole thing can be tested without reading a file.

    Attributes:
        reference_genome (str): 'GRCh37' or 'GRCh38', needed to know where PAR is.
        haploid_chroms (dict[str, set[str]]): Sample name -> chromosomes that sample has
            a haploid called genotype on, outside PAR. Sample names are the VCF's own
            column names, which for a single sample file read_vcf has renamed to 'SAMPLE'.
    """

    reference_genome: str
    haploid_chroms: dict[str, set[str]] = field(default_factory=lambda: defaultdict(set))

    def observe(self, chrom: str, pos: int, samples: list[str], fields: list[str]) -> None:
        """Record any haploid called genotype at a coordinate that is single copy.

        Called for every row read, so the cheap tests come first: everything that is not
        X or Y leaves on the first line, which is the overwhelming majority of a file.

        A no-call is not evidence - see gt_names_an_allele - and neither is a haploid GT
        inside PAR, which is a caller error rather than a ploidy statement.

        Args:
            chrom (str): Chromosome, 'chr' prefix already stripped.
            pos (int): 1 based position.
            samples (list[str]): The VCF's sample column names, in column order.
            fields (list[str]): The row's tab separated fields.
        """
        if chrom not in ("X", "Y"):
            return
        if not is_single_copy(chrom, pos, self.reference_genome):
            return
        for offset, sample in enumerate(samples):
            if chrom in self.haploid_chroms.get(sample, ()):
                continue
            index = 9 + offset
            if index >= len(fields):
                break
            GT = gt_of(fields[index])
            if is_haploid_gt(GT) and gt_names_an_allele(GT):
                self.haploid_chroms[sample].add(chrom)

    def for_sample(self, sample: str) -> frozenset[str]:
        """What was seen for one sample, as the frozenset VCF expects."""
        return frozenset(self.haploid_chroms.get(sample, ()))


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
    observed_haploid_chroms: frozenset[str] | None = None
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
        # Union rather than preference. The scan sees rows the frame never will, and the
        # frame is what every caller that does not scan still relies on, so neither is
        # allowed to overrule the other - either one finding a haploid genotype outside
        # PAR is positive evidence, which is the whole basis of this inference.
        haploid_chroms = self._infer_haploid_chroms() | (
            self.observed_haploid_chroms or frozenset()
        )
        # Reported here rather than inside _infer_haploid_chroms, because it has to say
        # what was concluded rather than what one of the two sources found. A sparse file
        # usually has no haploid rows left by the time the frame exists, so the scan is
        # the only source and logging from there would go quiet on exactly the samples
        # this is meant to announce.
        if haploid_chroms:
            logger.info(
                f"{self.sample}: haploid GTs found outside PAR on "
                f"{','.join(sorted(haploid_chroms))}, so those regions are reported "
                f"with one allele slot"
            )
        object.__setattr__(self, "haploid_chroms", haploid_chroms)
        haploid_loci, diploid_loci = self._infer_locus_ploidy()
        object.__setattr__(self, "haploid_loci", haploid_loci)
        object.__setattr__(self, "diploid_loci", diploid_loci)
        self.remove_home_ref()
        self.encode_variants()
        self.add_loci()
        # Before add_lane_variants and before get_variants, which is the point: both
        # pick a row per token and they pick different ones. After this there is one row
        # to pick.
        self.reconcile_duplicate_rows()
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
            GT = gt_of(sample_field)
            if not is_haploid_gt(GT) or not gt_names_an_allele(GT):
                continue
            try:
                position = int(pos)
            except (TypeError, ValueError):
                continue
            if is_single_copy(chrom, position, self.reference_genome):
                haploid.add(chrom)

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
            # A no-call belongs in neither set. '.' is haploid by shape and './.' is
            # diploid by shape, so recording them would put two spellings of the same
            # 'I do not know' on opposite sides of the copy number question - one voting
            # for a single copy and the other against it.
            if not GT or not gt_names_an_allele(GT):
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

    def reconcile_duplicate_rows(self) -> None:
        """Collapse the rows a file gives more than once for one small variant.

        A VCF can carry several rows for the same variant - a run configured to emit
        both a targeted caller's calls and the general caller's does exactly that, and
        marks in FILTER which row conflicts. Three places then pick a row, and they do
        not pick the same one: get_variants keeps the last it reads, so the genotype,
        the phase and every metric come from there; the FILTER lookups take the first;
        and add_lane_variants reads the first to decide whether to synthesise the '_ref'
        partner, then copies that row to make it. Nobody chose that, and it is why a
        sample can hold a lane locus whose two halves disagree about phase, the '_ref'
        side carrying a phase set from the row whose genotype was thrown away.

        None of those three sites is wrong on its own. What is wrong is that there is
        more than one row to pick from, so this removes the choice rather than
        arbitrating it three times. It runs before all of them.

        What it does, per token:

        - Rows that disagree about the call are left alone. Reconciling them would mean
          deciding which caller to believe, which is not something this layer can know.
          The duplicate genotype warning already names them.
        - Rows that agree collapse to one. The kept row is a phased one wherever any row
          is phased, since an unphased genotype does not contradict a phased one - it is
          the same call, less specific - and today's reading discards the phase in 825
          heterozygous calls for no reason but file order. Where the rows tie, the last
          is kept, which is what get_variants did.
        - The kept row's FILTER becomes every distinct value, ';' joined. Not a new
          convention: a FILTER field already carries several values that way, and
          filter_excludes_allele already splits on it and excludes if any one of them
          says the call is doubtful. Joining is what makes the verdict independent of
          which row the file lists first, and it errs towards excluding.

        Only single base substitutions are touched - see is_single_base_substitution for
        why the bar sits there rather than at "not structural". Rows the structural
        reader wants are never removed, whatever --min_size is set to. A token whose
        rows round together because their sizes do is a different problem: those are two
        real events sharing a name, not two readings of one.
        """
        if "variant" not in self.df.columns or self.df.empty:
            return

        droppable: list[int] = []
        for token, group in self.df.groupby("variant", sort=False):
            if len(group) < 2:
                continue
            rows = list(group.itertuples(index=True, name="Row"))
            if not all(
                is_single_base_substitution(row.REF, row.ALT, row.INFO) for row in rows
            ):
                continue
            gts = [gt_of(row.SAMPLE) for row in rows]
            if not rows_make_the_same_call(gts):
                continue

            keep = rows[-1]
            for row, gt in zip(rows, gts):
                if "|" in gt:
                    keep = row
                    break

            values: list[str] = []
            for row in rows:
                value = str(row.FILTER)
                if value not in values:
                    values.append(value)
            self.df.loc[keep.Index, "FILTER"] = FILTER_VALUE_SEPARATOR.join(values)
            droppable.extend(row.Index for row in rows if row.Index != keep.Index)

        if droppable:
            logger.info(
                f"{self.sample}: {len(droppable)} duplicate row/s were reconciled - a "
                f"variant reported more than once, by rows that agree about the call. "
                f"One row is kept, preferring a phased genotype, and it carries every "
                f"FILTER value the rows had"
            )
            self.df = self.df.drop(index=droppable)

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
                    rows_here = self.df.loc[self.df.loci == lane_loci]
                    alts_here = set(
                        zip(rows_here["REF"].tolist(), rows_here["ALT"].tolist())
                    )
                    if (
                        GT.startswith(("0/1", "0|1", "1/0", "1|0"))
                        and len(alts_here) == 1
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

        A token can be written more than once, because a file can carry more than one
        row for the same variant - two callers both reporting, one of them flagged in
        FILTER as the one that conflicts. The dict keeps the last row read, which is a
        decision nobody made: it is whichever row the file happens to end with. Where
        those rows agree on the genotype the choice does not matter, and where they do
        not the genotype used is a property of the file's order, so that case is warned
        about rather than resolved here. Resolving it means deciding which row to
        believe, and that is not something this function can know.

        Returns:
            dict[str, str]: A dictionary of variants and their associated metrics.
        """
        vcf_variants = {}
        genotypes_seen: dict[str, list[str]] = defaultdict(list)
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
                    genotypes_seen[alt_variant].append(per_alt_metrics["GT"])
                    vcf_variants[alt_variant] = per_alt_metrics
            else:
                genotypes_seen[variant].append(mapped_metrics["GT"])
                vcf_variants[variant] = mapped_metrics

        self._warn_about_contradicted_tokens(genotypes_seen)

        return vcf_variants

    def _warn_about_contradicted_tokens(
        self, genotypes_seen: dict[str, list[str]]
    ) -> None:
        """Say out loud where the file gave one variant more than one genotype.

        get_variants keeps the last row it reads for a token and discards the earlier
        ones. That is silent, and where the discarded row said something different it is
        also load bearing: the pool, the zygosity and the phase all come from whichever
        row the file ended with. Nothing else in the pipeline can see that a second row
        existed, so the only place it can be reported is here.

        Named per sample rather than per token, because the same handful of positions
        recur across a whole cohort and a line each would bury the answer. The genotypes
        are shown in the order the rows were read, so the last one listed is the one in
        use.

        A token every row agreed on is not reported. The overwrite happened, but it
        changed nothing, and a warning that fires on a file where nothing is at stake
        teaches the reader to ignore it.

        Args:
            genotypes_seen (dict[str, list[str]]): Token -> every genotype stored
            against it, in the order the rows were read.
        """
        contradicted = sorted(
            token for token, gts in genotypes_seen.items() if len(set(gts)) > 1
        )
        if not contradicted:
            return
        shown = ", ".join(
            f"{token} ({' then '.join(genotypes_seen[token])})"
            for token in contradicted[:MAX_CONTRADICTED_TOKENS_LOGGED]
        )
        if len(contradicted) > MAX_CONTRADICTED_TOKENS_LOGGED:
            shown += f" and {len(contradicted) - MAX_CONTRADICTED_TOKENS_LOGGED} more"
        logger.warning(
            f"{self.sample}: {len(contradicted)} variant/s were reported on more than "
            f"one row, with a different genotype each time. The last row in the file "
            f"wins, so the genotype used here is a property of the file's order rather "
            f"than of the data: {shown}"
        )


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
    ploidy_scan: PloidyScan | None = None,
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
        ploidy_scan (PloidyScan | None): Collects constitutional haploidy evidence as the
            file goes past, before either filter can discard it. None skips the scan
            entirely, which is what every caller written before it did.

    Returns:
        pl.DataFrame: DataFrame containing the VCF data."""

    open_func = gzip.open if vcf_path.endswith(".gz") else open
    header = None
    samples: list[str] = []
    rows: list[str] = []
    with open_func(vcf_path, "rt") as f:
        for line in f:  # TODO Pool
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                header = line.lstrip("#").strip().split("\t")
                if len(header) == 10:
                    header[-1] = "SAMPLE"  # for single sample
                samples = header[9:]
                continue

            # parse variant
            fields = line.split("\t")
            try:
                chrom, pos = fields[0].removeprefix("chr"), int(fields[1])
            except:
                raise
            # Before either filter, because both of them throw this evidence away - see
            # PloidyScan. Costs one string comparison on every non-X/Y row.
            if ploidy_scan is not None:
                ploidy_scan.observe(chrom, pos, samples, fields)
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
