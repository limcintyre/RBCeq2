"""Variant encoders for VCF data from various structural variant callers.

This module provides a flexible encoding system for converting VCF variant records
into human-readable string representations. It supports multiple SV callers:

Long-read callers:
    - Sniffles2
    - SVIM
    - CuteSV
    - NanoVar
    - pbsv (PacBio)

Short-read callers:
    - Manta
    - LUMPY
    - DELLY
    - GRIDSS
    - CNVnator
    - GATK-gCNV
    - DRAGEN
    - SurVIndel2

As well as small variants (SNPs/indels).

Example:
    >>> import pandas as pd
    >>> factory = VariantEncoderFactory()
    >>> row = pd.Series({
    ...     "CHROM": "chr1",
    ...     "POS": 1000,
    ...     "REF": "A",
    ...     "ALT": "<DEL>",
    ...     "INFO": "SVTYPE=DEL;SVLEN=-500;END=1500"
    ... })
    >>> factory.encode_variant(row)
    'chr1:1000_del_500bp'
"""

from __future__ import annotations

import logging
import re
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Protocol

if TYPE_CHECKING:
    import pandas as pd
else:
    try:
        import pandas as pd
    except ImportError:
        pd = None  # type: ignore[assignment]

# Set up logging - use loguru if available, otherwise standard logging
from loguru import logger



class VariantEncoder(Protocol):
    """Protocol defining the interface for variant encoders.

    All variant encoders must implement this protocol to be used with
    the VariantEncoderFactory. Encoders are responsible for detecting
    whether they can handle a specific variant type and encoding it
    into a standardized string format.
    """

    def can_encode(self, row: pd.Series) -> bool:
        """Check if this encoder can handle the given variant row.

        Args:
            row: A pandas Series representing one row from a VCF DataFrame.
                Expected to contain at minimum: CHROM, POS, REF, ALT columns.
                May also contain INFO column for structural variants.

        Returns:
            True if this encoder can handle the variant, False otherwise.
        """
        ...

    def encode(self, row: pd.Series) -> str:
        """Encode the variant into a human-readable string format.

        Args:
            row: A pandas Series representing one row from a VCF DataFrame.

        Returns:
            Encoded variant string in the format appropriate for the variant type.
            Typically: "chrom:pos_type_size" for SVs or "chrom:pos_ref_alt" for
            small variants.
        """
        ...


def parse_info_field(info: str) -> dict[str, str | bool]:
    """Parse a VCF INFO field into a dictionary.

    Handles both key=value pairs and flag fields (keys without values).

    Args:
        info: The INFO field string from a VCF record, with fields
            separated by semicolons (e.g., "SVTYPE=DEL;SVLEN=-500;IMPRECISE").

    Returns:
        Dictionary mapping INFO keys to their values. Flag fields (without
        values) are mapped to True.

    Example:
        >>> parse_info_field("SVTYPE=DEL;SVLEN=-500;IMPRECISE")
        {'SVTYPE': 'DEL', 'SVLEN': '-500', 'IMPRECISE': True}
    """
    result: dict[str, str | bool] = {}
    if not info or info == "." or pd.isna(info):
        return result

    for item in str(info).split(";"):
        if not item:
            continue
        if "=" in item:
            key, value = item.split("=", 1)
            result[key] = value
        else:
            result[item] = True

    return result


def format_size(bp: int) -> str:
    """Format a base pair size as a human-readable string.

    Converts base pair counts to kilobase or megabase notation when appropriate.

    Args:
        bp: Size in base pairs (should be a positive integer).

    Returns:
        Formatted size string (e.g., "500bp", "78kb", or "2mb").

    Example:
        >>> format_size(500)
        '500bp'
        >>> format_size(78874)
        '78kb'
        >>> format_size(2500000)
        '2mb'
    """
    if bp >= 1_000_000:
        return f"{bp // 1_000_000}mb"
    if bp >= 1000:
        return f"{bp // 1000}kb"
    return f"{bp}bp"


def get_svlen_value(svlen_str: str) -> int | None:
    """Extract the first SVLEN value, handling multi-allelic cases.

    SVLEN can be comma-separated for multi-allelic SVs (e.g., "-500,-1000").
    This function extracts the first value and returns its absolute value.

    Args:
        svlen_str: The SVLEN value string from the INFO field.

    Returns:
        Absolute value of the first SVLEN entry, or None if parsing fails.

    Example:
        >>> get_svlen_value("-500")
        500
        >>> get_svlen_value("-500,-1000")
        500
    """
    try:
        first_value = svlen_str.split(",")[0]
        return abs(int(first_value))
    except (ValueError, IndexError):
        return None


def validate_required_columns(row: pd.Series) -> None:
    """Validate that required VCF columns are present in the row.

    Args:
        row: A pandas Series representing one row from a VCF DataFrame.

    Raises:
        KeyError: If any required column (CHROM, POS, REF, ALT) is missing.
    """
    required = ["CHROM", "POS", "REF", "ALT"]
    missing = [col for col in required if col not in row.index]
    if missing:
        raise KeyError(f"Missing required VCF columns: {missing}")


def is_breakend_alt(alt: str) -> bool:
    """Check if an ALT field represents a breakend (BND) variant.

    Breakend notation uses brackets to indicate adjacency, e.g.:
    - ]13:123456]T (join after position on reverse strand)
    - T[13:123456[ (join before position on forward strand)

    Args:
        alt: The ALT field value from a VCF record.

    Returns:
        True if the ALT field contains breakend notation.
    """
    return bool(alt) and ("[" in alt or "]" in alt) and not alt.startswith("<")


def is_symbolic_alt(alt: str) -> bool:
    """Check if an ALT field is a symbolic allele (e.g., <DEL>, <INS>).

    Args:
        alt: The ALT field value from a VCF record.

    Returns:
        True if the ALT field is a symbolic structural variant allele.
    """
    return bool(alt) and alt.startswith("<") and alt.endswith(">")


def normalize_svtype(svtype: str) -> str:
    """Normalize an SVTYPE string to a consistent lowercase format.

    Handles subtypes like DUP:TANDEM -> dup_tandem.

    Args:
        svtype: The SVTYPE value from the INFO field.

    Returns:
        Normalized lowercase svtype with colons replaced by underscores.
    """
    return svtype.lower().replace(":", "_")


def compute_svlen(info: dict[str, str | bool], pos: int) -> int | None:
    """Compute SV length from INFO fields, trying SVLEN first then END.

    Args:
        info: Parsed INFO field dictionary.
        pos: The POS value from the VCF record.

    Returns:
        Absolute SV length in base pairs, or None if not determinable.
    """
    # Try SVLEN first
    if "SVLEN" in info:
        svlen = get_svlen_value(str(info["SVLEN"]))
        if svlen is not None and svlen > 0:
            return svlen

    # Fall back to END
    if "END" in info:
        try:
            end = int(info["END"])
            return abs(end - pos)
        except (ValueError, TypeError):
            pass

    return None


def encode_sv_standard(
    chrom: str, pos: int, svtype: str, svlen: int | None
) -> str:
    """Encode an SV in standard format: chrom:pos_svtype_size.

    Args:
        chrom: Chromosome name.
        pos: Position.
        svtype: SV type (will be normalized).
        svlen: SV length in base pairs (optional).

    Returns:
        Encoded SV string.
    """
    normalized_type = normalize_svtype(svtype)
    if svlen is not None and svlen > 0:
        size_str = format_size(svlen)
        return f"{chrom}:{pos}_{normalized_type}_{size_str}"
    return f"{chrom}:{pos}_{normalized_type}"


@dataclass(slots=True, frozen=True)
class SmallVariantEncoder:
    """Encoder for small variants (SNPs and small indels).

    Handles variants that are not structural variants, encoding them
    in the format: chrom:pos_ref_alt

    For multi-allelic sites, produces comma-separated encodings for
    each alternate allele.

    Example:
        Input: CHROM=chr1, POS=1000, REF=A, ALT=G
        Output: "chr1:1000_A_G"

        Input: CHROM=chr1, POS=1000, REF=A, ALT=G,T
        Output: "chr1:1000_A_G,chr1:1000_A_T"
    """

    def can_encode(self, row: pd.Series) -> bool:
        """Check if variant is a small variant (not a structural variant).

        A variant is considered "small" if it has explicit REF/ALT sequences
        (not symbolic alleles) and no SVTYPE in the INFO field.

        Args:
            row: A pandas Series representing one row from a VCF DataFrame.

        Returns:
            True if this is a small variant that can be encoded.
        """
        alt = str(row.get("ALT", ""))
        info_str = str(row.get("INFO", ""))
        info = parse_info_field(info_str)

        return (
            not is_symbolic_alt(alt)
            and not is_breakend_alt(alt)
            and "SVTYPE" not in info
            and "END" not in info
        )

    def encode(self, row: pd.Series) -> str:
        """Encode small variant as chrom:pos_ref_alt.

        Args:
            row: A pandas Series representing one row from a VCF DataFrame.

        Returns:
            Encoded variant string(s). For multi-allelic sites, returns
            comma-separated encodings.
        """
        chrom = str(row["CHROM"])
        pos = str(row["POS"])
        ref = str(row["REF"])
        alts = str(row["ALT"]).split(",")

        return ",".join([f"{chrom}:{pos}_{ref}_{alt}" for alt in alts])


@dataclass(slots=True, frozen=True)
class BreakendEncoder:
    """Encoder for breakend (BND/translocation) variants.

    Handles breakend notation where the ALT field contains bracket notation
    indicating joins between different genomic positions. This encoder is
    used by multiple callers including LUMPY, GRIDSS, Manta, SVIM, pbsv,
    NanoVar, and others.

    Example:
        Input: CHROM=chr1, POS=1000, ALT=]chr2:5000]T, INFO=SVTYPE=BND
        Output: "chr1:1000_bnd_chr2:5000"
    """

    def can_encode(self, row: pd.Series) -> bool:
        """Check if variant is a breakend variant.

        Args:
            row: A pandas Series representing one row from a VCF DataFrame.

        Returns:
            True if this is a breakend variant.
        """
        alt = str(row.get("ALT", ""))
        return is_breakend_alt(alt)

    def encode(self, row: pd.Series) -> str:
        """Encode breakend variant with partner location.

        Args:
            row: A pandas Series representing one row from a VCF DataFrame.

        Returns:
            Encoded breakend string showing both breakpoint positions.
        """
        chrom = str(row["CHROM"])
        pos = str(row["POS"])
        alt = str(row["ALT"])

        # Extract partner position from breakend notation
        # Formats: ]chr:pos]seq, [chr:pos[seq, seq]chr:pos], seq[chr:pos[
        match = re.search(r"[\[\]]([^:\[\]]+):(\d+)[\[\]]", alt)
        if match:
            partner_chrom = match.group(1)
            partner_pos = match.group(2)
            return f"{chrom}:{pos}_bnd_{partner_chrom}:{partner_pos}"

        return f"{chrom}:{pos}_bnd"


# =============================================================================
# Long-read SV Callers
# =============================================================================



@dataclass(slots=True, frozen=True)
class Sniffles2Encoder:
    """Encoder for Sniffles2 structural variants.

    Sniffles2 outputs SVs with SVTYPE and SVLEN in the INFO field.
    It may also include RE (read support) field.

    Example:
        Input: INFO=SVTYPE=DEL;SVLEN=-78874;RE=15
        Output: "4:143990187_del_78kb"
    """

    def can_encode(self, row: pd.Series) -> bool:
        """Check if variant is from Sniffles2 (has SVTYPE and SVLEN, may have RE).

        Args:
            row: A pandas Series representing one row from a VCF DataFrame.

        Returns:
            True if this appears to be a Sniffles2 structural variant.
        """
        alt = str(row.get("ALT", ""))
        info_str = str(row.get("INFO", ""))
        info = parse_info_field(info_str)

        # Sniffles2 has SVTYPE and SVLEN, may have RE (read evidence)
        # Exclude if has SUPPORT (SVIM) or SU/PE/SR (LUMPY)
        return (
            "SVTYPE" in info
            and "SVLEN" in info
            and "SUPPORT" not in info  # Not SVIM
            and "SU" not in info  # Not LUMPY
            and "PE" not in info
            and "SR" not in info
            and "SVANN" not in info  # Not pbsv
            and not is_breakend_alt(alt)
        )

    def encode(self, row: pd.Series) -> str:
        """Encode as chrom:pos_svtype_size.

        Args:
            row: A pandas Series representing one row from a VCF DataFrame.

        Returns:
            Encoded SV string (e.g., "4:143990187_del_78kb").
        """
        chrom = str(row["CHROM"])
        pos = int(row["POS"])
        info = parse_info_field(str(row.get("INFO", "")))

        svtype = str(info.get("SVTYPE", "unknown"))
        svlen = compute_svlen(info, pos)

        return encode_sv_standard(chrom, pos, svtype, svlen)


@dataclass(slots=True, frozen=True)
class SvimEncoder:
    """Encoder for SVIM structural variants.

    SVIM (Structural Variant Identification using Mapped long reads) outputs
    SVs with SVTYPE, END, SVLEN, and SUPPORT in the INFO field. SVIM uses
    symbolic ALT alleles like <DEL>, <INS>, <DUP>, <DUP:TANDEM>, <DUP:INT>,
    <INV>, and <BND>.

    SVIM is distinguished by having SUPPORT field (read support count).

    Example:
        Input: ALT=<DEL>, INFO=SVTYPE=DEL;END=1500;SVLEN=-500;SUPPORT=10
        Output: "chr1:1000_del_500bp"
    """

    def can_encode(self, row: pd.Series) -> bool:
        """Check if variant is from SVIM (has SVTYPE and SUPPORT).

        Args:
            row: A pandas Series representing one row from a VCF DataFrame.

        Returns:
            True if this appears to be a SVIM structural variant.
        """
        alt = str(row.get("ALT", ""))
        info_str = str(row.get("INFO", ""))
        info = parse_info_field(info_str)

        return (
            "SVTYPE" in info and "SUPPORT" in info and not is_breakend_alt(alt)
        )

    def encode(self, row: pd.Series) -> str:
        """Encode SVIM variant as chrom:pos_svtype_size.

        Args:
            row: A pandas Series representing one row from a VCF DataFrame.

        Returns:
            Encoded SV string (e.g., "chr1:1000_del_500bp").
        """
        chrom = str(row["CHROM"])
        pos = int(row["POS"])
        info = parse_info_field(str(row.get("INFO", "")))

        svtype = str(info.get("SVTYPE", "unknown"))
        svlen = compute_svlen(info, pos)

        return encode_sv_standard(chrom, pos, svtype, svlen)


@dataclass(slots=True, frozen=True)
class CuteSVEncoder:
    """Encoder for CuteSV structural variants.

    CuteSV outputs SVs with SVTYPE, SVLEN, END, and RE (read evidence) fields.
    It uses symbolic ALT alleles and supports DEL, INS, DUP, INV, TRA, and BND.

    Example:
        Input: INFO=SVTYPE=DEL;SVLEN=-500;END=1500;RE=8
        Output: "chr1:1000_del_500bp"
    """

    def can_encode(self, row: pd.Series) -> bool:
        """Check if variant is from CuteSV (has SVTYPE and RE).

        Args:
            row: A pandas Series representing one row from a VCF DataFrame.

        Returns:
            True if this appears to be a CuteSV structural variant.
        """
        alt = str(row.get("ALT", ""))
        info_str = str(row.get("INFO", ""))
        info = parse_info_field(info_str)

        # CuteSV uses RE for read evidence
        return (
            "SVTYPE" in info
            and "RE" in info
            and "SUPPORT" not in info  # Not SVIM
            and not is_breakend_alt(alt)
        )

    def encode(self, row: pd.Series) -> str:
        """Encode CuteSV variant as chrom:pos_svtype_size.

        Args:
            row: A pandas Series representing one row from a VCF DataFrame.

        Returns:
            Encoded SV string.
        """
        chrom = str(row["CHROM"])
        pos = int(row["POS"])
        info = parse_info_field(str(row.get("INFO", "")))

        svtype = str(info.get("SVTYPE", "unknown"))
        svlen = compute_svlen(info, pos)

        return encode_sv_standard(chrom, pos, svtype, svlen)


@dataclass(slots=True, frozen=True)
class NanoVarEncoder:
    """Encoder for NanoVar structural variants.

    NanoVar outputs SVs with SVTYPE, END, SVLEN, SR (supporting reads),
    and NN (neural network score). It also includes SV2 field for
    additional SV classification.

    Example:
        Input: INFO=SVTYPE=DEL;END=1500;SVLEN=-500;SR=5;NN=0.95
        Output: "chr1:1000_del_500bp"
    """

    def can_encode(self, row: pd.Series) -> bool:
        """Check if variant is from NanoVar (has SVTYPE and NN score).

        Args:
            row: A pandas Series representing one row from a VCF DataFrame.

        Returns:
            True if this appears to be a NanoVar structural variant.
        """
        alt = str(row.get("ALT", ""))
        info_str = str(row.get("INFO", ""))
        info = parse_info_field(info_str)

        # NanoVar uses NN for neural network confidence score
        return (
            "SVTYPE" in info and "NN" in info and not is_breakend_alt(alt)
        )

    def encode(self, row: pd.Series) -> str:
        """Encode NanoVar variant as chrom:pos_svtype_size.

        Args:
            row: A pandas Series representing one row from a VCF DataFrame.

        Returns:
            Encoded SV string.
        """
        chrom = str(row["CHROM"])
        pos = int(row["POS"])
        info = parse_info_field(str(row.get("INFO", "")))

        svtype = str(info.get("SVTYPE", "unknown"))
        svlen = compute_svlen(info, pos)

        return encode_sv_standard(chrom, pos, svtype, svlen)


@dataclass(slots=True, frozen=True)
class PbsvEncoder:
    """Encoder for pbsv (PacBio) structural variants.

    pbsv outputs SVs with SVTYPE, SVLEN, END, and often SVANN (annotation
    like TANDEM). It uses symbolic ALT alleles for DEL, INS, DUP, INV, BND.

    Example:
        Input: INFO=SVTYPE=DEL;END=904587;SVLEN=-97;SVANN=TANDEM
        Output: "chr1:904490_del_97bp"
    """

    def can_encode(self, row: pd.Series) -> bool:
        """Check if variant is from pbsv (has SVTYPE and often SVANN).

        Args:
            row: A pandas Series representing one row from a VCF DataFrame.

        Returns:
            True if this appears to be a pbsv structural variant.
        """
        alt = str(row.get("ALT", ""))
        info_str = str(row.get("INFO", ""))
        info = parse_info_field(info_str)

        # pbsv often has SVANN annotation
        return (
            "SVTYPE" in info
            and ("SVANN" in info or "MATEID" in info)
            and not is_breakend_alt(alt)
        )

    def encode(self, row: pd.Series) -> str:
        """Encode pbsv variant as chrom:pos_svtype_size.

        Args:
            row: A pandas Series representing one row from a VCF DataFrame.

        Returns:
            Encoded SV string.
        """
        chrom = str(row["CHROM"])
        pos = int(row["POS"])
        info = parse_info_field(str(row.get("INFO", "")))

        svtype = str(info.get("SVTYPE", "unknown"))
        svlen = compute_svlen(info, pos)

        return encode_sv_standard(chrom, pos, svtype, svlen)


# =============================================================================
# Short-read SV Callers
# =============================================================================


@dataclass(slots=True, frozen=True)
class LumpyEncoder:
    """Encoder for LUMPY structural variants.

    LUMPY outputs SVs with SVTYPE, SVLEN, END, and evidence counts:
    SU (total support), PE (paired-end), SR (split-read).
    Also includes STRANDS field.

    Example:
        Input: INFO=SVTYPE=DEL;SVLEN=-2887;END=93825700;SU=34;PE=22;SR=12
        Output: "chr1:93822813_del_2kb"
    """

    def can_encode(self, row: pd.Series) -> bool:
        """Check if variant is from LUMPY (has SVTYPE and SU/PE/SR).

        Args:
            row: A pandas Series representing one row from a VCF DataFrame.

        Returns:
            True if this appears to be a LUMPY structural variant.
        """
        alt = str(row.get("ALT", ""))
        info_str = str(row.get("INFO", ""))
        info = parse_info_field(info_str)

        # LUMPY uses SU (support), PE (paired-end), SR (split-read)
        return (
            "SVTYPE" in info
            and "SU" in info
            and not is_breakend_alt(alt)
        )

    def encode(self, row: pd.Series) -> str:
        """Encode LUMPY variant as chrom:pos_svtype_size.

        Args:
            row: A pandas Series representing one row from a VCF DataFrame.

        Returns:
            Encoded SV string.
        """
        chrom = str(row["CHROM"])
        pos = int(row["POS"])
        info = parse_info_field(str(row.get("INFO", "")))

        svtype = str(info.get("SVTYPE", "unknown"))
        svlen = compute_svlen(info, pos)

        return encode_sv_standard(chrom, pos, svtype, svlen)


@dataclass(slots=True, frozen=True)
class DellyEncoder:
    """Encoder for DELLY structural variants.

    DELLY outputs SVs with SVTYPE, END, PE (paired-end reads), MAPQ,
    and CT (connection type like 3to5). Does not include SVLEN.

    Example:
        Input: INFO=SVTYPE=DEL;END=144069060;PE=180;MAPQ=44;CT=3to5
        Output: "4:143990187_del_78kb"
    """

    def can_encode(self, row: pd.Series) -> bool:
        """Check if variant is from DELLY (has SVTYPE, END, CT but not SVLEN).

        Args:
            row: A pandas Series representing one row from a VCF DataFrame.

        Returns:
            True if this appears to be a DELLY structural variant.
        """
        alt = str(row.get("ALT", ""))
        info_str = str(row.get("INFO", ""))
        info = parse_info_field(info_str)

        # DELLY has CT (connection type) and often MAPQ
        return (
            "SVTYPE" in info
            and "END" in info
            and ("CT" in info or "SVMETHOD" in info)
            and "SVLEN" not in info
            and "SUPPORT" not in info  # Not SVIM
            and "SU" not in info  # Not LUMPY
            and not is_breakend_alt(alt)
        )

    def encode(self, row: pd.Series) -> str:
        """Encode DELLY variant as chrom:pos_svtype_size.

        Args:
            row: A pandas Series representing one row from a VCF DataFrame.

        Returns:
            Encoded SV string.
        """
        chrom = str(row["CHROM"])
        pos = int(row["POS"])
        info = parse_info_field(str(row.get("INFO", "")))

        svtype = str(info.get("SVTYPE", "unknown"))

        # DELLY doesn't have SVLEN, compute from END
        svlen = None
        if "END" in info:
            try:
                end = int(info["END"])
                svlen = abs(end - pos)
            except (ValueError, TypeError):
                pass

        return encode_sv_standard(chrom, pos, svtype, svlen)


@dataclass(slots=True, frozen=True)
class MantaEncoder:
    """Encoder for Manta structural variants.

    Manta outputs SVs with symbolic ALT alleles (<DEL>, <INS>, etc.)
    and SVTYPE in INFO. It may have either SVLEN or END (or both).
    Often includes CIPOS, CIEND for confidence intervals.

    This encoder is used as a general fallback for symbolic SV alleles
    that don't match more specific caller signatures.

    Example:
        Input: ALT=<DEL>, INFO=SVTYPE=DEL;SVLEN=-5000;END=6000
        Output: "chr1:1000_del_5kb"
    """

    def can_encode(self, row: pd.Series) -> bool:
        """Check if variant is a symbolic SV (generic handler for Manta-style).

        Args:
            row: A pandas Series representing one row from a VCF DataFrame.

        Returns:
            True if this is a symbolic structural variant.
        """
        alt = str(row.get("ALT", ""))
        info_str = str(row.get("INFO", ""))
        info = parse_info_field(info_str)

        return (
            is_symbolic_alt(alt)
            and "SVTYPE" in info
            and not is_breakend_alt(alt)
        )

    def encode(self, row: pd.Series) -> str:
        """Encode Manta-style variant as chrom:pos_svtype_size.

        Args:
            row: A pandas Series representing one row from a VCF DataFrame.

        Returns:
            Encoded SV string.
        """
        chrom = str(row["CHROM"])
        pos = int(row["POS"])
        info = parse_info_field(str(row.get("INFO", "")))

        svtype = str(info.get("SVTYPE", "unknown"))
        svlen = compute_svlen(info, pos)

        return encode_sv_standard(chrom, pos, svtype, svlen)


@dataclass(slots=True, frozen=True)
class GridssEncoder:
    """Encoder for GRIDSS structural variants.

    GRIDSS outputs all variants as BND (breakend) notation, even for simple
    events. It includes AS, RAS, CAS assembly fields and quality scores.
    The SIMPLE_TYPE annotation (if added by post-processing) indicates
    the actual SV type.

    Since GRIDSS uses BND for everything, breakends are handled by
    BreakendEncoder. This encoder handles annotated GRIDSS output.

    Example:
        Input: INFO=SVTYPE=BND;AS=1;EVENT=gridss1
        Output: Handled by BreakendEncoder
    """

    def can_encode(self, row: pd.Series) -> bool:
        """Check if variant is from GRIDSS (has AS/RAS/CAS assembly fields).

        Args:
            row: A pandas Series representing one row from a VCF DataFrame.

        Returns:
            True if this appears to be a GRIDSS structural variant.
        """
        info_str = str(row.get("INFO", ""))
        info = parse_info_field(info_str)

        # GRIDSS has assembly-specific fields
        return (
            "SVTYPE" in info
            and ("AS" in info or "RAS" in info or "EVENT" in info)
            and "SIMPLE_TYPE" in info  # Post-processed with SV type annotation
        )

    def encode(self, row: pd.Series) -> str:
        """Encode GRIDSS variant using SIMPLE_TYPE annotation.

        Args:
            row: A pandas Series representing one row from a VCF DataFrame.

        Returns:
            Encoded SV string.
        """
        chrom = str(row["CHROM"])
        pos = int(row["POS"])
        info = parse_info_field(str(row.get("INFO", "")))

        # Use SIMPLE_TYPE if available (from annotation), otherwise SVTYPE
        svtype = str(info.get("SIMPLE_TYPE", info.get("SVTYPE", "bnd")))
        svlen = compute_svlen(info, pos)

        return encode_sv_standard(chrom, pos, svtype, svlen)


@dataclass(slots=True, frozen=True)
class CnvnatorEncoder:
    """Encoder for CNVnator copy number variants.

    CNVnator outputs CNVs with SVTYPE (DEL or DUP), END, and often
    includes normalized read depth (natorRD) and other CNV-specific fields.

    Example:
        Input: INFO=SVTYPE=DEL;END=2000;natorRD=0.5
        Output: "chr1:1000_del_1kb"
    """

    def can_encode(self, row: pd.Series) -> bool:
        """Check if variant is from CNVnator (has natorRD or specific signature).

        Args:
            row: A pandas Series representing one row from a VCF DataFrame.

        Returns:
            True if this appears to be a CNVnator variant.
        """
        alt = str(row.get("ALT", ""))
        info_str = str(row.get("INFO", ""))
        info = parse_info_field(info_str)

        # CNVnator may have natorRD field or be identified by source
        return (
            "SVTYPE" in info
            and ("natorRD" in info or "natorP1" in info or "natorP2" in info)
            and not is_breakend_alt(alt)
        )

    def encode(self, row: pd.Series) -> str:
        """Encode CNVnator variant as chrom:pos_svtype_size.

        Args:
            row: A pandas Series representing one row from a VCF DataFrame.

        Returns:
            Encoded SV string.
        """
        chrom = str(row["CHROM"])
        pos = int(row["POS"])
        info = parse_info_field(str(row.get("INFO", "")))

        svtype = str(info.get("SVTYPE", "cnv"))
        svlen = compute_svlen(info, pos)

        return encode_sv_standard(chrom, pos, svtype, svlen)


@dataclass(slots=True, frozen=True)
class GatkGcnvEncoder:
    """Encoder for GATK-gCNV copy number variants.

    GATK-gCNV outputs CNVs with SVTYPE, END, and CN (copy number).
    Uses ALGORITHMS field to identify the caller.

    Example:
        Input: INFO=SVTYPE=DEL;END=2000;CN=1;ALGORITHMS=depth
        Output: "chr1:1000_del_1kb"
    """

    def can_encode(self, row: pd.Series) -> bool:
        """Check if variant is from GATK-gCNV (has ALGORITHMS field).

        Args:
            row: A pandas Series representing one row from a VCF DataFrame.

        Returns:
            True if this appears to be a GATK-gCNV variant.
        """
        alt = str(row.get("ALT", ""))
        info_str = str(row.get("INFO", ""))
        info = parse_info_field(info_str)

        # GATK-SV pipeline uses ALGORITHMS field
        return (
            "SVTYPE" in info
            and "ALGORITHMS" in info
            and not is_breakend_alt(alt)
        )

    def encode(self, row: pd.Series) -> str:
        """Encode GATK-gCNV variant as chrom:pos_svtype_size.

        Args:
            row: A pandas Series representing one row from a VCF DataFrame.

        Returns:
            Encoded SV string.
        """
        chrom = str(row["CHROM"])
        pos = int(row["POS"])
        info = parse_info_field(str(row.get("INFO", "")))

        svtype = str(info.get("SVTYPE", "cnv"))
        svlen = compute_svlen(info, pos)

        return encode_sv_standard(chrom, pos, svtype, svlen)


@dataclass(slots=True, frozen=True)
class DragenEncoder:
    """Encoder for DRAGEN CNV/SV variants.

    DRAGEN outputs CNVs with SVTYPE, END, SVLEN, and REFLEN (reference length).
    Often uses symbolic alleles like <DEL>, <DUP>, <CNV>.

    Example:
        Input: INFO=SVTYPE=DEL;END=2000;SVLEN=-1000;REFLEN=1000
        Output: "chr1:1000_del_1kb"
    """

    def can_encode(self, row: pd.Series) -> bool:
        """Check if variant is from DRAGEN (has REFLEN field).

        Args:
            row: A pandas Series representing one row from a VCF DataFrame.

        Returns:
            True if this appears to be a DRAGEN variant.
        """
        alt = str(row.get("ALT", ""))
        info_str = str(row.get("INFO", ""))
        info = parse_info_field(info_str)

        # DRAGEN uses REFLEN field
        return (
            "SVTYPE" in info
            and "REFLEN" in info
            and not is_breakend_alt(alt)
        )

    def encode(self, row: pd.Series) -> str:
        """Encode DRAGEN variant as chrom:pos_svtype_size.

        Args:
            row: A pandas Series representing one row from a VCF DataFrame.

        Returns:
            Encoded SV string.
        """
        chrom = str(row["CHROM"])
        pos = int(row["POS"])
        info = parse_info_field(str(row.get("INFO", "")))

        svtype = str(info.get("SVTYPE", "cnv"))
        # DRAGEN has REFLEN which is the reference length
        svlen = None
        if "REFLEN" in info:
            try:
                svlen = abs(int(info["REFLEN"]))
            except (ValueError, TypeError):
                pass
        if svlen is None:
            svlen = compute_svlen(info, pos)

        return encode_sv_standard(chrom, pos, svtype, svlen)


@dataclass(slots=True, frozen=True)
class SurVIndel2Encoder:
    """Encoder for SurVIndel2 structural variants.

    SurVIndel2 outputs SVs with standard VCF fields and may include
    specific annotations. Uses SVTYPE, SVLEN, END.

    Example:
        Input: INFO=SVTYPE=DEL;SVLEN=-500;END=1500
        Output: "chr1:1000_del_500bp"
    """

    def can_encode(self, row: pd.Series) -> bool:
        """Check if variant is from SurVIndel2.

        Note: SurVIndel2 uses standard VCF format, so this encoder
        checks for the absence of other caller-specific signatures.

        Args:
            row: A pandas Series representing one row from a VCF DataFrame.

        Returns:
            True if this might be a SurVIndel2 variant.
        """
        # SurVIndel2 uses standard format - this is a placeholder
        # In practice, you'd check for source header or specific fields
        return False  # Let MantaEncoder handle as fallback

    def encode(self, row: pd.Series) -> str:
        """Encode SurVIndel2 variant as chrom:pos_svtype_size.

        Args:
            row: A pandas Series representing one row from a VCF DataFrame.

        Returns:
            Encoded SV string.
        """
        chrom = str(row["CHROM"])
        pos = int(row["POS"])
        info = parse_info_field(str(row.get("INFO", "")))

        svtype = str(info.get("SVTYPE", "unknown"))
        svlen = compute_svlen(info, pos)

        return encode_sv_standard(chrom, pos, svtype, svlen)


@dataclass(slots=True, frozen=True)
class ChallengerEncoder:
    """Encoder for CHALLENGER structural variants.

    CHALLENGER is a short-read SV caller that outputs standard VCF format.

    Example:
        Input: INFO=SVTYPE=DEL;SVLEN=-500;END=1500
        Output: "chr1:1000_del_500bp"
    """

    def can_encode(self, row: pd.Series) -> bool:
        """Check if variant is from CHALLENGER.

        Note: CHALLENGER uses standard VCF format. Detection would
        require checking the VCF header for source.

        Args:
            row: A pandas Series representing one row from a VCF DataFrame.

        Returns:
            False - let MantaEncoder handle as fallback.
        """
        return False  # Let MantaEncoder handle as fallback

    def encode(self, row: pd.Series) -> str:
        """Encode CHALLENGER variant as chrom:pos_svtype_size.

        Args:
            row: A pandas Series representing one row from a VCF DataFrame.

        Returns:
            Encoded SV string.
        """
        chrom = str(row["CHROM"])
        pos = int(row["POS"])
        info = parse_info_field(str(row.get("INFO", "")))

        svtype = str(info.get("SVTYPE", "unknown"))
        svlen = compute_svlen(info, pos)

        return encode_sv_standard(chrom, pos, svtype, svlen)


# =============================================================================
# Factory
# =============================================================================


@dataclass(slots=True, frozen=False)
class VariantEncoderFactory:
    """Factory to manage and apply variant encoders.

    Encoders are tried in order until one matches. The order matters:
    more specific encoders (caller-specific) should come before
    generic ones (like MantaEncoder).

    The default encoder order is:
    1. BreakendEncoder - for BND/translocation variants (all callers)
    2. SvimEncoder - SVIM-specific (has SUPPORT field)
    3. CuteSVEncoder - CuteSV-specific (has RE field)
    4. NanoVarEncoder - NanoVar-specific (has NN field)
    5. PbsvEncoder - pbsv-specific (has SVANN field)
    6. LumpyEncoder - LUMPY-specific (has SU/PE/SR fields)
    7. Sniffles2Encoder - Sniffles2-specific (has SVLEN, no other signatures)
    8. DellyEncoder - DELLY-specific (has CT field, no SVLEN)
    9. GridssEncoder - GRIDSS-specific (has AS/RAS fields)
    10. CnvnatorEncoder - CNVnator-specific (has natorRD)
    11. GatkGcnvEncoder - GATK-specific (has ALGORITHMS)
    12. DragenEncoder - DRAGEN-specific (has REFLEN)
    13. MantaEncoder - generic symbolic SV handler
    14. SmallVariantEncoder - fallback for SNPs/indels

    Attributes:
        encoders: List of encoder instances in priority order.

    Example:
        >>> factory = VariantEncoderFactory()
        >>> encoded = factory.encode_variant(vcf_row)

        # Add a custom encoder with highest priority
        >>> factory.add_encoder(MyCustomEncoder(), priority=0)
    """

    encoders: list[VariantEncoder] = field(
        default_factory=lambda: [
            BreakendEncoder(),
            # Long-read callers
            SvimEncoder(),
            CuteSVEncoder(),
            NanoVarEncoder(),
            PbsvEncoder(),
            # Short-read callers
            LumpyEncoder(),
            Sniffles2Encoder(),
            DellyEncoder(),
            GridssEncoder(),
            CnvnatorEncoder(),
            GatkGcnvEncoder(),
            DragenEncoder(),
            # Generic fallbacks
            MantaEncoder(),
            SmallVariantEncoder(),
        ]
    )

    def add_encoder(self, encoder: VariantEncoder, priority: int = 0) -> None:
        """Add a custom encoder at specified priority position.

        Args:
            encoder: The encoder instance to add. Must implement the
                VariantEncoder protocol.
            priority: Position in the encoder list where 0 is highest
                priority (tried first). Use len(encoders) to add at end.

        Example:
            >>> factory = VariantEncoderFactory()
            >>> factory.add_encoder(MyEncoder(), priority=0)  # Highest priority
            >>> factory.add_encoder(MyEncoder(), priority=len(factory.encoders))
        """
        self.encoders.insert(priority, encoder)

    def encode_variant(self, row: pd.Series) -> str:
        """Encode a variant using the first matching encoder.

        Iterates through encoders in priority order, using the first
        encoder that can handle the variant.

        Args:
            row: A pandas Series representing one row from the VCF DataFrame.
                Must contain CHROM, POS, REF, ALT columns. May contain INFO.

        Returns:
            Encoded variant string in the format appropriate for the variant
            type.

        Raises:
            KeyError: If required columns (CHROM, POS, REF, ALT) are missing.

        Example:
            >>> row = pd.Series({"CHROM": "chr1", "POS": 1000, "REF": "A",
            ...                  "ALT": "G", "INFO": "."})
            >>> factory.encode_variant(row)
            'chr1:1000_A_G'
        """
        validate_required_columns(row)

        for encoder in self.encoders:
            if encoder.can_encode(row):
                return encoder.encode(row)

        # Ultimate fallback if no encoder matches
        logger.warning(
            f"No encoder matched variant at {row['CHROM']}:{row['POS']}. "
            "Using raw format."
        )
        return f"{row['CHROM']}:{row['POS']}_{row['REF']}_{row['ALT']}"
    
# """Variant encoders for VCF data from various structural variant callers.

# This module provides a flexible encoding system for converting VCF variant records
# into human-readable string representations. It supports multiple SV callers including
# Sniffles2, Delly, Manta, and SVIM, as well as small variants (SNPs/indels).

# Example:
#     >>> import pandas as pd
#     >>> factory = VariantEncoderFactory()
#     >>> row = pd.Series({
#     ...     "CHROM": "chr1",
#     ...     "POS": 1000,
#     ...     "REF": "A",
#     ...     "ALT": "<DEL>",
#     ...     "INFO": "SVTYPE=DEL;SVLEN=-500;END=1500"
#     ... })
#     >>> factory.encode_variant(row)
#     'chr1:1000_del_500bp'
# """

# from __future__ import annotations

# import logging
# from dataclasses import dataclass, field
# from typing import TYPE_CHECKING, Protocol

# if TYPE_CHECKING:
#     import pandas as pd
# else:
#     try:
#         import pandas as pd
#     except ImportError:
#         pd = None  # type: ignore[assignment]

# # Set up logging - use loguru if available, otherwise standard logging
# try:
#     from loguru import logger
# except ImportError:
#     logging.basicConfig(level=logging.DEBUG)
#     logger = logging.getLogger(__name__)  # type: ignore[assignment]


# class VariantEncoder(Protocol):
#     """Protocol defining the interface for variant encoders.

#     All variant encoders must implement this protocol to be used with
#     the VariantEncoderFactory. Encoders are responsible for detecting
#     whether they can handle a specific variant type and encoding it
#     into a standardized string format.
#     """

#     def can_encode(self, row: pd.Series) -> bool:
#         """Check if this encoder can handle the given variant row.

#         Args:
#             row: A pandas Series representing one row from a VCF DataFrame.
#                 Expected to contain at minimum: CHROM, POS, REF, ALT columns.
#                 May also contain INFO column for structural variants.

#         Returns:
#             True if this encoder can handle the variant, False otherwise.
#         """
#         ...

#     def encode(self, row: pd.Series) -> str:
#         """Encode the variant into a human-readable string format.

#         Args:
#             row: A pandas Series representing one row from a VCF DataFrame.

#         Returns:
#             Encoded variant string in the format appropriate for the variant type.
#             Typically: "chrom:pos_type_size" for SVs or "chrom:pos_ref_alt" for small variants.
#         """
#         ...


# def parse_info_field(info: str) -> dict[str, str | bool]:
#     """Parse a VCF INFO field into a dictionary.

#     Handles both key=value pairs and flag fields (keys without values).

#     Args:
#         info: The INFO field string from a VCF record, with fields
#             separated by semicolons (e.g., "SVTYPE=DEL;SVLEN=-500;IMPRECISE").

#     Returns:
#         Dictionary mapping INFO keys to their values. Flag fields (without
#         values) are mapped to True.

#     Example:
#         >>> parse_info_field("SVTYPE=DEL;SVLEN=-500;IMPRECISE")
#         {'SVTYPE': 'DEL', 'SVLEN': '-500', 'IMPRECISE': True}
#     """
#     result: dict[str, str | bool] = {}
#     if not info or info == "." or pd.isna(info):
#         return result

#     for item in str(info).split(";"):
#         if not item:
#             continue
#         if "=" in item:
#             key, value = item.split("=", 1)
#             result[key] = value
#         else:
#             result[item] = True

#     return result


# def format_size(bp: int) -> str:
#     """Format a base pair size as a human-readable string.

#     Converts base pair counts to kilobase notation when >= 1000bp.

#     Args:
#         bp: Size in base pairs (should be a positive integer).

#     Returns:
#         Formatted size string (e.g., "500bp" or "78kb").

#     Example:
#         >>> format_size(500)
#         '500bp'
#         >>> format_size(78874)
#         '78kb'
#     """
#     if bp >= 1000:
#         return f"{bp // 1000}kb"
#     return f"{bp}bp"


# def get_svlen_value(svlen_str: str) -> int | None:
#     """Extract the first SVLEN value, handling multi-allelic cases.

#     SVLEN can be comma-separated for multi-allelic SVs (e.g., "-500,-1000").
#     This function extracts the first value and returns its absolute value.

#     Args:
#         svlen_str: The SVLEN value string from the INFO field.

#     Returns:
#         Absolute value of the first SVLEN entry, or None if parsing fails.

#     Example:
#         >>> get_svlen_value("-500")
#         500
#         >>> get_svlen_value("-500,-1000")
#         500
#     """
#     try:
#         first_value = svlen_str.split(",")[0]
#         return abs(int(first_value))
#     except (ValueError, IndexError):
#         return None


# def validate_required_columns(row: pd.Series) -> None:
#     """Validate that required VCF columns are present in the row.

#     Args:
#         row: A pandas Series representing one row from a VCF DataFrame.

#     Raises:
#         KeyError: If any required column (CHROM, POS, REF, ALT) is missing.
#     """
#     required = ["CHROM", "POS", "REF", "ALT"]
#     missing = [col for col in required if col not in row.index]
#     if missing:
#         raise KeyError(f"Missing required VCF columns: {missing}")


# def is_breakend_alt(alt: str) -> bool:
#     """Check if an ALT field represents a breakend (BND) variant.

#     Breakend notation uses brackets to indicate adjacency, e.g.:
#     - ]13:123456]T (join after position on reverse strand)
#     - T[13:123456[ (join before position on forward strand)

#     Args:
#         alt: The ALT field value from a VCF record.

#     Returns:
#         True if the ALT field contains breakend notation.
#     """
#     return bool(alt) and (
#         ("[" in alt or "]" in alt) and not alt.startswith("<")
#     )


# def is_symbolic_alt(alt: str) -> bool:
#     """Check if an ALT field is a symbolic allele (e.g., <DEL>, <INS>).

#     Args:
#         alt: The ALT field value from a VCF record.

#     Returns:
#         True if the ALT field is a symbolic structural variant allele.
#     """
#     return bool(alt) and alt.startswith("<") and alt.endswith(">")


# @dataclass(slots=True, frozen=True)
# class SmallVariantEncoder:
#     """Encoder for small variants (SNPs and small indels).

#     Handles variants that are not structural variants, encoding them
#     in the format: chrom:pos_ref_alt

#     For multi-allelic sites, produces comma-separated encodings for
#     each alternate allele.

#     Example:
#         Input: CHROM=chr1, POS=1000, REF=A, ALT=G
#         Output: "chr1:1000_A_G"

#         Input: CHROM=chr1, POS=1000, REF=A, ALT=G,T
#         Output: "chr1:1000_A_G,chr1:1000_A_T"
#     """

#     def can_encode(self, row: pd.Series) -> bool:
#         """Check if variant is a small variant (not a structural variant).

#         A variant is considered "small" if it has explicit REF/ALT sequences
#         (not symbolic alleles) and no SVTYPE in the INFO field.

#         Args:
#             row: A pandas Series representing one row from a VCF DataFrame.

#         Returns:
#             True if this is a small variant that can be encoded.
#         """
#         alt = str(row.get("ALT", ""))
#         info_str = str(row.get("INFO", ""))
#         info = parse_info_field(info_str)

#         return (
#             not is_symbolic_alt(alt)
#             and not is_breakend_alt(alt)
#             and "SVTYPE" not in info
#             and "END" not in info
#         )

#     def encode(self, row: pd.Series) -> str:
#         """Encode small variant as chrom:pos_ref_alt.

#         Args:
#             row: A pandas Series representing one row from a VCF DataFrame.

#         Returns:
#             Encoded variant string(s). For multi-allelic sites, returns
#             comma-separated encodings.
#         """
#         chrom = str(row["CHROM"])
#         pos = str(row["POS"])
#         ref = str(row["REF"])
#         alts = str(row["ALT"]).split(",")

#         return ",".join([f"{chrom}:{pos}_{ref}_{alt}" for alt in alts])


# @dataclass(slots=True, frozen=True)
# class BreakendEncoder:
#     """Encoder for breakend (BND/translocation) variants.

#     Handles breakend notation where the ALT field contains bracket notation
#     indicating joins between different genomic positions.

#     Example:
#         Input: CHROM=chr1, POS=1000, ALT=]chr2:5000]T, INFO=SVTYPE=BND
#         Output: "chr1:1000_bnd_chr2:5000"
#     """

#     def can_encode(self, row: pd.Series) -> bool:
#         """Check if variant is a breakend variant.

#         Args:
#             row: A pandas Series representing one row from a VCF DataFrame.

#         Returns:
#             True if this is a breakend variant.
#         """
#         alt = str(row.get("ALT", ""))
#         return is_breakend_alt(alt)

#     def encode(self, row: pd.Series) -> str:
#         """Encode breakend variant with partner location.

#         Args:
#             row: A pandas Series representing one row from a VCF DataFrame.

#         Returns:
#             Encoded breakend string showing both breakpoint positions.
#         """
#         chrom = str(row["CHROM"])
#         pos = str(row["POS"])
#         alt = str(row["ALT"])

#         # Extract partner position from breakend notation
#         # Formats: ]chr:pos]seq, [chr:pos[seq, seq]chr:pos], seq[chr:pos[
#         import re
#         match = re.search(r"[\[\]]([^:\[\]]+):(\d+)[\[\]]", alt)
#         if match:
#             partner_chrom = match.group(1)
#             partner_pos = match.group(2)
#             return f"{chrom}:{pos}_bnd_{partner_chrom}:{partner_pos}"

#         return f"{chrom}:{pos}_bnd"


# @dataclass(slots=True, frozen=True)
# class Sniffles2Encoder:
#     """Encoder for Sniffles2 structural variants.

#     Sniffles2 outputs SVs with SVTYPE and SVLEN in the INFO field.
#     This encoder specifically identifies Sniffles2 output by the presence
#     of both SVTYPE and SVLEN fields.

#     Example:
#         Input: INFO=SVTYPE=DEL;SVLEN=-78874
#         Output: "4:143990187_del_78kb"
#     """

#     def can_encode(self, row: pd.Series) -> bool:
#         """Check if variant is from Sniffles2 (has SVTYPE and SVLEN).

#         Args:
#             row: A pandas Series representing one row from a VCF DataFrame.

#         Returns:
#             True if this appears to be a Sniffles2 structural variant.
#         """
#         alt = str(row.get("ALT", ""))
#         info_str = str(row.get("INFO", ""))
#         info = parse_info_field(info_str)

#         # Sniffles2 has SVTYPE and SVLEN, may or may not have symbolic ALT
#         return (
#             "SVTYPE" in info
#             and "SVLEN" in info
#             and not is_breakend_alt(alt)
#         )

#     def encode(self, row: pd.Series) -> str:
#         """Encode as chrom:pos_svtype_size.

#         Args:
#             row: A pandas Series representing one row from a VCF DataFrame.

#         Returns:
#             Encoded SV string (e.g., "4:143990187_del_78kb").
#         """
#         chrom = str(row["CHROM"])
#         pos = str(row["POS"])
#         info = parse_info_field(str(row.get("INFO", "")))

#         svtype = str(info.get("SVTYPE", "unknown")).lower()
#         svlen = get_svlen_value(str(info.get("SVLEN", "0")))

#         if svlen is not None and svlen > 0:
#             size_str = format_size(svlen)
#             return f"{chrom}:{pos}_{svtype}_{size_str}"

#         # Fallback if SVLEN parsing fails
#         return f"{chrom}:{pos}_{svtype}"


# @dataclass(slots=True, frozen=True)
# class SvimEncoder:
#     """Encoder for SVIM structural variants.

#     SVIM (Structural Variant Identification using Mapped long reads) outputs
#     SVs with SVTYPE, END, and optionally SVLEN in the INFO field. SVIM uses
#     symbolic ALT alleles like <DEL>, <INS>, <DUP>, <DUP:TANDEM>, <DUP:INT>,
#     <INV>, and <BND>.

#     SVIM is distinguished by having SUPPORT field in INFO (read support count)
#     along with SVTYPE.

#     Example:
#         Input: ALT=<DEL>, INFO=SVTYPE=DEL;END=1500;SVLEN=-500;SUPPORT=10
#         Output: "chr1:1000_del_500bp"
#     """

#     def can_encode(self, row: pd.Series) -> bool:
#         """Check if variant is from SVIM (has SVTYPE and SUPPORT).

#         SVIM includes a SUPPORT field indicating the number of reads
#         supporting the variant call, which distinguishes it from other callers.

#         Args:
#             row: A pandas Series representing one row from a VCF DataFrame.

#         Returns:
#             True if this appears to be a SVIM structural variant.
#         """
#         alt = str(row.get("ALT", ""))
#         info_str = str(row.get("INFO", ""))
#         info = parse_info_field(info_str)

#         return (
#             "SVTYPE" in info
#             and "SUPPORT" in info
#             and not is_breakend_alt(alt)
#         )

#     def encode(self, row: pd.Series) -> str:
#         """Encode SVIM variant as chrom:pos_svtype_size.

#         Args:
#             row: A pandas Series representing one row from a VCF DataFrame.

#         Returns:
#             Encoded SV string (e.g., "chr1:1000_del_500bp").
#         """
#         chrom = str(row["CHROM"])
#         pos = int(row["POS"])
#         info = parse_info_field(str(row.get("INFO", "")))

#         # Handle subtypes like DUP:TANDEM -> dup_tandem
#         svtype_raw = str(info.get("SVTYPE", "unknown"))
#         svtype = svtype_raw.lower().replace(":", "_")

#         # Try SVLEN first, then calculate from END
#         svlen: int | None = None
#         if "SVLEN" in info:
#             svlen = get_svlen_value(str(info["SVLEN"]))

#         if svlen is None and "END" in info:
#             try:
#                 end = int(info["END"])
#                 svlen = abs(end - pos)
#             except (ValueError, TypeError):
#                 pass

#         if svlen is not None and svlen > 0:
#             size_str = format_size(svlen)
#             return f"{chrom}:{pos}_{svtype}_{size_str}"

#         return f"{chrom}:{pos}_{svtype}"


# @dataclass(slots=True, frozen=True)
# class DellyEncoder:
#     """Encoder for Delly structural variants.

#     Delly outputs SVs with SVTYPE and END but typically not SVLEN.
#     This encoder identifies Delly output by the presence of SVTYPE
#     and END without SVLEN.

#     Example:
#         Input: INFO=SVTYPE=DEL;END=144069060
#         Output: "4:143990187_del_78kb"
#     """

#     def can_encode(self, row: pd.Series) -> bool:
#         """Check if variant is from Delly (has SVTYPE and END but not SVLEN).

#         Args:
#             row: A pandas Series representing one row from a VCF DataFrame.

#         Returns:
#             True if this appears to be a Delly structural variant.
#         """
#         alt = str(row.get("ALT", ""))
#         info_str = str(row.get("INFO", ""))
#         info = parse_info_field(info_str)

#         return (
#             "SVTYPE" in info
#             and "END" in info
#             and "SVLEN" not in info
#             and "SUPPORT" not in info  # Exclude SVIM
#             and not is_breakend_alt(alt)
#         )

#     def encode(self, row: pd.Series) -> str:
#         """Encode Delly variant as chrom:pos_svtype_size.

#         Args:
#             row: A pandas Series representing one row from a VCF DataFrame.

#         Returns:
#             Encoded SV string (e.g., "4:143990187_del_78kb").
#         """
#         chrom = str(row["CHROM"])
#         pos = int(row["POS"])
#         info = parse_info_field(str(row.get("INFO", "")))

#         svtype = str(info.get("SVTYPE", "unknown")).lower()

#         try:
#             end = int(info["END"])
#             svlen = abs(end - pos)
#             size_str = format_size(svlen)
#             return f"{chrom}:{pos}_{svtype}_{size_str}"
#         except (KeyError, ValueError, TypeError):
#             return f"{chrom}:{pos}_{svtype}"


# @dataclass(slots=True, frozen=True)
# class MantaEncoder:
#     """Encoder for Manta structural variants.

#     Manta outputs SVs with symbolic ALT alleles (<DEL>, <INS>, etc.)
#     and SVTYPE in INFO. It may have either SVLEN or END (or both).

#     This encoder is used as a fallback for symbolic SV alleles that
#     don't match more specific caller signatures.

#     Example:
#         Input: ALT=<DEL>, INFO=SVTYPE=DEL;SVLEN=-5000
#         Output: "chr1:1000_del_5kb"
#     """

#     def can_encode(self, row: pd.Series) -> bool:
#         """Check if variant is a symbolic SV (generic handler for Manta-style).

#         Args:
#             row: A pandas Series representing one row from a VCF DataFrame.

#         Returns:
#             True if this is a symbolic structural variant.
#         """
#         alt = str(row.get("ALT", ""))
#         info_str = str(row.get("INFO", ""))
#         info = parse_info_field(info_str)

#         return (
#             is_symbolic_alt(alt)
#             and "SVTYPE" in info
#             and not is_breakend_alt(alt)
#         )

#     def encode(self, row: pd.Series) -> str:
#         """Encode Manta-style variant as chrom:pos_svtype_size.

#         Args:
#             row: A pandas Series representing one row from a VCF DataFrame.

#         Returns:
#             Encoded SV string (e.g., "chr1:1000_del_5kb").
#         """
#         chrom = str(row["CHROM"])
#         pos = int(row["POS"])
#         info = parse_info_field(str(row.get("INFO", "")))

#         svtype = str(info.get("SVTYPE", "unknown")).lower()

#         # Try SVLEN first
#         svlen: int | None = None
#         if "SVLEN" in info:
#             svlen = get_svlen_value(str(info["SVLEN"]))

#         # Fall back to END if no SVLEN
#         if svlen is None and "END" in info:
#             try:
#                 end = int(info["END"])
#                 svlen = abs(end - pos)
#             except (ValueError, TypeError):
#                 pass

#         if svlen is not None and svlen > 0:
#             size_str = format_size(svlen)
#             return f"{chrom}:{pos}_{svtype}_{size_str}"

#         return f"{chrom}:{pos}_{svtype}"


# @dataclass(slots=True, frozen=False)
# class VariantEncoderFactory:
#     """Factory to manage and apply variant encoders.

#     Encoders are tried in order until one matches. The order matters:
#     more specific encoders (like SVIM, Sniffles2) should come before
#     generic ones (like MantaEncoder).

#     The default encoder order is:
#     1. BreakendEncoder - for BND/translocation variants
#     2. SvimEncoder - SVIM-specific (has SUPPORT field)
#     3. Sniffles2Encoder - Sniffles2-specific (has SVLEN)
#     4. DellyEncoder - Delly-specific (has END but not SVLEN)
#     5. MantaEncoder - generic symbolic SV handler
#     6. SmallVariantEncoder - fallback for SNPs/indels

#     Attributes:
#         encoders: List of encoder instances in priority order.

#     Example:
#         >>> factory = VariantEncoderFactory()
#         >>> encoded = factory.encode_variant(vcf_row)

#         # Add a custom encoder with highest priority
#         >>> factory.add_encoder(MyCustomEncoder(), priority=0)
#     """

#     encoders: list[VariantEncoder] = field(
#         default_factory=lambda: [
#             BreakendEncoder(),
#             SvimEncoder(),
#             Sniffles2Encoder(),
#             DellyEncoder(),
#             MantaEncoder(),
#             SmallVariantEncoder(),  # Keep this last as fallback
#         ]
#     )

#     def add_encoder(self, encoder: VariantEncoder, priority: int = 0) -> None:
#         """Add a custom encoder at specified priority position.

#         Args:
#             encoder: The encoder instance to add. Must implement the
#                 VariantEncoder protocol.
#             priority: Position in the encoder list where 0 is highest
#                 priority (tried first). Use len(encoders) to add at end.

#         Example:
#             >>> factory = VariantEncoderFactory()
#             >>> factory.add_encoder(MyEncoder(), priority=0)  # Highest priority
#             >>> factory.add_encoder(MyEncoder(), priority=len(factory.encoders))  # Lowest
#         """
#         self.encoders.insert(priority, encoder)

#     def encode_variant(self, row: pd.Series) -> str:
#         """Encode a variant using the first matching encoder.

#         Iterates through encoders in priority order, using the first
#         encoder that can handle the variant.

#         Args:
#             row: A pandas Series representing one row from the VCF DataFrame.
#                 Must contain CHROM, POS, REF, ALT columns. May contain INFO.

#         Returns:
#             Encoded variant string in the format appropriate for the variant type.

#         Raises:
#             KeyError: If required columns (CHROM, POS, REF, ALT) are missing.

#         Example:
#             >>> row = pd.Series({"CHROM": "chr1", "POS": 1000, "REF": "A",
#             ...                  "ALT": "G", "INFO": "."})
#             >>> factory.encode_variant(row)
#             'chr1:1000_A_G'
#         """
#         validate_required_columns(row)

#         for encoder in self.encoders:
#             if encoder.can_encode(row):
#                 logger.debug(f"Using encoder: {encoder.__class__.__name__}")
#                 return encoder.encode(row)

#         # Ultimate fallback if no encoder matches (should not happen with SmallVariantEncoder)
#         logger.warning(
#             f"No encoder matched variant at {row['CHROM']}:{row['POS']}. "
#             "Using raw format."
#         )
#         return f"{row['CHROM']}:{row['POS']}_{row['REF']}_{row['ALT']}"
    
# from dataclasses import dataclass, field
# from typing import Protocol
# import pandas as pd
# from loguru import logger
# from icecream import ic

# class VariantEncoder(Protocol):
#     """Protocol defining the interface for variant encoders."""

#     def can_encode(self, row: pd.Series) -> bool:
#         """Check if this encoder can handle the given variant row."""
#         ...

#     def encode(self, row: pd.Series) -> str:
#         """Encode the variant into the desired format."""
#         ...


# @dataclass(slots=True, frozen=True)
# class SmallVariantEncoder:
#     """Encoder for small variants (SNPs, small indels)."""

#     def can_encode(self, row: pd.Series) -> bool:
#         """Check if variant is a small variant (not SV)."""
#         ref = str(row["REF"])
#         alt = str(row["ALT"])
#         info = str(row.get("INFO", ""))

#         # Check if it's NOT a structural variant
#         return (
#             not alt.startswith("<") and not alt.endswith(">") and "SVTYPE=" not in info
#         )

#     def encode(self, row: pd.Series) -> str:
#         """Encode small variant as chrom:pos_ref_alt."""
#         chrom = str(row["CHROM"])
#         pos = str(row["POS"])
#         ref = str(row["REF"])
#         alts = str(row["ALT"]).split(",")

#         return ",".join([f"{chrom}:{pos}_{ref}_{alt}" for alt in alts])


# @dataclass(slots=True, frozen=True)
# class Sniffles2Encoder:
#     """Encoder for Sniffles2 structural variants.

#     Sniffles2 outputs SVs with SVTYPE and SVLEN in INFO field.
#     Example: SVTYPE=DEL;SVLEN=-78874
#     """

#     def can_encode(self, row: pd.Series) -> bool:
#         """Check if variant is from Sniffles2 (has SVTYPE and SVLEN)."""
#         info = str(row.get("INFO", ""))
#         return "SVTYPE=" in info and "SVLEN=" in info

#     def encode(self, row: pd.Series) -> str:
#         """Encode as chrom:pos_svtype_size (e.g., 4:143990187_del_78kb)."""
#         chrom = str(row["CHROM"])
#         pos = str(row["POS"])
#         info = str(row["INFO"])

#         # Extract SVTYPE
#         svtype = None
#         for field in info.split(";"):
#             if field.startswith("SVTYPE="):
#                 svtype = field.split("=")[1].lower()
#                 break

#         # Extract SVLEN
#         svlen = None
#         for field in info.split(";"):
#             if field.startswith("SVLEN="):
#                 svlen_str = field.split("=")[1]
#                 svlen = abs(int(svlen_str))  # Take absolute value
#                 break

#         if svtype and svlen:
#             # Convert to kb format
#             kb = svlen / 1000
#             if kb >= 1:
#                 size_str = f"{int(kb)}kb"
#             else:
#                 size_str = f"{svlen}bp"

#             return f"{chrom}:{pos}_{svtype}_{size_str}"

#         # Fallback if parsing fails
#         return f"{chrom}:{pos}_{row['REF']}_{row['ALT']}"


# @dataclass(slots=True, frozen=True)
# class DellyEncoder:
#     """Encoder for Delly structural variants.

#     Delly outputs SVs with SVTYPE and END (but not SVLEN).
#     Example: SVTYPE=DEL;END=144069060
#     """

#     def can_encode(self, row: pd.Series) -> bool:
#         """Check if variant is from Delly (has SVTYPE and END but not SVLEN)."""
#         info = str(row.get("INFO", ""))
#         return "SVTYPE=" in info and "END=" in info and "SVLEN=" not in info

#     def encode(self, row: pd.Series) -> str:
#         """Encode as chrom:pos_svtype_size."""
#         chrom = str(row["CHROM"])
#         pos = int(row["POS"])
#         info = str(row["INFO"])

#         # Extract SVTYPE
#         svtype = None
#         for field in info.split(";"):
#             if field.startswith("SVTYPE="):
#                 svtype = field.split("=")[1].lower()
#                 break

#         # Extract END position
#         end = None
#         for field in info.split(";"):
#             if field.startswith("END="):
#                 end = int(field.split("=")[1])
#                 break

#         if svtype and end:
#             svlen = abs(end - pos)
#             kb = svlen / 1000
#             if kb >= 1:
#                 size_str = f"{int(kb)}kb"
#             else:
#                 size_str = f"{svlen}bp"

#             return f"{chrom}:{pos}_{svtype}_{size_str}"

#         # Fallback
#         return f"{chrom}:{pos}_{row['REF']}_{row['ALT']}"


# @dataclass(slots=True, frozen=True)
# class MantaEncoder:
#     """Encoder for Manta structural variants.

#     Manta outputs SVs with SVTYPE and may have SVLEN or END.
#     Example: SVTYPE=DEL;SVLEN=-5000 or SVTYPE=DEL;END=144069060
#     """

#     def can_encode(self, row: pd.Series) -> bool:
#         """Check if variant is from Manta (bracketed ALT with SVTYPE)."""
#         info = str(row.get("INFO", ""))
#         alt = str(row.get("ALT", ""))
#         return (alt.startswith("<") and alt.endswith(">")) and "SVTYPE=" in info

#     def encode(self, row: pd.Series) -> str:
#         """Encode as chrom:pos_svtype_size."""
#         chrom = str(row["CHROM"])
#         pos = str(row["POS"])
#         info = str(row["INFO"])

#         # Extract SVTYPE
#         svtype = None
#         for field in info.split(";"):
#             if field.startswith("SVTYPE="):
#                 svtype = field.split("=")[1].lower()
#                 break

#         # Try to extract SVLEN (Manta may have it)
#         svlen = None
#         for field in info.split(";"):
#             if field.startswith("SVLEN="):
#                 svlen_str = field.split("=")[1]
#                 svlen = abs(int(svlen_str))
#                 break

#         # If no SVLEN, try END
#         if svlen is None:
#             end = None
#             for field in info.split(";"):
#                 if field.startswith("END="):
#                     end = int(field.split("=")[1])
#                     break
#             if end:
#                 svlen = abs(end - int(pos))

#         if svtype and svlen:
#             kb = svlen / 1000
#             if kb >= 1:
#                 size_str = f"{int(kb)}kb"
#             else:
#                 size_str = f"{svlen}bp"

#             return f"{chrom}:{pos}_{svtype}_{size_str}"

#         # Fallback
#         return f"{chrom}:{pos}_{row['REF']}_{row['ALT']}"


# @dataclass(slots=True, frozen=False)
# class VariantEncoderFactory:
#     """Factory to manage and apply variant encoders.

#     Encoders are tried in order until one matches. Order matters:
#     more specific encoders should come before general ones.
#     """

#     encoders: list[VariantEncoder] = field(
#         default_factory=lambda: [
#             Sniffles2Encoder(),
#             DellyEncoder(),
#             MantaEncoder(),
#             SmallVariantEncoder(),  # Keep this last as fallback
#         ]
#     )

#     def add_encoder(self, encoder: VariantEncoder, priority: int = 0) -> None:
#         """Add a custom encoder at specified priority (lower index = higher priority).

#         Args:
#             encoder: The encoder instance to add.
#             priority: Position in the encoder list (0 = highest priority).
#         """
#         self.encoders.insert(priority, encoder)

#     def encode_variant(self, row: pd.Series) -> str:
#         """Encode a variant using the first matching encoder.

#         Args:
#             row: A pandas Series representing one row from the VCF DataFrame.

#         Returns:
#             Encoded variant string in the format appropriate for the variant type.
#         """
#         for encoder in self.encoders:
#             if encoder.can_encode(row):
#                 ic(encoder)
#                 return encoder.encode(row)
#         logger.warning("Variants encoded in an unsupported way")
#         # Ultimate fallback if no encoder matches
#         return f"{row['CHROM']}:{row['POS']}_{row['REF']}_{row['ALT']}"
