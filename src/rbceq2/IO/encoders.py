from dataclasses import dataclass, field
from typing import Protocol
import pandas as pd
from loguru import logger


class VariantEncoder(Protocol):
    """Protocol defining the interface for variant encoders."""

    def can_encode(self, row: pd.Series) -> bool:
        """Check if this encoder can handle the given variant row."""
        ...

    def encode(self, row: pd.Series) -> str:
        """Encode the variant into the desired format."""
        ...


@dataclass(slots=True, frozen=True)
class SmallVariantEncoder:
    """Encoder for small variants (SNPs, small indels)."""

    def can_encode(self, row: pd.Series) -> bool:
        """Check if variant is a small variant (not SV)."""
        ref = str(row["REF"])
        alt = str(row["ALT"])
        info = str(row.get("INFO", ""))

        # Check if it's NOT a structural variant
        return (
            not alt.startswith("<") and not alt.endswith(">") and "SVTYPE=" not in info
        )

    def encode(self, row: pd.Series) -> str:
        """Encode small variant as chrom:pos_ref_alt."""
        chrom = str(row["CHROM"])
        pos = str(row["POS"])
        ref = str(row["REF"])
        alts = str(row["ALT"]).split(",")

        return ",".join([f"{chrom}:{pos}_{ref}_{alt}" for alt in alts])


@dataclass(slots=True, frozen=True)
class Sniffles2Encoder:
    """Encoder for Sniffles2 structural variants.

    Sniffles2 outputs SVs with SVTYPE and SVLEN in INFO field.
    Example: SVTYPE=DEL;SVLEN=-78874
    """

    def can_encode(self, row: pd.Series) -> bool:
        """Check if variant is from Sniffles2 (has SVTYPE and SVLEN)."""
        info = str(row.get("INFO", ""))
        return "SVTYPE=" in info and "SVLEN=" in info

    def encode(self, row: pd.Series) -> str:
        """Encode as chrom:pos_svtype_size (e.g., 4:143990187_del_78kb)."""
        chrom = str(row["CHROM"])
        pos = str(row["POS"])
        info = str(row["INFO"])

        # Extract SVTYPE
        svtype = None
        for field in info.split(";"):
            if field.startswith("SVTYPE="):
                svtype = field.split("=")[1].lower()
                break

        # Extract SVLEN
        svlen = None
        for field in info.split(";"):
            if field.startswith("SVLEN="):
                svlen_str = field.split("=")[1]
                svlen = abs(int(svlen_str))  # Take absolute value
                break

        if svtype and svlen:
            # Convert to kb format
            kb = svlen / 1000
            if kb >= 1:
                size_str = f"{int(kb)}kb"
            else:
                size_str = f"{svlen}bp"

            return f"{chrom}:{pos}_{svtype}_{size_str}"

        # Fallback if parsing fails
        return f"{chrom}:{pos}_{row['REF']}_{row['ALT']}"


@dataclass(slots=True, frozen=True)
class DellyEncoder:
    """Encoder for Delly structural variants.

    Delly outputs SVs with SVTYPE and END (but not SVLEN).
    Example: SVTYPE=DEL;END=144069060
    """

    def can_encode(self, row: pd.Series) -> bool:
        """Check if variant is from Delly (has SVTYPE and END but not SVLEN)."""
        info = str(row.get("INFO", ""))
        return "SVTYPE=" in info and "END=" in info and "SVLEN=" not in info

    def encode(self, row: pd.Series) -> str:
        """Encode as chrom:pos_svtype_size."""
        chrom = str(row["CHROM"])
        pos = int(row["POS"])
        info = str(row["INFO"])

        # Extract SVTYPE
        svtype = None
        for field in info.split(";"):
            if field.startswith("SVTYPE="):
                svtype = field.split("=")[1].lower()
                break

        # Extract END position
        end = None
        for field in info.split(";"):
            if field.startswith("END="):
                end = int(field.split("=")[1])
                break

        if svtype and end:
            svlen = abs(end - pos)
            kb = svlen / 1000
            if kb >= 1:
                size_str = f"{int(kb)}kb"
            else:
                size_str = f"{svlen}bp"

            return f"{chrom}:{pos}_{svtype}_{size_str}"

        # Fallback
        return f"{chrom}:{pos}_{row['REF']}_{row['ALT']}"


@dataclass(slots=True, frozen=True)
class MantaEncoder:
    """Encoder for Manta structural variants.

    Manta outputs SVs with SVTYPE and may have SVLEN or END.
    Example: SVTYPE=DEL;SVLEN=-5000 or SVTYPE=DEL;END=144069060
    """

    def can_encode(self, row: pd.Series) -> bool:
        """Check if variant is from Manta (bracketed ALT with SVTYPE)."""
        info = str(row.get("INFO", ""))
        alt = str(row.get("ALT", ""))
        return (alt.startswith("<") and alt.endswith(">")) and "SVTYPE=" in info

    def encode(self, row: pd.Series) -> str:
        """Encode as chrom:pos_svtype_size."""
        chrom = str(row["CHROM"])
        pos = str(row["POS"])
        info = str(row["INFO"])

        # Extract SVTYPE
        svtype = None
        for field in info.split(";"):
            if field.startswith("SVTYPE="):
                svtype = field.split("=")[1].lower()
                break

        # Try to extract SVLEN (Manta may have it)
        svlen = None
        for field in info.split(";"):
            if field.startswith("SVLEN="):
                svlen_str = field.split("=")[1]
                svlen = abs(int(svlen_str))
                break

        # If no SVLEN, try END
        if svlen is None:
            end = None
            for field in info.split(";"):
                if field.startswith("END="):
                    end = int(field.split("=")[1])
                    break
            if end:
                svlen = abs(end - int(pos))

        if svtype and svlen:
            kb = svlen / 1000
            if kb >= 1:
                size_str = f"{int(kb)}kb"
            else:
                size_str = f"{svlen}bp"

            return f"{chrom}:{pos}_{svtype}_{size_str}"

        # Fallback
        return f"{chrom}:{pos}_{row['REF']}_{row['ALT']}"


@dataclass(slots=True, frozen=False)
class VariantEncoderFactory:
    """Factory to manage and apply variant encoders.

    Encoders are tried in order until one matches. Order matters:
    more specific encoders should come before general ones.
    """

    encoders: list[VariantEncoder] = field(
        default_factory=lambda: [
            Sniffles2Encoder(),
            DellyEncoder(),
            MantaEncoder(),
            SmallVariantEncoder(),  # Keep this last as fallback
        ]
    )

    def add_encoder(self, encoder: VariantEncoder, priority: int = 0) -> None:
        """Add a custom encoder at specified priority (lower index = higher priority).

        Args:
            encoder: The encoder instance to add.
            priority: Position in the encoder list (0 = highest priority).
        """
        self.encoders.insert(priority, encoder)

    def encode_variant(self, row: pd.Series) -> str:
        """Encode a variant using the first matching encoder.

        Args:
            row: A pandas Series representing one row from the VCF DataFrame.

        Returns:
            Encoded variant string in the format appropriate for the variant type.
        """
        for encoder in self.encoders:
            if encoder.can_encode(row):
                return encoder.encode(row)
        logger.warning("Variants encoded in an unsupported way")
        # Ultimate fallback if no encoder matches
        return f"{row['CHROM']}:{row['POS']}_{row['REF']}_{row['ALT']}"
