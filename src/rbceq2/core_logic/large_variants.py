from __future__ import annotations

from pathlib import Path

from collections import defaultdict
import gzip
import io
from dataclasses import dataclass, field
from typing import Iterator, Protocol, runtime_checkable

from loguru import logger
import csv

import re
from typing import Iterable

from icecream import ic


@dataclass(slots=True, frozen=True)
class SvDef:
    """A structural variant definition parsed from a database token.

    Supports mixed encodings commonly found in allele-definition tables:
      - ``"<pos>_del_<59kb>"``     → type=DEL, length from unitized number.
      - ``"<pos>_<REFseq>_<ALTseq>"`` → type inferred from ``len(ALT)-len(REF)``.
      - ``"<pos>_ins_<12kb>"``     → type=INS.
      - ``"<pos>_dup_<Xkb|Xbp>"``  → type=DUP.

    Attributes:
        chrom (str): Chromosome name (no ``chr`` prefix).
        pos (int): Left-most 1-based position from the DB entry (approximate).
        svtype (str): SV type (e.g., DEL, INS, DUP, INV, CNV, INDEL).
        length (int): Absolute length in bp (approximate for unitized inputs).
        raw (str): Raw token from the DB for traceability.
        id (str): Allele/db identifier if available.
        tol_pos (int): Allowed ± positional tolerance in bp when matching.
        tol_len (int): Allowed absolute length tolerance in bp.
        tol_ratio (float): Allowed fractional tolerance on length for large SVs.
    """

    chrom: str
    pos: int
    svtype: str
    length: int
    raw: str
    id: str = "-"
    tol_pos: int = 25_000
    tol_len: int = 10_000
    tol_ratio: float = 0.25

    @property
    def interval(self) -> tuple[int, int]:
        """Compute an approximate affected interval.

        For deletions/duplications/inversions, returns a span of ``length``.
        For insertions, returns a 1-bp span starting at ``pos``.

        Returns:
            tuple[int, int]: Inclusive 1-based start and end.
        """
        if self.svtype in {"DEL", "INV", "CNV", "DUP"}:
            return (self.pos, self.pos + max(self.length, 1))
        return (self.pos, self.pos + 1)


_UNIT = re.compile(r"^\s*(\d+(?:\.\d+)?)\s*(bp|b|kb|mb)\s*$", re.I)


def _parse_length_unit(tok: str) -> int:
    """Parse a unitized length token into base pairs.

    Args:
        tok (str): Token like ``"59kb"``, ``"400bp"``, or a bare integer.

    Returns:
        int: Length in base pairs.

    Raises:
        ValueError: If the token cannot be interpreted.
    """
    m = _UNIT.match(tok)
    if not m:
        if tok.isdigit():
            return int(tok)
        raise ValueError(f"Unrecognized length unit: {tok}")
    val = float(m.group(1))
    unit = m.group(2).lower()
    if unit in {"bp", "b"}:
        return int(round(val))
    if unit == "kb":
        return int(round(val * 1_000))
    if unit == "mb":
        return int(round(val * 1_000_000))
    raise ValueError(f"Unhandled unit {unit}")


def parse_db_token(chrom: str, token: str, db_id: str = "-") -> SvDef | None:
    """Parse a single DB token into an :class:`SvDef`.

    Handles both sequence and word forms. Returns ``None`` if the token does
    not encode a structural event.

    Args:
        chrom (str): Chromosome (no ``chr`` prefix).
        token (str): DB token, e.g., ``"25272547_del_59kb"`` or
            ``"126690214_<REFseq>_<ALTseq>"``.
        db_id (str): Optional allele identifier for traceability.

    Returns:
        SvDef | None: Parsed definition or ``None`` when unsupported.
    """
    raw = token.strip()
    parts = raw.split("_")
    if not parts or not parts[0].isdigit():
        return None

    pos = int(parts[0])

    # Sequence form: <pos>_<REF>_<ALT>
    if len(parts) >= 3 and (set(parts[1]) <= set("ACGTNacgtn") or len(parts[1]) > 50):
        ref = parts[1]
        alt = parts[2]
        delta = len(alt) - len(ref)
        svtype = "DEL" if delta < 0 else ("INS" if delta > 0 else "INDEL")
        return SvDef(chrom=chrom, pos=pos, svtype=svtype, length=abs(delta), raw=raw, id=db_id)

    # Word form: <pos>_<type>_<len>
    if len(parts) >= 3 and parts[1].lower() in {"del", "dup", "ins", "inv", "cnv"}:
        svt = parts[1].upper()
        ln = _parse_length_unit(parts[2])
        return SvDef(chrom=chrom, pos=pos, svtype=svt, length=ln, raw=raw, id=db_id)

    return None


@dataclass(slots=True)
class MatchResult:
    """A single best match between a DB definition and a VCF event.

    Attributes:
        db (SvDef): Database definition that was matched.
        vcf (SvEvent): VCF structural variant event selected as the match.
        score (float): Composite score (lower is better).
        pos_delta (int): Absolute difference in positions (bp).
        len_delta (int): Absolute difference in lengths (bp).
    """

    db: SvDef
    vcf: SvEvent
    score: float
    pos_delta: int
    len_delta: int


@dataclass(slots=True)
class SvMatcher:
    """Fuzzy matcher for DB SV definitions vs VCF events.

    Implements tolerant matching suitable for imprecise breakpoints, mixing
    absolute and fractional length tolerances, with optional interval overlap
    reinforcement.

    Attributes:
        pos_bonus_overlap (float): Score bonus subtracted when intervals overlap.
        require_same_type (bool): If True, require compatible SV types.
    """

    pos_bonus_overlap: float = 0.5
    require_same_type: bool = True

    def compatible(self, db: SvDef, ev: SvEvent) -> bool:
        """Check basic type compatibility between DB and VCF events.

        Args:
            db (SvDef): DB definition.
            ev (SvEvent): VCF event.

        Returns:
            bool: True if types are considered compatible.
        """
        if not self.require_same_type:
            return True
        dt = db.svtype
        et = ev.svtype
        if dt == "INDEL":
            return et in {"DEL", "INS", "INDEL"}
        if et == "INDEL":
            return dt in {"DEL", "INS", "INDEL"}
        return dt == et

    def _length(self, ev: SvEvent) -> int:
        """Compute absolute length for a VCF event.

        Args:
            ev (SvEvent): VCF event.

        Returns:
            int: Absolute length (prefers SVLEN, else end-pos).
        """
        return abs(ev.svlen or ev.size)


    def _intervals_overlap(self, db: SvDef, ev: SvEvent) -> bool:
        """Test reciprocal overlap between DB and VCF intervals.

        Requires at least 10% overlap in both directions.

        Args:
            db (SvDef): DB definition.
            ev (SvEvent): VCF event.

        Returns:
            bool: True if intervals reciprocally overlap.
        """
        db_s, db_e = db.interval
        ev_s, ev_e = ev.pos, ev.end
        if ev_e < db_s or ev_s > db_e:
            return False
        inter = min(db_e, ev_e) - max(db_s, ev_s)
        if inter <= 0:
            return False
        db_len = max(1, db_e - db_s)
        ev_len = max(1, ev_e - ev_s)
        return (inter / db_len) >= 0.1 and (inter / ev_len) >= 0.1

    def match(self, db_defs: Iterable[SvDef], events: Iterable[SvEvent]) -> list[MatchResult]:
        """Match DB definitions to VCF events and return best hits per DB token.

        Args:
            db_defs (Iterable[SvDef]): Iterable of parsed DB definitions.
            events (Iterable[SvEvent]): Iterable of VCF SV events.

        Returns:
            list[MatchResult]: Best match per DB token, ordered by
            (chrom, pos, score).
        """
        results: list[MatchResult] = []
        ev_by_chrom: dict[str, list[SvEvent]] = {}
        for ev in events:
            ev_by_chrom.setdefault(ev.chrom, []).append(ev)

        for db in db_defs:
            for ev in ev_by_chrom.get(db.chrom, []):
                if not self.compatible(db, ev):
                    continue
                s, pd, ld = self.score(db, ev)
                if s != float("inf"):
                    results.append(MatchResult(db=db, vcf=ev, score=s, pos_delta=pd, len_delta=ld))

        # Keep best per (allele id, raw token, chrom)
        best: dict[tuple[str, str, str], MatchResult] = {}
        for r in results:
            key = (r.db.id, r.db.raw, r.db.chrom)
            if key not in best or r.score < best[key].score:
                best[key] = r

        return sorted(best.values(), key=lambda r: (r.db.chrom, r.db.pos, r.score))
    
    
    def _adaptive_pos_tol(self, db: SvDef, ev: SvEvent, overlap: bool) -> int:
        """Compute a size-aware positional tolerance.

        Uses the larger of:
          * a small absolute floor (guards tiny SVs),
          * a fraction of the SV length (looser for bigger SVs),
          * any reported CI width (CIPOS/CIEND) if present.
        Becomes more permissive when DB/VCF intervals already overlap.

        Args:
            db (SvDef): DB definition.
            ev (SvEvent): VCF event.
            overlap (bool): Whether DB and VCF intervals overlap.

        Returns:
            int: Allowed ±bp shift for breakpoint comparisons.
        """
        L = max(1, max(db.length, self._length(ev)))
        ci_pos = abs(ev.cipos.left) + abs(ev.cipos.right)
        ci_end = abs(ev.ciend.left) + abs(ev.ciend.right)
        ci_span = max(ci_pos, ci_end)

        POS_FLOOR = 200          # tiny SVs
        POS_FRAC = 0.50          # no-overlap: 50% of size
        POS_FRAC_OVERLAP = 0.75  # overlap: 75% of size
        POS_CAP = 50_000

        frac = POS_FRAC_OVERLAP if overlap else POS_FRAC
        tol = max(POS_FLOOR, int(frac * L), ci_span)
        return min(tol, POS_CAP)

    def _adaptive_len_tol(self, db: SvDef, ev: SvEvent, overlap: bool) -> int:
        """Compute a size-aware length tolerance (HARD gate).

        Uses the larger of:
          * a small absolute floor (for measurement noise),
          * a fraction of the larger of DB/VCF lengths,
          * (optionally) a slightly higher fraction if intervals overlap.

        Args:
            db (SvDef): DB definition.
            ev (SvEvent): VCF event.
            overlap (bool): Whether DB and VCF intervals overlap.

        Returns:
            int: Allowed absolute difference in length (bp).
        """
        Ldb = max(1, db.length)
        Lev = max(1, self._length(ev))
        L = max(Ldb, Lev)

        LEN_FLOOR = 50            # tolerate small caller jitter
        LEN_FRAC = 0.35           # no-overlap: 35% of size
        LEN_FRAC_OVERLAP = 0.50   # overlap: 50% of size
        LEN_CAP = 100_000         # safety cap

        frac = LEN_FRAC_OVERLAP if overlap else LEN_FRAC
        tol = max(LEN_FLOOR, int(frac * L))
        return min(tol, LEN_CAP)

    def score(self, db: SvDef, ev: SvEvent) -> tuple[float, int, int]:
        """Score a DB/VCF pair using adaptive POS and LEN tolerances.

        Position:
            - Adaptive tolerance; if exceeded and there is no overlap → reject.
        Length:
            - Adaptive tolerance; if exceeded → reject (hard), even if overlapping.

        Args:
            db (SvDef): DB definition.
            ev (SvEvent): VCF event.

        Returns:
            tuple[float, int, int]: (score, Δpos, Δlen). ``score`` is ``inf`` if rejected.
        """
        pos_delta = abs(ev.pos - db.pos)
        ev_len = self._length(ev)
        len_delta = abs(ev_len - db.length)

        ov = self._intervals_overlap(db, ev)
        pos_tol = self._adaptive_pos_tol(db, ev, overlap=ov)
        len_tol = self._adaptive_len_tol(db, ev, overlap=ov)

        # Gates
        if pos_delta > pos_tol and not ov:
            return (float("inf"), pos_delta, len_delta)
        if len_delta > len_tol:
            return (float("inf"), pos_delta, len_delta)

        # Score normalized by adaptive tolerances
        s = (pos_delta / (pos_tol + 1)) + (len_delta / (len_tol + 1))
        if ov:
            s -= self.pos_bonus_overlap

        return (max(s, 0.0), pos_delta, len_delta)



def _sniff_delimiter(path: Path) -> str:
    """Guess file delimiter using csv.Sniffer, fallback to tab.

    Args:
        path (Path): TSV/CSV path.

    Returns:
        str: Detected delimiter.
    """
    with open(path, "rt", encoding="utf-8-sig", newline="") as fh:
        sample = fh.read(8192)
    try:
        dialect = csv.Sniffer().sniff(sample, delimiters=[",", "\t", ";", "|"])
        return dialect.delimiter
    except Exception:
        return "\t"


def _ci_lookup(names: list[str]) -> dict[str, str]:
    """Case-insensitive header map: lower->original."""
    return {n.lower(): n for n in names}


def _looks_like_sv_token(tok: str) -> bool:
    """Heuristically check if a string looks like our SV token."""
    if not tok or "_" not in tok:
        return False
    s = tok.strip()
    # word form: 25272547_del_59kb (or dup/ins/inv/cnv)
    if re.match(r"^\d+_(del|dup|ins|inv|cnv)_", s, flags=re.I):
        return True
    # sequence form: 126690214_<REF>_<ALT> (long REF preferred)
    parts = s.split("_")
    if len(parts) >= 3:
        p0, p1, p2 = parts[0], parts[1], parts[2]
        if p0.isdigit() and (set(p1) <= set("ACGTNacgtn") and len(p1) >= 20) and (set(p2) <= set("ACGTNacgtn")):
            return True
    return False


def load_db_defs(db_tsv: Path,
                 chrom_col: str | None = None,
                 svtoken_col: str | None = None,
                 ) -> list[SvDef]:
    """Auto-detect and load SV definitions from a DB file.

    Tries hard to cope with variant headers and delimiters. If `chrom_col` or
    `svtoken_col` are provided, they are honored (case-insensitive). Otherwise:
      - Chromosome column is chosen from {'chrom','chr'} (any case).
      - SV tokens are collected from *all* columns that look like SV tokens.

    Args:
        db_tsv (Path): Path to DB TSV/CSV.
        chrom_col (str | None): Column name for chromosome (case-insensitive).
        svtoken_col (str | None): Column with SV token(s); if None, scan all.

    Returns:
        list[SvDef]: Parsed structural variant definitions.
    """
    delim = _sniff_delimiter(db_tsv)
    defs: list[SvDef] = []
    token_hits = 0
    row_count = 0

    with open(db_tsv, "rt", encoding="utf-8-sig", newline="") as fh:
        reader = csv.reader(fh, delimiter=delim)
        header = next(reader, None)
        if header is None:
            return defs

        ci = _ci_lookup(header)

        # Resolve chrom col
        if chrom_col:
            chrom_key = ci.get(chrom_col.lower())
        else:
            chrom_key = ci.get("chrom") or ci.get("chr")

        if chrom_key is None:
            # Try GRCh38 scaffold columns (very repo-specific fallback)
            chrom_key = ci.get("grch38") or ci.get("grch37")
            # If still none, we cannot proceed; tokens don’t have chrom embedded.
            if chrom_key is None:
                return defs

        if svtoken_col:
            token_key = ci.get(svtoken_col.lower())
        else:
            token_key = None  # trigger scan-all

        # Determine allele id column (nice to have)
        allele_key = ci.get("allele") or ci.get("id") or ci.get("genotype") or ci.get("name")

        # Scan rows
        for parts in reader:
            row_count += 1
            row = dict(zip(header, parts))

            chrom = (row.get(chrom_key) or "").strip()
            chrom = chrom.removeprefix("chr").removeprefix("CHR")
            if not chrom:
                continue

            candidates: list[str] = []
            if token_key:
                # Single specified token column
                candidates = [row.get(token_key, "")]
            else:
                # Scan all columns for likely tokens; split by comma/semicolon if present
                for k, v in row.items():
                    if not v:
                        continue
                    # quick path: if it obviously looks like SV token or contains commas with such tokens
                    if _looks_like_sv_token(v):
                        candidates.append(v)
                    elif any(_looks_like_sv_token(t.strip()) for t in re.split(r"[;,]", v)):
                        candidates.append(v)

            if not candidates:
                continue

            allele_id = (row.get(allele_key) or row.get("Genotype") or "-") if allele_key else (row.get("Genotype") or "-")

            for raw in candidates:
                for tok in re.split(r"[;,]", raw):
                    tok = tok.strip()
                    if not tok:
                        continue
                    if not _looks_like_sv_token(tok):
                        continue
                    svd = parse_db_token(chrom, tok, db_id=str(allele_id))
                    if svd:
                        defs.append(svd)
                        token_hits += 1


    return defs


def match_db_to_vcf(db_tsv: Path, vcf_path: Path, min_size: int = 50) -> list[MatchResult]:
    """High-level helper: match DB TSV definitions to VCF SV events.

    Args:
        db_tsv (Path): Path to DB TSV containing allele SV tokens.
        vcf_path (Path): Path to sample VCF (``.vcf`` or ``.vcf.gz``).
        min_size (int): Minimum event size in bp to consider from the VCF.

    Returns:
        list[MatchResult]: Best matches per DB SV token.
    """
    events = list(VcfSvReader(path=vcf_path, min_size=min_size).events())
    db_defs = load_db_defs(db_tsv=db_tsv)
    matcher = SvMatcher()
    return matcher.match(db_defs, events)



@runtime_checkable
class VariantSource(Protocol):
    """Protocol for any source that yields structural variant events."""

    def events(self) -> Iterator["SvEvent"]:
        """Yield structural variant events.

        Returns:
            Iterator[SvEvent]: Sequence of parsed structural variant events.
        """
        ...


@dataclass(slots=True, frozen=True)
class ConfidenceInterval:
    """Confidence intervals for SV breakpoints.

    Attributes:
        left (int): Lower CI bound (inclusive) relative to POS/END.
        right (int): Upper CI bound (inclusive) relative to POS/END.
    """
    left: int = 0
    right: int = 0


@dataclass(slots=True, frozen=True)
class SvEvent:
    """Normalized structural variant description.

    Attributes:
        chrom (str): Chromosome name (no 'chr' prefix).
        pos (int): 1-based left-most position.
        end (int): 1-based end position.
        svtype (str): SV type (DEL, DUP, INS, INV, CNV, BND, INDEL, etc.).
        svlen (int): Signed SV length (DEL negative, DUP/INS positive). 0 if unknown.
        alt (str): ALT field from VCF.
        id (str): Variant ID from VCF.
        qual (str): QUAL column from VCF.
        info (dict[str, str]): Parsed INFO field as strings.
        cipos (ConfidenceInterval): CI around POS.
        ciend (ConfidenceInterval): CI around END.
        sample_fmt (str): Raw FORMAT column for single-sample VCFs.
        sample_value (str): Raw sample value column for single-sample VCFs.
    """

    chrom: str
    pos: int
    end: int
    svtype: str
    svlen: int
    alt: str
    id: str
    qual: str
    info: dict[str, str]
    cipos: ConfidenceInterval = field(default_factory=ConfidenceInterval)
    ciend: ConfidenceInterval = field(default_factory=ConfidenceInterval)
    sample_fmt: str = "."
    sample_value: str = "."

    @property
    def size(self) -> int:
        """Absolute event size.

        Returns:
            int: Absolute size from `svlen` if present, else from `end - pos`.
        """
        if self.svlen != 0:
            return abs(self.svlen)
        return abs(self.end - self.pos)


def _open_text_auto(path: Path) -> io.TextIOBase:
    """Open a plain or gzipped text file in text mode.

    Args:
        path (Path): Path to file.

    Returns:
        io.TextIOBase: Opened file handle.
    """
    if str(path).endswith(".gz"):
        return io.TextIOWrapper(gzip.open(path, mode="rb"), encoding="utf-8", newline="")
    return open(path, "rt", encoding="utf-8", newline="")


def _parse_info(info: str) -> dict[str, str]:
    """Parse a VCF INFO field into a dictionary.

    Args:
        info (str): Raw INFO column.

    Returns:
        dict[str, str]: Parsed key-value pairs, with flags set to "True".
    """
    if info == "." or not info:
        return {}
    out: dict[str, str] = {}
    for field in info.split(";"):
        if not field:
            continue
        if "=" in field:
            k, v = field.split("=", 1)
            out[k] = v
        else:
            out[field] = "True"
    return out


def _parse_ci(val: str | None) -> ConfidenceInterval:
    """Parse a CIPOS/CIEND style value into a ConfidenceInterval.

    Args:
        val (str | None): Comma-separated 'a,b' string.

    Returns:
        ConfidenceInterval: Parsed confidence interval.
    """
    if not val:
        return ConfidenceInterval()
    try:
        a, b = val.split(",")
        return ConfidenceInterval(left=int(a), right=int(b))
    except Exception:
        return ConfidenceInterval()


def _is_large_indel(ref: str, alt: str, threshold: int) -> bool:
    """Determine if REF/ALT implies a large indel.

    Args:
        ref (str): REF allele sequence.
        alt (str): ALT allele sequence(s), comma-separated.
        threshold (int): Minimum size threshold in bp.

    Returns:
        bool: True if large indel detected.
    """
    if alt == "." or alt == "*" or alt.startswith("<"):
        return False
    for a in alt.split(","):
        if abs(len(a) - len(ref)) >= threshold or len(ref) >= threshold or len(a) >= threshold:
            return True
    return False


@dataclass(slots=True, frozen=True)
class VcfSvReader:
    """Portable, minimal structural variant reader for VCF.

    Attributes:
        path (Path): Path to VCF/VCF.GZ file.
        min_size (int): Minimum size threshold for emitting events.
    """

    path: Path
    min_size: int = 50

    def events(self) -> Iterator[SvEvent]:
        """Iterate over structural variant events in a VCF.

        Returns:
            Iterator[SvEvent]: Yielded SV events.
        """
        bnd_cache: dict[str, SvEvent] = {}

        with _open_text_auto(self.path) as fh:
            header_cols: list[str] | None = None
            for line in fh:
                if line.startswith("##"):
                    continue
                if line.startswith("#CHROM"):
                    header_cols = line.rstrip("\n").split("\t")
                    continue
                if not line.strip():
                    continue

                if header_cols is None:
                    raise RuntimeError("VCF header not found before records")

                cols = line.rstrip("\n").split("\t")
                chrom, pos_s, id_s, ref, alt, qual, filt, info_s = cols[:8]
                fmt = "."
                sample = "."
                if len(cols) >= 10:
                    fmt = cols[8]
                    sample = cols[9]
                if chrom.startswith("chr"):
                    chrom = chrom[3:]

                info = _parse_info(info_s)
                pos = int(pos_s)
                end = int(info.get("END", pos_s))
                alt_is_symbolic = alt.startswith("<") and alt.endswith(">")

                svtype = info.get("SVTYPE")
                svlen = int(info["SVLEN"]) if "SVLEN" in info and info["SVLEN"].lstrip("-").isdigit() else 0

                if svtype is None and _is_large_indel(ref, alt, self.min_size):
                    first_alt = alt.split(",")[0]
                    delta = len(first_alt) - len(ref)
                    inferred_type = "DEL" if delta < 0 else ("INS" if delta > 0 else "INDEL")
                    svtype = inferred_type
                    svlen = delta
                    end = pos + max(len(ref), 1)

                if svtype is None and alt_is_symbolic:
                    token = alt.strip("<>")
                    svtype = token.split(":")[0].upper()

                if svtype is None:
                    continue

                cipos = _parse_ci(info.get("CIPOS"))
                ciend = _parse_ci(info.get("CIEND"))

                event = SvEvent(
                    chrom=chrom,
                    pos=pos,
                    end=end,
                    svtype=svtype,
                    svlen=svlen,
                    alt=alt,
                    id=id_s,
                    qual=qual,
                    info=info,
                    cipos=cipos,
                    ciend=ciend,
                    sample_fmt=fmt,
                    sample_value=sample,
                )

                if svtype == "BND":
                    mate_id = event.info.get("MATEID") or event.info.get("MATE") or ""
                    if mate_id:
                        if mate_id in bnd_cache:
                            yield bnd_cache.pop(mate_id)
                            yield event
                        else:
                            bnd_cache[event.id] = event
                    else:
                        yield event
                else:
                    if event.size >= self.min_size:
                        yield event



def select_best_per_vcf(matches: Iterable[MatchResult],
                        tie_tol: float = 1e-6
                        ) -> list[MatchResult]:
    """Select the best match per VCF event, breaking ties by Δpos then Δlen.

    Args:
        matches (Iterable[MatchResult]): All DB↔VCF matches.
        tie_tol (float): Scores within this of the minimum are considered ties.

    Returns:
        list[MatchResult]: Filtered matches; ≤1 per VCF event unless perfect tie after Δpos/Δlen.
    """
    by_vcf: dict[tuple[str, int, int, str], list[MatchResult]] = defaultdict(list)
    for m in matches:
        key = (m.vcf.chrom, m.vcf.pos, m.vcf.end, m.vcf.svtype)
        by_vcf[key].append(m)

    filtered: list[MatchResult] = []
    for group in by_vcf.values():
        # Sort primarily by score, then Δpos, then Δlen
        group.sort(key=lambda r: (r.score, r.pos_delta, r.len_delta))
        best_score = group[0].score
        # Keep only matches within score tolerance
        tied = [g for g in group if abs(g.score - best_score) <= tie_tol]

        if len(tied) > 1:
            # Break ties by Δpos
            tied.sort(key=lambda r: (r.pos_delta, r.len_delta))
            best_pos = tied[0].pos_delta
            tied = [t for t in tied if t.pos_delta == best_pos]

        if len(tied) > 1:
            # Break remaining ties by Δlen
            tied.sort(key=lambda r: r.len_delta)
            best_len = tied[0].len_delta
            tied = [t for t in tied if t.len_delta == best_len]

        # If still tied after Δlen → keep all of them
        filtered.extend(tied)

    # Stable-ish ordering: by VCF, then score, then DB id
    filtered.sort(key=lambda r: (r.vcf.chrom, r.vcf.pos, r.vcf.end, r.score, r.db.id))
    return filtered

def find_sv_events(vcf_path: Path, min_size: int = 50) -> list[SvEvent]:
    """Parse all SV events from a VCF.

    Args:
        vcf_path (Path): Path to VCF/VCF.GZ file.
        min_size (int): Minimum size threshold in bp.

    Returns:
        list[SvEvent]: Parsed SV events.
    """
    reader = VcfSvReader(path=vcf_path, min_size=min_size)
    return list(reader.events())

