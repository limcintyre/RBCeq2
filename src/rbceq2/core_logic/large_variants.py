from __future__ import annotations


from collections import defaultdict
from dataclasses import dataclass, field
from typing import Iterator, Protocol, runtime_checkable

import pandas as pd
import re
from typing import Iterable


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
        return SvDef(
            chrom=chrom, pos=pos, svtype=svtype, length=abs(delta), raw=raw, id=db_id
        )

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
    variant: str

    def __repr__(self) -> str:
        return (
            f"MatchResult(score={self.score}, pos_delta={self.pos_delta}, "
            f"len_delta={self.len_delta}, db_id='{self.db.id}', "
            f"db_raw='{self.db.raw}', vcf_variant='{self.vcf.variant}')"
        )


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

    def match(
        self, db_defs: Iterable[SvDef], events: Iterable[SvEvent]
    ) -> list[MatchResult]:
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
                s, pdelta, ld = self.score(db, ev)
                if s != float("inf"):
                    results.append(
                        MatchResult(
                            db=db,
                            vcf=ev,
                            score=s,
                            pos_delta=pdelta,
                            len_delta=ld,
                            variant=ev.variant,
                        )
                    )

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

        POS_FLOOR = 200  # tiny SVs
        POS_FRAC = 0.50  # no-overlap: 50% of size
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

        LEN_FLOOR = 50  # tolerate small caller jitter
        LEN_FRAC = 0.35  # no-overlap: 35% of size
        LEN_FRAC_OVERLAP = 0.50  # overlap: 50% of size
        LEN_CAP = 100_000  # safety cap

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


def _ci_lookup(names: list[str]) -> dict[str, str]:
    """Case-insensitive header map: lower->original.

    Args:
        names (list[str]): Column names.

    Returns:
        dict[str, str]: Mapping from lowercase column names to originals.
    """
    return {n.lower(): n for n in names}


def _looks_like_sv_token(tok: str, min_delta: int = 10) -> bool:
    """Heuristically check if a string looks like our SV token.

    Accepts:
      - word form: 25272547_del_59kb
      - sequence form: <pos>_<REF>_<ALT>, where sequences are DNA and either
        |len(ALT) - len(REF)| >= min_delta or max(len(REF), len(ALT)) >= min_delta.
    """
    if not tok or "_" not in tok:
        return False
    s = tok.strip()

    # word form
    if re.match(r"^\d+_(del|dup|ins|inv|cnv)_", s, flags=re.I):
        return True

    # sequence form
    parts = s.split("_")
    if len(parts) >= 3:
        p0, p1, p2 = parts[0], parts[1], parts[2]
        if p0.isdigit():
            dna = set("ACGTNacgtn")
            if set(p1) <= dna and set(p2) <= dna:
                delta = abs(len(p2) - len(p1))
                if delta >= min_delta or max(len(p1), len(p2)) >= min_delta:
                    return True
    return False


# def _looks_like_sv_token(tok: str) -> bool:
#     """Heuristically check if a string looks like our SV token.

#     Args:
#         tok (str): Candidate string.

#     Returns:
#         bool: True if `tok` looks like an SV token.
#     """
#     if not tok or "_" not in tok:
#         return False
#     s = tok.strip()
#     # word form: 25272547_del_59kb (or dup/ins/inv/cnv)
#     if re.match(r"^\d+_(del|dup|ins|inv|cnv)_", s, flags=re.I):
#         return True
#     # sequence form: 126690214_<REF>_<ALT> (long REF preferred)
#     parts = s.split("_")
#     if len(parts) >= 3:
#         p0, p1, p2 = parts[0], parts[1], parts[2]
#         if (
#             p0.isdigit()
#             and (set(p1) <= set("ACGTNacgtn") and len(p1) >= 20)
#             and (set(p2) <= set("ACGTNacgtn"))
#         ):
#             return True
#     return False


def load_db_defs(
    df: pd.DataFrame,
    chrom_col: str | None = None,
    svtoken_col: str | None = None,
) -> list["SvDef"]:
    """Auto-detect and load SV definitions from a pre-read DB DataFrame.

    Args:
        df (pd.DataFrame): Pre-loaded DB as a DataFrame.
        chrom_col (str | None): Column name for chromosome (case-insensitive).
        svtoken_col (str | None): Column with SV token(s); if None, scan all.

    Returns:
        list[SvDef]: Parsed structural variant definitions.
    """
    defs: list["SvDef"] = []
    token_hits = 0
    ci = _ci_lookup(df.columns.tolist())

    # Resolve chrom col
    if chrom_col:
        chrom_key = ci.get(chrom_col.lower())
    else:
        chrom_key = ci.get("chrom") or ci.get("chr")

    if chrom_key is None:
        # Try GRCh38 scaffold columns (repo-specific fallback)
        chrom_key = ci.get("grch38") or ci.get("grch37")
        if chrom_key is None:
            return defs

    if svtoken_col:
        token_key = ci.get(svtoken_col.lower())
    else:
        token_key = None  # trigger scan-all

    # Determine allele id column (nice to have)
    allele_key = (
        ci.get("allele") or ci.get("id") or ci.get("genotype") or ci.get("name")
    )

    # Iterate rows
    for _, row in df.iterrows():
        chrom = (str(row.get(chrom_key) or "")).strip()
        chrom = chrom.removeprefix("chr").removeprefix("CHR")
        if not chrom:
            continue

        candidates: list[str] = []
        if token_key:
            # Single specified token column
            candidates = [str(row.get(token_key, ""))]
        else:
            # Scan all columns for likely tokens
            for k, v in row.items():
                if not v or not isinstance(v, str):
                    continue
                if _looks_like_sv_token(v):
                    candidates.append(v)
                elif any(_looks_like_sv_token(t.strip()) for t in re.split(r"[;,]", v)):
                    candidates.append(v)

        if not candidates:
            continue

        allele_id = (
            str(row.get(allele_key) or row.get("Genotype") or "-")
            if allele_key
            else (str(row.get("Genotype") or "-"))
        )

        for raw in candidates:
            for tok in re.split(r"[;,]", raw):
                tok = tok.strip()
                if not tok or not _looks_like_sv_token(tok):
                    continue
                svd = parse_db_token(chrom, tok, db_id=str(allele_id))
                if svd:
                    defs.append(svd)
                    token_hits += 1

    return defs


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
    variant: str
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
        if (
            abs(len(a) - len(ref)) >= threshold
            or len(ref) >= threshold
            or len(a) >= threshold
        ):
            return True
    return False


@dataclass(slots=True, frozen=True)
class SnifflesVcfSvReader:
    """Portable, minimal structural variant reader for VCF.

    Attributes:
        df (df): df of to VCF file.
        min_size (int): Minimum size threshold for emitting events.
    """

    df: pd.DataFrame
    min_size: int = 10

    def events(self) -> Iterator[SvEvent]:
        """Iterate over structural variant events in a VCF.
        6
                Returns:
                    Iterator[SvEvent]: Yielded SV events.
        """
        bnd_cache: dict[str, SvEvent] = {}
        for row in self.df.itertuples(index=True, name="Row"):
            assert not row.CHROM.startswith("chr")
            info = _parse_info(row.INFO)
            pos = int(row.POS)
            end = int(info.get("END", row.POS))
            alt_is_symbolic = row.ALT.startswith("<") and row.ALT.endswith(">")

            svtype = info.get("SVTYPE")
            svlen = (
                int(info["SVLEN"])
                if "SVLEN" in info and info["SVLEN"].lstrip("-").isdigit()
                else 0
            )

            if svtype is None and _is_large_indel(row.REF, row.ALT, self.min_size):
                first_alt = row.ALT.split(",")[0]
                delta = len(first_alt) - len(row.REF)
                inferred_type = (
                    "DEL" if delta < 0 else ("INS" if delta > 0 else "INDEL")
                )
                svtype = inferred_type
                svlen = delta
                end = pos + max(len(row.REF), 1)

            if svtype is None and alt_is_symbolic:
                token = row.ALT.strip("<>")
                svtype = token.split(":")[0].upper()

            if svtype is None:
                continue

            cipos = _parse_ci(info.get("CIPOS"))
            ciend = _parse_ci(info.get("CIEND"))

            event = SvEvent(
                chrom=row.CHROM,
                pos=pos,
                end=end,
                svtype=svtype,
                svlen=svlen,
                alt=row.ALT,
                id=row.ID,
                qual=row.QUAL,
                info=info,
                variant=row.variant,
                cipos=cipos,
                ciend=ciend,
                sample_fmt=row.FORMAT,
                sample_value=row.SAMPLE,
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


def select_best_per_vcf(
    matches: Iterable[MatchResult], tie_tol: float = 1e-6, delta_weight: float = 0.5
) -> list[MatchResult]:
    """Select the best match per VCF event, breaking ties by combined Δpos and Δlen.

    Args:
        matches (Iterable[MatchResult]): All DB↔VCF matches.
        tie_tol (float): Scores within this of the minimum are considered ties.
        delta_weight (float): Weight for pos_delta vs len_delta (0.5 = equal weight).
            Higher values (e.g., 0.7) favor positional accuracy over length accuracy.

    Returns:
        list[MatchResult]: Filtered matches; ≤1 per VCF event unless perfect tie.
    """
    by_vcf: dict[tuple[str, int, int, str], list[MatchResult]] = defaultdict(list)
    for m in matches:
        key = (m.vcf.chrom, m.vcf.pos, m.vcf.end, m.vcf.svtype)
        by_vcf[key].append(m)

    filtered: list[MatchResult] = []
    for group in by_vcf.values():
        # Sort primarily by score
        group.sort(key=lambda r: r.score)
        best_score = group[0].score
        # Keep only matches within score tolerance
        tied = [g for g in group if abs(g.score - best_score) <= tie_tol]

        if len(tied) > 1:
            # Normalize deltas within this group for fair comparison
            max_pos = max(t.pos_delta for t in tied)
            max_len = max(t.len_delta for t in tied)

            # Avoid division by zero
            max_pos = max(max_pos, 1)
            max_len = max(max_len, 1)

            # Compute combined delta score
            def combined_delta(r: MatchResult) -> float:
                norm_pos = r.pos_delta / max_pos
                norm_len = r.len_delta / max_len
                return delta_weight * norm_pos + (1 - delta_weight) * norm_len

            tied.sort(key=combined_delta)
            best_combined = combined_delta(tied[0])
            # Keep matches within a small tolerance of best combined score
            tied = [t for t in tied if abs(combined_delta(t) - best_combined) < 1e-9]

        # If still tied, break by DB id for deterministic results
        if len(tied) > 1:
            tied.sort(key=lambda r: r.db.id)
            tied = [tied[0]]

        filtered.extend(tied)

    # Stable ordering: by VCF, then score, then DB id
    filtered.sort(key=lambda r: (r.vcf.chrom, r.vcf.pos, r.vcf.end, r.score, r.db.id))
    return filtered


# def select_best_per_vcf(
#     matches: Iterable[MatchResult], tie_tol: float = 1e-6
# ) -> list[MatchResult]:
#     """Select the best match per VCF event, breaking ties by Δpos then Δlen.

#     Args:
#         matches (Iterable[MatchResult]): All DB↔VCF matches.
#         tie_tol (float): Scores within this of the minimum are considered ties.

#     Returns:
#         list[MatchResult]: Filtered matches; ≤1 per VCF event unless perfect tie after Δpos/Δlen.
#     """
#     by_vcf: dict[tuple[str, int, int, str], list[MatchResult]] = defaultdict(list)
#     for m in matches:
#         key = (m.vcf.chrom, m.vcf.pos, m.vcf.end, m.vcf.svtype)
#         by_vcf[key].append(m)

#     filtered: list[MatchResult] = []
#     for group in by_vcf.values():
#         # Sort primarily by score, then Δpos, then Δlen
#         group.sort(key=lambda r: (r.score, r.pos_delta, r.len_delta))
#         best_score = group[0].score
#         # Keep only matches within score tolerance
#         tied = [g for g in group if abs(g.score - best_score) <= tie_tol]

#         if len(tied) > 1:
#             # Break ties by Δpos
#             tied.sort(key=lambda r: (r.pos_delta, r.len_delta))
#             best_pos = tied[0].pos_delta
#             tied = [t for t in tied if t.pos_delta == best_pos]

#         if len(tied) > 1:
#             # Break remaining ties by Δlen
#             tied.sort(key=lambda r: r.len_delta)
#             best_len = tied[0].len_delta
#             tied = [t for t in tied if t.len_delta == best_len]

#         # If still tied after Δlen → keep all of them
#         filtered.extend(tied)

#     # Stable-ish ordering: by VCF, then score, then DB id
#     filtered.sort(key=lambda r: (r.vcf.chrom, r.vcf.pos, r.vcf.end, r.score, r.db.id))
#     return filtered
