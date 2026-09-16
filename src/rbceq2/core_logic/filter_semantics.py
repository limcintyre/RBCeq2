"""What a VCF FILTER value is a statement about.

FILTER is one value per row, and a row is a site, so the column is a grab-bag: some values
say the caller doubts this sample's genotype, some report a statistic computed across a whole
cohort, and some are neither - they mark which of several rows for a marker is the recommended
one, or complain about a breakpoint coordinate this tool deliberately rounds away.

Only the first kind is grounds to drop an allele. Treating all three the same is what made a
sample's allele disappear because of other samples' data on a jointly called cohort, and what
would revert a sample to the reference allele on an array where PASS/FAIL is probeset
selection rather than call quality.

The judgement is made once, by a person, when a value is added to resources/filter_values.tsv,
and the Evidence column records the caller's own description of the value so the
classification can be checked. Nothing is decided at run time - this module is a lookup, and
the tool stays deterministic.

An unrecognised value excludes, as before. A new caller's vocabulary is not assumed benign;
the value is named in a warning so it can be classified.

One limitation worth knowing. A few values are site-scoped on a jointly called file and
effectively sample-scoped on a single-sample one, because a 'site' with one sample in it is
that sample. Only values whose statistic is unambiguously computed across samples belong in
site_statistic - genotyping rate, for instance, is arithmetic over the cohort and cannot be
about any one member of it. Anything less clear-cut is better left out of the table, where it
keeps the old behaviour.
"""

import csv
import importlib.resources
from functools import cache
from io import StringIO

from loguru import logger

PASS = "PASS"
NO_FILTER_APPLIED = frozenset({".", ""})
FILTER_VALUE_SEPARATOR = ";"


class FilterMeaning:
    """What kind of statement a FILTER value makes."""

    SAMPLE_CALL = "sample_call"
    SITE_STATISTIC = "site_statistic"
    NOT_CORRECTNESS = "not_correctness"


EXCLUDING = frozenset({FilterMeaning.SAMPLE_CALL})


def load_filter_values() -> str:
    """Load the filter_values.tsv file from package resources.

    Returns:
        str: The raw contents of the table.

    Raises:
        Exception: If the resource cannot be read, matching load_db - a missing table is a
            broken install, not something to paper over with an empty default that would
            silently change every filtering decision.
    """
    try:
        resource_path = importlib.resources.files("rbceq2").joinpath(
            "resources", "filter_values.tsv"
        )
        logger.debug(f"Attempting to load FILTER values from: {resource_path}")
        return resource_path.read_text(encoding="utf-8")
    except Exception as e:
        logger.error(
            "Failed to load resource 'resources/filter_values.tsv' from package "
            f"'rbceq2': {e}"
        )
        raise


@cache
def filter_value_meanings() -> dict[str, str]:
    """Map each known FILTER value to what it is a statement about.

    Cached because the table is static for the life of a run and this is called once per
    defining variant per allele per sample.

    Returns:
        dict[str, str]: FILTER value -> one of FilterMeaning's members.

    Raises:
        ValueError: If the table carries a classification that is not a FilterMeaning
            member. A typo there would silently change filtering for every sample, so it
            is caught on load rather than read as 'unknown'.
    """
    known = {
        FilterMeaning.SAMPLE_CALL,
        FilterMeaning.SITE_STATISTIC,
        FilterMeaning.NOT_CORRECTNESS,
    }
    meanings: dict[str, str] = {}
    for row in csv.DictReader(StringIO(load_filter_values()), delimiter="\t"):
        value = (row["Value"] or "").strip()
        classification = (row["Classification"] or "").strip()
        if not value:
            continue
        if classification not in known:
            raise ValueError(
                f"filter_values.tsv gives '{value}' the classification "
                f"'{classification}', which is not one of {sorted(known)}"
            )
        meanings[value] = classification

    return meanings


def filter_excludes_allele(filter_field: str) -> tuple[bool, tuple[str, ...]]:
    """Decide whether a row's FILTER field is grounds to drop an allele.

    PASS and an absent filter ('.' or empty, which the VCF specification uses for 'no
    filtering applied') never exclude. Otherwise the field is split on ';' - a row may carry
    several values, and 28 distinct combinations appear at once in a jointly called cohort -
    and the allele is dropped if any one of them says this sample's call is doubtful, or if
    any one of them is unrecognised.

    Args:
        filter_field (str): The FILTER column of the row, verbatim.

    Returns:
        tuple[bool, tuple[str, ...]]: Whether to drop the allele, and any values that were
            not in the table. The caller names those in a warning - they are the signal
            that a new input type needs classifying.

    Example:
        >>> filter_excludes_allele("PASS")
        (False, ())
        >>> filter_excludes_allele("LowQual")
        (True, ())
        >>> filter_excludes_allele("LowMLSQ")
        (False, ())
        >>> filter_excludes_allele("cnvLength;LowQual")
        (True, ())
    """
    field = (filter_field or "").strip()
    if field == PASS or field in NO_FILTER_APPLIED:
        return False, ()

    meanings = filter_value_meanings()
    unclassified: list[str] = []
    excludes = False
    for value in field.split(FILTER_VALUE_SEPARATOR):
        value = value.strip()
        if not value or value == PASS:
            continue
        meaning = meanings.get(value)
        if meaning is None:
            unclassified.append(value)
            excludes = True
        elif meaning in EXCLUDING:
            excludes = True

    return excludes, tuple(unclassified)
