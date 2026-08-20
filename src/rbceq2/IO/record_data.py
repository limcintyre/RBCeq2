#!/usr/bin/env python3

import uuid
from collections import defaultdict
from functools import cache
from typing import Any

import pandas as pd
from loguru import logger
import argparse

from rbceq2.core_logic.constants import DB_VERSION, VERSION, AlleleState, GENOMIC_TO_TRANSCRIPT_GRCh37, GENOMIC_TO_TRANSCRIPT_GRCh38
from rbceq2.IO.validation import validate_vcf

from rbceq2.core_logic.utils import collapse_variant


def configure_logging(args: argparse.Namespace) -> str:
    """
    Configures the logging for the application and logs arguments line by line.

    One run writes exactly one log file, whatever --processes is set to.

    Note:
        main() forks a multiprocessing Pool and every worker inherits this sink, so the
        sink has to be safe to share. enqueue=True gives each process a queue instead of
        the file handle: records are pickled onto the queue and a single writer thread,
        in this process, is the only thing that ever opens the file. loguru registers
        os.register_at_fork hooks that hold its locks across the fork, so this is safe
        with the default fork start method on Linux.

        Without it every worker wrote to - and rotated - the same path. One
        --processes 12 run of ALL_just_genes.vcf.gz produced 13 log artefacts (a dozen
        truncated .zip fragments plus the .txt), two FileNotFoundError tracebacks on the
        console out of loguru's _terminate_file when two workers rotated at once, and
        only 1256 of the 2318 per sample debug traces survived anywhere on disk. Keep
        enqueue=True for as long as find_hits runs in a Pool.

        Rotation is deliberately not set. A run is one logical log and --debug makes it
        the primary debugging tool, so it stays one greppable .txt rather than 50 MB
        slices that have to be reassembled before anything can be counted. Compression
        goes with it: loguru only compresses at exit when no rotation is configured (see
        _file_sink._terminate_file), so leaving compression="zip" here would start
        zipping the sub 50 MB logs that are plain text today.

    Args:
        args: Command-line arguments (typically from argparse.parse_args()).

    Returns:
        str: The session UUID, also stamped into the header of each output TSV.
    """
    UUID = str(uuid.uuid4())
    log_level = "DEBUG" if args.debug else "INFO"
    log_file_path = f"{args.out}_{UUID}_log.txt"

    logger.remove()
    logger.add(
        log_file_path,
        level=log_level,
        format="{time:YYYY-MM-DD HH:mm:ss.SSS} | {level: <8} | {message}",
        enqueue=True,  # multiprocessing safe - see Note above before changing
    )

    logger.info("=" * 20 + " SESSION START " + "=" * 20)
    logger.info("NOT FOR CLINICAL USE")
    logger.info(f"RBCeq2 Version: {VERSION}")
    logger.info(f"RBCeq2 database Version: {DB_VERSION}")
    logger.info(f"Session UUID: {UUID}")

    logger.info("Command-line arguments provided:")
    args_dict = vars(args)
    if not args_dict:
        logger.info("  (No arguments captured)")
    else:
        max_key_len = max(len(key) for key in args_dict.keys())
        for key, value in args_dict.items():
            logger.info(f"  {key:<{max_key_len}} : {value}")

    logger.info("=" * 20 + " LOGGING STARTED " + "=" * 20)

    return UUID


@cache
def _reference_coding_by_position(ref_genome: str) -> dict[str, str]:
    """Coding notation for the reference base at each position, where it is certain.

    Built once per genome build and cached, because the map is a module level constant
    and the debug trace asks for it three times per blood group per sample.

    A position earns an entry only if every alternate mapped there agrees about what
    the reference is. 1:3775836 has 'c.152T>G' and 'c.152T>A', which both say c.152 is
    T, so it gets 'c.152T'. 4:144120554 has 'c.72G>T' and 'c.71_72delinsGT'; the indel
    names no single reference base so it is skipped, and the SNV settles it as 'c.72G'.
    9:133257521 has only 'c.261del', so nothing is claimed and the position is absent.

    An '=' form is already a statement about the reference - 1:25420739 maps to
    'c.48C=' - so the '=' is dropped rather than the string being rebuilt, which keeps
    the column reading the same way throughout.

    Args:
        ref_genome (str): 'GRCh37' or 'GRCh38'.

    Returns:
        dict[str, str]: 'chrom:pos' -> coding notation for the reference, e.g.
            {'18:45739554': 'c.838G'}. Positions the map cannot settle are absent.
    """
    transcripts = (
        GENOMIC_TO_TRANSCRIPT_GRCh37
        if ref_genome == "GRCh37"
        else GENOMIC_TO_TRANSCRIPT_GRCh38
    )
    forms_by_position: dict[str, set[str]] = defaultdict(set)
    for token, coding in transcripts.items():
        if token.endswith("_ref"):
            continue
        position = token.split("_")[0]
        if coding.endswith("="):
            forms_by_position[position].add(coding[:-1])
        elif ">" in coding:
            forms_by_position[position].add(coding.split(">")[0])

    return {
        position: next(iter(forms))
        for position, forms in forms_by_position.items()
        if len(forms) == 1
    }


def reference_coding(variant: str, ref_genome: str) -> str | None:
    """Coding notation for a '_ref' token, taken from the alternate beside it.

    A '_ref' token asserts the reference at a locus, so it has no entry of its own in
    the genomic to transcript map - only the alternates at that position do. That left
    a lane variant rendering as '(None)' in the debug trace, directly under the
    alternate that already names it:

        18:45739554_G_A : c.838G>A : Heterozygous
        18:45739554_ref : (None) : Heterozygous

    The reference base is read out of the alternate's *coding* string rather than out
    of the genomic token, because the two disagree on the minus strand: 1:25420739_G_G
    maps to 'c.48C=', G in the genome and C in the transcript. Deriving from the
    genomic REF would print the wrong base there.

    Only '_ref' tokens qualify. An unmapped alternate must not borrow from a sibling at
    the same position, because the sibling describes a different change - 1:3775836_T_C
    is 'c.152T>C', not the 'c.152T>G' mapped next to it.

    Args:
        variant (str): A variant token, e.g. '18:45739554_ref'.
        ref_genome (str): 'GRCh37' or 'GRCh38'.

    Returns:
        str | None: 'c.838G' for the example above. None when the token is not a '_ref'
            token, when nothing is mapped at that position, or when the alternates
            there disagree - all of which the caller renders as before. Guessing would
            put a coding string in the trace that the database never said.
    """
    if not variant.endswith("_ref"):
        return None
    return _reference_coding_by_position(ref_genome).get(variant.split("_")[0])


def record_filtered_data(results: tuple[Any], ref: str) -> None:
    """Record filtered data by logging debug information for each blood group.

    This function unpacks the results tuple into sample identifier, numeric and
    alphanumeric phenotypes, and a mapping of blood group names to BloodGroup
    objects. For each blood group with filtered out data, it logs details including
    genotypes, numeric and alphanumeric phenotypes, variant pool, raw allele data,
    and the filters applied.

    Args:
        results (tuple[Any]): A tuple containing the following elements:
            - sample: The sample identifier.
            - _: An unused placeholder.
            - numeric_phenos: A dict mapping blood group names to numeric phenotypes.
            - alphanumeric_phenos: A dict mapping blood group names to
              alphanumeric phenotypes.
            - res: A dict mapping blood group names to BloodGroup objects.
    """

    def format_vars(pool):
        transcripts = GENOMIC_TO_TRANSCRIPT_GRCh37 if ref == 'GRCh37' else GENOMIC_TO_TRANSCRIPT_GRCh38

        def coding(variant):
            curated = transcripts.get(variant)
            if curated is not None:
                return curated
            return reference_coding(variant, ref) or '(None)'

        return "\n" + "\n".join(
            [" : ".join([collapse_variant(k), coding(k), v]) for k, v in sorted(pool.items())]
        )

    sample, genos, numeric_phenos, alphanumeric_phenos, res, var_map = results

    logger.debug("\n### Blood group allele info ###:\n")

    for bg_name, bg_data in res.items():
        if bg_data.filtered_out:
            logger.debug(
                f"Sample: {sample} BG Name: {bg_name}\n"
                f"\n#Results:\n"
                f"Genotypes count: {len(genos.get(bg_name, '').split(','))}\n"
                f"Genotypes: {'\n'.join(genos.get(bg_name, '').split(','))}\n"
                f"Phenotypes (numeric): {'\n'.join(numeric_phenos.get(bg_name, '').split(' | '))}\n"
                f"Phenotypes (alphanumeric): {'\n'.join(alphanumeric_phenos.get(bg_name, '').split(' | '))}\n"
                f"\n#Data:\n"
                f"Vars: {format_vars(bg_data.variant_pool)}\n"
                f"Vars_phase: {format_vars(bg_data.variant_pool_phase)}\n"
                f"Vars_phase_set: {format_vars(bg_data.variant_pool_phase_set)}\n"
                f"Raw: {'\n' + '\n'.join(map(str, bg_data.alleles[AlleleState.RAW]))}\n"
            )

            for variant_db, variant_vcf in var_map.items():
                if variant_db in bg_data.variant_pool:
                    logger.debug(
                        f"BIG VARIANT fuzzy matching map:\n"
                        f"\tDatabase_variant: {collapse_variant(variant_db)}\n"
                        f"\tVCF_variant: {collapse_variant(variant_vcf)}\n"
                    )

            logger.debug("### Filters applied ###:\n")
            no_filters = True
            for k, v in bg_data.filtered_out.items():
                if v:
                    logger.debug(f"\n{k}: {'\n'.join(sorted(map(str, v)))}\n")
                    no_filters = False
            if no_filters:
                logger.debug("No filters were applied\n")
            logger.debug("\n__________________________________________\n")


def check_VCF(VCF_file):
    return validate_vcf(VCF_file), VCF_file


def log_validation(result, VCF_file):
    """Log the validation result for a VCF file.

    Args:
        result (Any): An object representing the validation result. It must have
            attributes 'is_valid' (bool) and 'errors' (iterable of str).
        VCF_file (str): The path or identifier of the VCF file being validated.

    Returns:
        None
    """
    if result.is_valid:
        logger.info(f"VCF file {VCF_file} passed all checks. Proceed with analysis.")
    else:
        logger.error(f"VCF file {VCF_file} failed validation:")
        for error in result.errors:
            logger.warning(f" - {error}")


def save_df(df: pd.DataFrame, name: str, UUID: str) -> None:
    """Sorts the columns of a DataFrame in alphabetical order then writes

    Args:
        df (pd.DataFrame): Data to reorder.

    Returns:
        A DataFrame with columns sorted alphabetically.
    """
    df = df.reindex(sorted(df.columns), axis=1)
    df.index.name = f"UUID: {UUID}"
    df.to_csv(name, sep="\t")


def report_no_call_summary(
    no_calls: dict[tuple[str, str], set[str]], samples_processed: int
) -> None:
    """Warn, once per run, about database loci the caller did not call.

    A no-call means the alleles that need that locus could not be confirmed, so they were
    excluded (see remove_alleles_with_no_call_variants) and the call fell back to the
    reference allele. That is the same contract every other missing-data path in RBCeq2
    follows - a hom ref row is dropped, a non-PASS variant is dropped, a low depth variant
    is dropped, and in each case the result reverts to reference. Applying it silently is
    the problem: a probe that fails in most samples produces confident wildtype calls for
    all of them with nothing in the normal log to say so.

    Worth reading carefully at an antithetical locus, where 'reverted to reference' is a
    positive claim about an antigen rather than an absence - GYPB*04 asserts s+, it does
    not mean 'no S/s result'.

    Real example, from a 1000 Genomes microarray call set of 2318 samples, where the FY
    GATA probe had failed almost completely:

        13 database loci were not called in at least one sample ...
          FY 1:159174683: no call in 2141/2318 samples (92.4%) <- failed in the majority
          HPA3 17:42453065: no call in 413/2318 samples (17.8%)
          VEL 1:3691528: no call in 252/2318 samples (10.9%)

    Args:
        no_calls (dict[tuple[str, str], set[str]]): Maps (blood group, variant) to the set
            of samples in which that variant was Zygosity.NO_DATA.
        samples_processed (int): Total samples in the run, used as the denominator.

    Returns:
        None. Emits loguru warnings; emits nothing at all if there were no no-calls.
    """
    if not no_calls or not samples_processed:
        return

    logger.warning(
        f"{len(no_calls)} database loci were not called in at least one sample. The "
        f"alleles that need them were excluded (no_call_at_defining_variant) and the "
        f"call reverted to reference, so missing or low quality data is reported as "
        f"wildtype. Check these before trusting the affected blood groups:"
    )
    ranked = sorted(no_calls.items(), key=lambda kv: (-len(kv[1]), kv[0]))
    for (bg_name, variant), samples in ranked:
        pct = 100 * len(samples) / samples_processed
        note = " <- failed in the majority of samples" if pct > 50 else ""
        logger.warning(
            f"  {bg_name} {variant}: no call in {len(samples)}/{samples_processed} "
            f"samples ({pct:.1f}%){note}"
        )


def stamps(start: pd.Timestamp) -> str:
    delta = pd.Timestamp.now() - start
    total_seconds = delta.total_seconds()

    # Calculate minutes and remaining seconds
    minutes = int(total_seconds // 60)  # Get whole minutes
    remaining_seconds = total_seconds % 60  # Get the remainder seconds

    # Format the output string conditionally (optional, but nice)
    if minutes > 0:
        time_str = f"{minutes} minutes and {remaining_seconds:.2f} seconds"
    else:
        time_str = f"{remaining_seconds:.2f} seconds"  # Or just total_seconds:.2f

    return time_str
