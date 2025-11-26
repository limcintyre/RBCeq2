#!/usr/bin/env python3

import uuid
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

    Args:
        args: Command-line arguments (typically from argparse.parse_args()).
    """
    UUID = str(uuid.uuid4())
    log_level = "DEBUG" if args.debug else "INFO"
    log_file_path = f"{args.out}_{UUID}_log.txt"

    logger.remove()
    logger.add(
        log_file_path,
        level=log_level,
        format="{time:YYYY-MM-DD HH:mm:ss.SSS} | {level: <8} | {message}",
        rotation="50 MB",
        compression="zip",
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
        return "\n" + "\n".join(
            [" : ".join([collapse_variant(k), transcripts.get(k, '(None)'), v]) for k, v in pool.items()]
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
                    logger.debug(f"\n{k}: {'\n'.join(map(str, v))}\n")
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
