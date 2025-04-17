#!/usr/bin/env python3

import uuid
from typing import Any

import pandas as pd
from loguru import logger

from rbceq2.core_logic.constants import DB_VERSION, VERSION, AlleleState
from rbceq2.IO.validation import validate_vcf


def configure_logging(args) -> str:
    """
    Configures the logging for the application.

    Args:
        args: Command-line arguments containing debug flag and log file path.
    """
    log_level = "DEBUG" if args.debug else "INFO"
    logger.remove()  # Remove default logger configuration
    logger.add(
        f"{args.out}_log.txt",  # Log file path from arguments
        level=log_level,
        format="{time} {level} {message}",
        rotation="50 MB",  # Rotate the file when it reaches 10 MB
        compression="zip",  # Compress old logs
    )
    UUID = str(uuid.uuid4())
    logger.info("NOT FOR CLINICAL USE")
    logger.info(f"RBCeq2 Version: {VERSION}")
    logger.info(f"RBCeq2 database Version: {DB_VERSION}")
    logger.info(f"Session UUID: {UUID}")

    return UUID


def record_filtered_data(results: tuple[Any]) -> None:
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
    sample, _, numeric_phenos, alphanumeric_phenos, res = results
    for bg_name, bg_data in res.items():
        if bg_data.filtered_out:
            logger.debug(
                f"Sample: {sample} BG Name: {bg_name}\n"
                f"\n#Results:\n"
                f"Genotypes: {'\n'.join(bg_data.genotypes)}\n"
                f"Phenotypes (numeric): {numeric_phenos.get(bg_name, '')}\n"
                f"Phenotypes (alphanumeric): {alphanumeric_phenos.get(bg_name, '')}\n"
                f"\n#Data:\n"
                f"Vars: {bg_data.variant_pool}\n"
                f"Raw: {'\n' + '\n'.join(map(str, bg_data.alleles[AlleleState.RAW]))}\n"
                f"\n#Filters applied:\n"
            )
            no_filters = True
            for k, v in bg_data.filtered_out.items():
                if v:
                    logger.debug(f"\n{k}: {'\n'.join(map(str, v))}\n")
                    no_filters = False
            if no_filters:
                logger.debug('No filters were applied\n')
            logger.debug('\n______\n')


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
