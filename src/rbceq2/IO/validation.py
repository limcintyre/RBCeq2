#!/usr/bin/env python3
"""
Optimized VCF Validator Module with Faster Control Character Check
------------------------------
This module validates VCF files for correct formatting and scans for malicious content in a single pass.
It replaces the regex-based control character check with a set-based lookup, which is faster.
"""

import os
import gzip
from dataclasses import dataclass

# Precompute a set of disallowed control characters (allowed: tab, newline, carriage return)
SUSPICIOUS_CHARS = {chr(c) for c in list(range(0, 9)) + list(range(11, 13)) + list(range(14, 32))}


@dataclass
class VCFValidationResult:
    """
    A data class representing the result of VCF validation.
    ----------
    Attributes
    ----------
    is_valid (bool):
        Indicates if the validation passed.
    errors (list[str]):
        List of error messages encountered during validation.
    """
    is_valid: bool
    errors: list[str]


def open_vcf_file(file_path: str, mode: str = "rt"):
    """
    Open a VCF file that may be plain text or gzip compressed.
    
    Args:
        file_path (str): Path to the VCF file.
        mode (str, optional): File open mode. Defaults to "rt".
    
    Returns:
        A file object.
    """
    if str(file_path).endswith(".gz"):
        return gzip.open(file_path, mode, encoding="utf-8")
    return open(file_path, mode, encoding="utf-8")


class VCFValidator:
    """
    A class for validating VCF files in a single pass.
    """
    MAX_FILE_SIZE = 2 * 1024 * 1024 * 1024  # 2GB file size limit

    @staticmethod
    def validate_file(file_path: str) -> VCFValidationResult:
        """
        Validate the VCF file for correct format and malicious content in a single pass.
        
        Args:
            file_path (str): Path to the VCF file.
        
        Returns:
            VCFValidationResult: Result of validation.
        """
        errors: list[str] = []

        if not os.path.exists(file_path):
            return VCFValidationResult(False, ["File does not exist."])

        file_size: int = os.path.getsize(file_path)
        if file_size > VCFValidator.MAX_FILE_SIZE:
            return VCFValidationResult(False, ["File size exceeds allowed limit."])

        found_fileformat = False
        header_line: str | None = None

        try:
            with open_vcf_file(file_path, "rt") as f:
                for line in f:
                    # Check for malicious content.
                    if "\0" in line:
                        errors.append("File contains null byte(s).")
                    # Replace regex search with set intersection for control characters.
                    if set(line) & SUSPICIOUS_CHARS:
                        errors.append("File contains suspicious control characters.")
                    if "<script" in line.lower():
                        errors.append("File contains embedded script tags.")

                    # Check for meta-information header.
                    if not found_fileformat and line.startswith("##fileformat=VCF"):
                        found_fileformat = True
                    # Find header line: first line starting with '#' but not '##'.
                    if header_line is None and line.startswith("#") and not line.startswith("##"):
                        header_line = line.strip()
        except UnicodeDecodeError:
            return VCFValidationResult(False, ["File is not valid UTF-8 text."])
        except Exception as e:
            return VCFValidationResult(False, [f"Error reading file: {str(e)}"])

        if not found_fileformat:
            errors.append("Missing '##fileformat=VCF' header.")

        if header_line is None:
            errors.append("Missing header line starting with '#CHROM'.")
        else:
            columns = header_line.split("\t")
            if len(columns) < 8:
                errors.append("Header line has insufficient columns.")
            else:
                required = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]
                for i, col in enumerate(required):
                    if columns[i] != col:
                        errors.append(f"Header column {i + 1} expected '{col}', found '{columns[i]}'.")
                        break

        # Remove duplicate errors.
        unique_errors = list(dict.fromkeys(errors))
        return VCFValidationResult(len(unique_errors) == 0, unique_errors)


def validate_vcf(file_path: str) -> VCFValidationResult:
    """
    Validate a VCF file for proper format and absence of malicious content.
    
    Args:
        file_path (str): Path to the VCF file.
    
    Returns:
        VCFValidationResult: Combined result of format and security checks.
    """
    return VCFValidator.validate_file(file_path)