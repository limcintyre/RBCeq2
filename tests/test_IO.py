#!/usr/bin/env python3
"""
Unit tests for the VCF Validator Module.
"""

import gzip
import os
import tempfile
import unittest
from unittest.mock import patch

# Import functions and classes from the vcf validator module.
# Adjust the import path if the module is located elsewhere.
from rbceq2.IO.validation import VCFValidator, validate_vcf


class TestVCFValidator(unittest.TestCase):
    """Unit tests for the VCFValidator module."""

    def setUp(self) -> None:
        """Initialize list for temporary file paths."""
        self.temp_files: list[str] = []

    def tearDown(self) -> None:
        """Remove temporary files created during tests."""
        for file_path in self.temp_files:
            try:
                os.remove(file_path)
            except Exception:
                pass

    def create_temp_file(self, content: str, suffix: str = ".vcf") -> str:
        """Create a temporary text file with the given content.

        Args:
            content (str): Text content to write to file.
            suffix (str): Suffix for the temporary file.

        Returns:
            str: Path to the temporary file.
        """
        tmp = tempfile.NamedTemporaryFile(
            delete=False, suffix=suffix, mode="w", encoding="utf-8"
        )
        tmp.write(content)
        tmp.close()
        self.temp_files.append(tmp.name)
        return tmp.name

    def create_temp_file_bytes(self, content: bytes, suffix: str = ".vcf") -> str:
        """Create a temporary binary file with the given content.

        Args:
            content (bytes): Binary content to write to file.
            suffix (str): Suffix for the temporary file.

        Returns:
            str: Path to the temporary file.
        """
        tmp = tempfile.NamedTemporaryFile(delete=False, suffix=suffix)
        tmp.write(content)
        tmp.close()
        self.temp_files.append(tmp.name)
        return tmp.name

    def create_temp_gz_file(self, content: str, suffix: str = ".vcf.gz") -> str:
        """Create a temporary gzipped file with the given content.

        Args:
            content (str): Text content to compress and write.
            suffix (str): Suffix for the temporary file.

        Returns:
            str: Path to the temporary gzipped file.
        """
        tmp = tempfile.NamedTemporaryFile(delete=False, suffix=suffix)
        with gzip.open(tmp.name, "wt", encoding="utf-8") as f:
            f.write(content)
        self.temp_files.append(tmp.name)
        return tmp.name

    def test_nonexistent_file(self) -> None:
        """Test validation for a nonexistent file."""
        result = validate_vcf("nonexistent_file.vcf")
        self.assertFalse(result.is_valid)
        self.assertIn("File does not exist.", result.errors)

    # def test_empty_file(self) -> None:
    #     """Test validation for an empty file."""
    #     file_path = self.create_temp_file("")
    #     result = validate_vcf(file_path)
    #     self.assertFalse(result.is_valid)
    #     self.assertIn("File is empty.", result.errors)
    def test_empty_file(self) -> None:
        """Test validation for an empty file."""
        file_path = self.create_temp_file("")
        result = validate_vcf(file_path)
        self.assertFalse(result.is_valid)
        # Expect header-related errors instead of a "File is empty." message.
        self.assertIn("Missing '##fileformat=VCF' header.", result.errors)
        self.assertIn("Missing header line starting with '#CHROM'.", result.errors)

    def test_invalid_utf8(self) -> None:
        """Test validation for a file with invalid UTF-8 encoding."""
        file_path = self.create_temp_file_bytes(b"\xff\xff\xff")
        result = validate_vcf(file_path)
        self.assertFalse(result.is_valid)
        self.assertIn("File is not valid UTF-8 text.", result.errors)

    def test_missing_fileformat_header(self) -> None:
        """Test validation for a file missing the '##fileformat=VCF' header."""
        content = (
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "chr1\t12345\trs123\tA\tG\t.\tPASS\t.\n"
        )
        file_path = self.create_temp_file(content)
        result = validate_vcf(file_path)
        self.assertFalse(result.is_valid)
        self.assertIn("Missing '##fileformat=VCF' header.", result.errors)

    def test_missing_chrom_header(self) -> None:
        """Test validation for a file missing the '#CHROM' header line."""
        content = "##fileformat=VCFv4.2\nchr1\t12345\trs123\tA\tG\t.\tPASS\t.\n"
        file_path = self.create_temp_file(content)
        result = validate_vcf(file_path)
        self.assertFalse(result.is_valid)
        self.assertIn("Missing header line starting with '#CHROM'.", result.errors)

    def test_invalid_header_columns(self) -> None:
        """Test validation for a file with invalid header columns."""
        content = (
            "##fileformat=VCFv4.2\n"
            "#CHR\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "chr1\t12345\trs123\tA\tG\t.\tPASS\t.\n"
        )
        file_path = self.create_temp_file(content)
        result = validate_vcf(file_path)
        self.assertFalse(result.is_valid)
        self.assertTrue(any("expected '#CHROM'" in error for error in result.errors))

    def test_valid_vcf(self) -> None:
        """Test validation for a correctly formatted VCF file."""
        content = (
            "##fileformat=VCFv4.2\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "chr1\t12345\trs123\tA\tG\t.\tPASS\t.\n"
        )
        file_path = self.create_temp_file(content)
        result = validate_vcf(file_path)
        self.assertTrue(result.is_valid)

    def test_malicious_null_bytes(self) -> None:
        """Test detection of null bytes in the file."""
        content = (
            "##fileformat=VCFv4.2\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "chr1\t12345\trs123\tA\tG\t.\tPASS\t.\n"
            "\0\n"
        )
        file_path = self.create_temp_file(content)
        result = validate_vcf(file_path)
        self.assertFalse(result.is_valid)
        self.assertIn("File contains null byte(s).", result.errors)

    def test_malicious_control_chars(self) -> None:
        """Test detection of suspicious control characters in the file."""
        content = (
            "##fileformat=VCFv4.2\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "chr1\t12345\trs123\tA\tG\t.\tPASS\t.\n"
            "\x07\n"  # Bell character.
        )
        file_path = self.create_temp_file(content)
        result = validate_vcf(file_path)
        self.assertFalse(result.is_valid)
        self.assertIn("File contains suspicious control characters.", result.errors)

    def test_malicious_script(self) -> None:
        """Test detection of embedded script tags in the file."""
        content = (
            "##fileformat=VCFv4.2\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "<script>alert('malicious')</script>\n"
        )
        file_path = self.create_temp_file(content)
        result = validate_vcf(file_path)
        self.assertFalse(result.is_valid)
        self.assertIn("File contains embedded script tags.", result.errors)

    def test_gzipped_file(self) -> None:
        """Test validation for a gzipped VCF file."""
        content = (
            "##fileformat=VCFv4.2\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "chr1\t12345\trs123\tA\tG\t.\tPASS\t.\n"
        )
        file_path = self.create_temp_gz_file(content)
        result = validate_vcf(file_path)
        self.assertTrue(result.is_valid)

    def test_file_size_exceeded(self) -> None:
        """Test validation when file size exceeds the allowed limit."""
        content = (
            "##fileformat=VCFv4.2\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "chr1\t12345\trs123\tA\tG\t.\tPASS\t.\n"
        )
        file_path = self.create_temp_file(content)
        # Monkey patch os.path.getsize to simulate an oversized file.
        with patch("os.path.getsize", return_value=VCFValidator.MAX_FILE_SIZE + 1):
            result = validate_vcf(file_path)
            self.assertFalse(result.is_valid)
            self.assertIn("File size exceeds allowed limit.", result.errors)


if __name__ == "__main__":
    unittest.main()
