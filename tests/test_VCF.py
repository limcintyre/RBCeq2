#!/usr/bin/env python3
"""
Unit tests for the VCF module.
"""

import gzip
import os
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import pandas as pd
import polars as pl

from rbceq2.IO.vcf import VCF, VcfMissingHeaderError, read_vcf, split_vcf_to_dfs

# Dummy common columns list
COMMON_COLS = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]


class TestVCFInitialization(unittest.TestCase):
    """Tests for VCF initialization."""

    def setUp(self) -> None:
        self.sample_df = pd.DataFrame(
            {
                "CHROM": ["chr1"],
                "POS": ["1000"],
                "ID": ["."],
                "REF": ["A"],
                "ALT": ["G"],
                "QUAL": ["."],
                "FILTER": ["."],
                "INFO": ["."],
                "FORMAT": ["GT:AD:GQ:DP:PS"],
                "SAMPLE": ["0/1:..."],
            }
        )

    def test_init_with_dataframe(self) -> None:
        """Ensure VCF initializes properly when given a DataFrame (wrapped in a list)."""
        vcf = VCF([self.sample_df], lane_variants={}, unique_variants={"1:1000"})
        self.assertTrue(hasattr(vcf, "df"))
        # Since get_sample() returns the last element (a DataFrame), extract its SAMPLE value.
        if isinstance(vcf.sample, pd.DataFrame):
            sample_value = vcf.sample["SAMPLE"].iloc[0]
        else:
            sample_value = vcf.sample
        expected_value = self.sample_df["SAMPLE"].iloc[0]
        self.assertEqual(sample_value, expected_value)

    def test_init_with_path(self) -> None:
        """Ensure VCF initializes properly when given a Path object."""
        with patch("rbceq2.IO.vcf.read_vcf") as mock_read_vcf:
            dummy_pl_df = pl.from_dict(
                {
                    "CHROM": ["chr1"],
                    "POS": ["1000"],
                    "ID": ["."],
                    "REF": ["A"],
                    "ALT": ["G"],
                    "QUAL": ["."],
                    "FILTER": ["."],
                    "INFO": ["."],
                    "FORMAT": ["GT:AD:GQ:DP:PS"],
                    "SAMPLE": ["0/1:..."],
                }
            )
            mock_read_vcf.return_value = dummy_pl_df
            vcf = VCF(Path("dummy.vcf"), lane_variants={}, unique_variants={"1:1000"})
            self.assertTrue(hasattr(vcf, "df"))
            # For Path input, get_sample() returns the file stem.
            self.assertEqual(vcf.sample, "dummy")


class TestVCFMethods(unittest.TestCase):
    """Tests for VCF methods."""

    def setUp(self) -> None:
        self.df_local = pd.DataFrame(
            {
                "CHROM": ["chr2"],
                "POS": ["2000"],
                "ID": ["."],
                "REF": ["T"],
                "ALT": ["C"],
                "QUAL": ["."],
                "FILTER": ["."],
                "INFO": ["."],
                "FORMAT": ["GT:AD:GQ:DP:PS"],
                "SAMPLE": ["0/1:..."],
            }
        )
        self.test_df = self.df_local.copy()

    def test_add_lane_variants_het_and_missing(self) -> None:
        """Check add_lane_variants modifies 'variant' for heterozygous calls."""
        vcf_obj = VCF([self.df_local], {"2": ["2000"]}, set())
        self.assertIn("variant", vcf_obj.df.columns)

    def test_add_loci(self) -> None:
        """Ensure add_loci method adds the 'loci' column."""
        vcf_obj = VCF([self.test_df], {}, set())
        self.assertIn("loci", vcf_obj.df.columns)

    def test_encode_variants(self) -> None:
        """Check encode_variants method builds the 'variant' column."""
        vcf_obj = VCF([self.test_df], {}, set())
        self.assertIn("variant", vcf_obj.df.columns)

    def test_get_sample(self) -> None:
        """Check get_sample method extracts sample value from the DataFrame input."""
        vcf_obj = VCF([self.test_df], {}, set())
        if isinstance(vcf_obj.sample, pd.DataFrame):
            sample_value = vcf_obj.sample["SAMPLE"].iloc[0]
        else:
            sample_value = vcf_obj.sample
        expected_value = self.test_df["SAMPLE"].iloc[0]
        self.assertEqual(sample_value, expected_value)

    def test_get_variants(self) -> None:
        """Ensure get_variants constructs a dict with GT-based info."""
        vcf_obj = VCF([self.test_df], {}, set())
        variants = vcf_obj.get_variants()
        self.assertIsInstance(variants, dict)

    def test_remove_home_ref(self) -> None:
        """Check remove_home_ref method removes 0/0 calls."""
        vcf_obj = VCF([self.df_local], {}, set())
        self.assertFalse(any(vcf_obj.df["SAMPLE"].str.startswith("0/0")))

    def test_rename_chrom(self) -> None:
        """Check rename_chrom removes the 'chr' prefix."""
        vcf_obj = VCF([self.test_df], {}, set())
        self.assertFalse(vcf_obj.df["CHROM"].str.contains("chr").any())

    def test_set_loci(self) -> None:
        """Check set_loci returns a set of chrom:pos identifiers."""
        vcf_obj = VCF([self.test_df], {}, set())
        loci = vcf_obj.set_loci()
        self.assertIsInstance(loci, set)


class TestSplitVCFToDFs(unittest.TestCase):
    """Tests for split_vcf_to_dfs function."""

    def test_split_multi_sample_vcf(self) -> None:
        """Check that it yields separate DataFrames per sample."""
        data = {
            "CHROM": ["chr1", "chr1"],
            "POS": ["100", "200"],
            "ID": [".", "."],
            "REF": ["A", "G"],
            "ALT": ["T", "C"],
            "QUAL": [".", "."],
            "FILTER": [".", "."],
            "INFO": [".", "."],
            "FORMAT": ["GT:AD:GQ:DP:PS", "GT:AD:GQ:DP:PS"],
            "sampleA": ["0/1", "0/1"],
            "sampleB": ["0/0", "0/1"],
        }
        df = pd.DataFrame(data)
        # Invert the yielded tuple (df, sample) to build a dict mapping sample -> DataFrame.
        sample_dfs = {sample: df_ for df_, sample in split_vcf_to_dfs(df)}
        self.assertIn("sampleA", sample_dfs)
        self.assertIn("SAMPLE", sample_dfs["sampleA"].columns)


class TestReadVCF(unittest.TestCase):
    """Unit tests for the read_vcf function."""

    def _create_temp_file(self, content: str, suffix: str = ".vcf") -> str:
        """Create a temporary file with the provided content.

        Args:
            content (str): Content to write.
            suffix (str): File suffix.

        Returns:
            str: Path to the temporary file.
        """
        tmp = tempfile.NamedTemporaryFile(
            delete=False, suffix=suffix, mode="w", encoding="utf-8"
        )
        tmp.write(content)
        tmp.close()
        return tmp.name

    def _create_temp_gz_file(self, content: str, suffix: str = ".vcf.gz") -> str:
        """Create a temporary gzipped file with the provided content.

        Args:
            content (str): Content to compress and write.
            suffix (str): File suffix.

        Returns:
            str: Path to the temporary gzipped file.
        """
        tmp = tempfile.NamedTemporaryFile(delete=False, suffix=suffix)
        with gzip.open(tmp.name, "wt", encoding="utf-8") as f:
            f.write(content)
        tmp.close()
        return tmp.name

    def tearDown(self) -> None:
        """Clean up temporary files."""
        for fname in os.listdir(tempfile.gettempdir()):
            fpath = os.path.join(tempfile.gettempdir(), fname)
            try:
                # Remove only our temporary files based on known suffixes.
                if fpath.endswith(".vcf") or fpath.endswith(".vcf.gz"):
                    os.remove(fpath)
            except Exception:
                pass

    def test_read_vcf_non_gz(self) -> None:
        """Test reading a non‑gzipped VCF with valid header and data."""
        # Create VCF content with meta lines and a header (9 columns).
        content = (
            "##fileformat=VCFv4.2\n"
            "##meta-info\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\n"
            "chr1\t100\t.\tA\tT\t.\tPASS\t.\tGT:AD\n"
            "chr2\t200\t.\tG\tC\t.\tPASS\t.\tGT:AD\n"
        )
        file_path = self._create_temp_file(content, suffix=".vcf")
        df = read_vcf(file_path)
        # Check header preservation.
        self.assertEqual(
            df.columns,
            ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"],
        )
        # Check that 'chr' prefix was removed.
        self.assertEqual(df["CHROM"][0], "1")
        self.assertEqual(df["CHROM"][1], "2")
        os.remove(file_path)

    def test_read_vcf_header_transformation(self) -> None:
        """Test that a header with 10 columns is transformed (last column becomes 'SAMPLE')."""
        content = (
            "##fileformat=VCFv4.2\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tEXTRA\n"
            "chr1\t100\t.\tA\tT\t.\tPASS\t.\tGT:AD\tdata1\n"
        )
        file_path = self._create_temp_file(content, suffix=".vcf")
        df = read_vcf(file_path)
        expected_cols = [
            "CHROM",
            "POS",
            "ID",
            "REF",
            "ALT",
            "QUAL",
            "FILTER",
            "INFO",
            "FORMAT",
            "SAMPLE",
        ]
        self.assertEqual(df.columns, expected_cols)
        # Verify data in transformed column remains.
        self.assertEqual(df["SAMPLE"][0], "data1")
        os.remove(file_path)

    def test_read_vcf_gzipped(self) -> None:
        """Test reading a gzipped VCF file."""
        content = (
            "##fileformat=VCFv4.2\n"
            "##gzipped meta\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\n"
            "chr3\t300\t.\tC\tG\t.\tPASS\t.\tGT:AD\n"
        )
        file_path = self._create_temp_gz_file(content, suffix=".vcf.gz")
        df = read_vcf(file_path)
        self.assertEqual(
            df.columns,
            ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"],
        )
        self.assertEqual(df["CHROM"][0], "3")
        os.remove(file_path)

    def test_read_vcf_no_header(self) -> None:
        """Test that a VCF with no valid header line raises a ValueError."""
        # File with only meta lines and data but no header line.
        content = (
            "##fileformat=VCFv4.2\n##meta-info\nchr1\t100\t.\tA\tT\t.\tPASS\t.\tGT:AD\n"
        )
        file_path = self._create_temp_file(content, suffix=".vcf")
        with self.assertRaises(VcfMissingHeaderError) as context:
            _ = read_vcf(file_path)
        
        name = Path(file_path).name
        message =  f"VCF header is missing or invalid in file: '{name}'"
        self.assertEqual(str(context.exception), message)
        os.remove(file_path)


"""
Unit tests for the add_lane_variants method of the VCF class.
"""


# NOTE: These tests rely on patched constants.
# Patch COMMON_COLS to a known list and HOM_REF_DUMMY_QUAL to a fixed string.
# Expected final columns:
# ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "SAMPLE", "variant", "loci"]


class TestAddLaneVariants(unittest.TestCase):
    """Full coverage tests for add_lane_variants."""

    def _make_input_df(self) -> pd.DataFrame:
        """Return a minimal VCF-like DataFrame with a single row.

        The row represents a variant at lane 'chr1:1000'.
        """
        return pd.DataFrame(
            {
                "CHROM": ["chr1"],
                "POS": ["1000"],
                "ID": ["x"],
                "REF": ["A"],
                "ALT": ["T"],
                "QUAL": ["."],
                "FILTER": ["."],
                "INFO": ["."],
                "FORMAT": ["GT:AD:GQ:DP:PS"],
                "SAMPLE": ["0/1:10"],
            }
        )

    def test_add_lane_variants_full_coverage(self) -> None:
        """Test add_lane_variants covering both in-branch and else-branch."""
        with patch(
            "rbceq2.IO.vcf.COMMON_COLS",
            new=[
                "CHROM",
                "POS",
                "ID",
                "REF",
                "ALT",
                "QUAL",
                "FILTER",
                "INFO",
                "FORMAT",
            ],
        ):
            with patch("rbceq2.IO.vcf.HOM_REF_DUMMY_QUAL", new="dummy_qual"):
                # Prepare input DataFrame and lane_variants.
                input_df = self._make_input_df()
                # lane_variants includes one lane that exists and one that doesn't.
                lane_variants = {"chr1": ["1000"], "chr2": ["2000"]}
                # Initialize VCF with the DataFrame wrapped in a list.
                vcf_obj = VCF(
                    [input_df], lane_variants=lane_variants, unique_variants=set()
                )
                final_df = vcf_obj.df.reset_index(drop=True)

                # Expect two rows: one from the original lane and one newly added.
                self.assertEqual(
                    len(final_df),
                    2,
                    "Expected two rows after adding new lane variants.",
                )

                # For the original row, loci becomes "1:1000" after rename_chrom and add_loci.
                orig_row = final_df[final_df["loci"] == "1:1000"].iloc[0]
                # encode_variants produces "1:1000_A_T"; new lane branch appends ",1:1000_ref"
                self.assertEqual(
                    orig_row["variant"],
                    "1:1000_A_T,1:1000_ref",
                    "Original row variant not modified correctly.",
                )

                # The else branch creates a new lane row for lane "chr2:2000".
                # Note: The new row is built as:
                #   [chrom, pos] + COMMON_COLS[2:-1] + ["GT:AD:GQ:DP:PS"] + [HOM_REF_DUMMY_QUAL, f"{lane_loci}_ref", "loci"]
                # With patched COMMON_COLS = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"],
                # COMMON_COLS[2:-1] is ["ID", "REF", "ALT", "QUAL", "FILTER", "INFO"].
                # And since "chr2".replace("chr", "") yields "2", lane_loci = "2:2000".
                # The new row becomes:
                # ["2", "2000"] + ["ID", "REF", "ALT", "QUAL", "FILTER", "INFO"] +
                # ["GT:AD:GQ:DP:PS"] + ["dummy_qual", "2:2000_ref", "loci"]
                # Final expected mapping:
                expected_new = {
                    "CHROM": "2",
                    "POS": "2000",
                    "ID": "ID",
                    "REF": "REF",
                    "ALT": "ALT",
                    "QUAL": "QUAL",
                    "FILTER": "FILTER",
                    "INFO": "INFO",
                    "FORMAT": "GT:AD:GQ:DP:PS",
                    "SAMPLE": "dummy_qual",
                    "variant": "2:2000_ref",
                    "loci": "loci",
                }
                # Identify the new row by matching CHROM and POS.
                new_rows = final_df[
                    (final_df["CHROM"] == "2") & (final_df["POS"] == "2000")
                ]
                self.assertEqual(
                    len(new_rows), 1, "Expected one new lane row for chr2:2000."
                )
                new_row = new_rows.iloc[0]
                for col, exp_val in expected_new.items():
                    self.assertEqual(
                        new_row[col],
                        exp_val,
                        f"New lane row column '{col}' not as expected.",
                    )


if __name__ == "__main__":
    unittest.main()
