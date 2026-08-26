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
        vcf = VCF([self.sample_df], lane_variants={}, unique_variants={"1:1000"}, sample="test_sample")
        self.assertTrue(hasattr(vcf, "df"))
        # Verify the sample name is stored correctly
        self.assertEqual(vcf.sample, "test_sample")


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
        vcf_obj = VCF([self.df_local], {"9": ["2000"]}, set(), sample="test_sample")
        self.assertIn("variant", vcf_obj.df.columns)

    def test_add_loci(self) -> None:
        """Ensure add_loci method adds the 'loci' column."""
        vcf_obj = VCF([self.test_df], {}, set(), sample="test_sample")
        self.assertIn("loci", vcf_obj.df.columns)

    def test_encode_variants(self) -> None:
        """Check encode_variants method builds the 'variant' column."""
        vcf_obj = VCF([self.test_df], {}, set(), sample="test_sample")
        self.assertIn("variant", vcf_obj.df.columns)
    
    def test_get_sample(self) -> None:
        """Check that the sample attribute is stored correctly."""
        # We pass "test_sample" as the sample name
        vcf_obj = VCF([self.test_df], {}, set(), sample="test_sample")
        
        # Assert that vcf_obj.sample holds the name string, not the genotype data
        self.assertEqual(vcf_obj.sample, "test_sample")

    def test_get_variants(self) -> None:
        """Ensure get_variants constructs a dict with GT-based info."""
        vcf_obj = VCF([self.test_df], {}, set(), sample="test_sample")
        variants = vcf_obj.get_variants()
        self.assertIsInstance(variants, dict)

    def test_remove_home_ref(self) -> None:
        """Check remove_home_ref method removes 0/0 calls."""
        vcf_obj = VCF([self.df_local], {}, set(), sample="test_sample")
        self.assertFalse(any(vcf_obj.df["SAMPLE"].str.startswith("0/0")))

    @staticmethod
    def _df_with_gts(gts: list[str]) -> pd.DataFrame:
        """Return a VCF-like DataFrame with one row per genotype, at distinct positions.

        Args:
            gts (list[str]): GT strings, ie ['0/0', '0|0', '0/1'].

        Returns:
            pd.DataFrame: One row per genotype, all T>C on chr2.
        """
        n = len(gts)
        return pd.DataFrame(
            {
                "CHROM": ["chr2"] * n,
                "POS": [str(2000 + i * 10) for i in range(n)],
                "ID": ["."] * n,
                "REF": ["T"] * n,
                "ALT": ["C"] * n,
                "QUAL": ["."] * n,
                "FILTER": ["."] * n,
                "INFO": ["."] * n,
                "FORMAT": ["GT:AD:GQ:DP:PS"] * n,
                "SAMPLE": [f"{gt}:1,1:30:30:1" for gt in gts],
            }
        )

    def test_remove_home_ref_drops_phased_hom_ref(self) -> None:
        """A phased hom ref ('0|0') must be dropped, exactly as '0/0' is.

        Regression guard for the A4 row of the ploidy state table. Until v2.4.3 only
        '0/0' was matched, so a '0|0' row survived, was encoded as an ALT token and then
        read as Homozygous - asserting the sample carries the variant on both chromosomes
        when the GT says it carries it on neither.
        """
        vcf_obj = VCF(
            [self._df_with_gts(["0/0", "0|0", "0/1"])],
            {},
            set(),
            sample="test_sample",
        )
        surviving_gts = [s.split(":")[0] for s in vcf_obj.df["SAMPLE"]]
        self.assertEqual(surviving_gts, ["0/1"])

    def test_get_variants_excludes_phased_hom_ref(self) -> None:
        """No hom ref call of either separator may reach the variant dict.

        get_variants has its own hom ref skip. It is a safety net behind
        remove_home_ref, and it matched '0/0' only for the same reason.
        """
        vcf_obj = VCF(
            [self._df_with_gts(["0/0", "0|0", "0/1"])],
            {},
            set(),
            sample="test_sample",
        )
        leaked = {
            variant: metrics["GT"]
            for variant, metrics in vcf_obj.variants.items()
            if metrics["GT"] in ("0/0", "0|0")
        }
        self.assertEqual(leaked, {})

    @staticmethod
    def _df_with_two_rows_for_one_variant(
        gts: list[str], filters: list[str]
    ) -> pd.DataFrame:
        """Return a VCF-like DataFrame whose rows are all the same variant.

        The shape a file takes when two callers both report a position and FILTER marks
        which of them is the one that conflicts. One position, one REF/ALT, so every row
        encodes to the same token.

        Args:
            gts (list[str]): One GT per row, in file order.
            filters (list[str]): One FILTER value per row, in file order.

        Returns:
            pd.DataFrame: One row per genotype, all T>C at chr2:2000.
        """
        n = len(gts)
        return pd.DataFrame(
            {
                "CHROM": ["chr2"] * n,
                "POS": ["2000"] * n,
                "ID": ["."] * n,
                "REF": ["T"] * n,
                "ALT": ["C"] * n,
                "QUAL": ["."] * n,
                "FILTER": filters,
                "INFO": ["."] * n,
                "FORMAT": ["GT:AD:GQ:DP:PS"] * n,
                "SAMPLE": [f"{gt}:1,1:30:30:1" for gt in gts],
            }
        )

    def test_get_variants_warns_when_rows_contradict_each_other(self) -> None:
        """Two rows, one token, different genotypes - the discarded one gets named.

        The dict keeps the last row read, so the genotype in use is whichever the file
        ended with. Nothing downstream can see the earlier row, so this is the only
        place the loss can be reported.
        """
        with patch("rbceq2.IO.vcf.logger") as mock_logger:
            vcf_obj = VCF(
                [
                    self._df_with_two_rows_for_one_variant(
                        ["0|1", "0/1"], ["PASS", "TargetedConflict"]
                    )
                ],
                {},
                set(),
                sample="test_sample",
            )
        self.assertEqual(vcf_obj.variants["2:2000_T_C"]["GT"], "0/1")
        warnings = [call.args[0] for call in mock_logger.warning.call_args_list]
        contradiction = [w for w in warnings if "more than one row" in w]
        self.assertEqual(len(contradiction), 1)
        self.assertIn("2:2000_T_C (0|1 then 0/1)", contradiction[0])
        self.assertIn("test_sample", contradiction[0])

    def test_get_variants_silent_when_the_rows_agree(self) -> None:
        """A duplicated row that says the same thing twice changes nothing.

        The overwrite still happens and is still invisible, but there is nothing at
        stake, and a warning that fires on those would teach the reader to skip it.
        """
        with patch("rbceq2.IO.vcf.logger") as mock_logger:
            VCF(
                [
                    self._df_with_two_rows_for_one_variant(
                        ["0/1", "0/1"], ["PASS", "TargetedConflict"]
                    )
                ],
                {},
                set(),
                sample="test_sample",
            )
        warnings = [call.args[0] for call in mock_logger.warning.call_args_list]
        self.assertEqual([w for w in warnings if "more than one row" in w], [])

    def test_remove_home_ref_drops_haploid_zero(self) -> None:
        """Haploid '0' IS hom ref, wherever it is.

        Reversed in v2.4.5. It used to be kept so it would fail loudly downstream, on the
        grounds that whether the region is single copy was undecided. That turned out to
        be the wrong thing to wait for: the caller has said the copies it can see are
        reference, so the ALT token has zero copies under every reading of what the
        haploidy means, and raising about a row that carries no allele helped nobody.

        The ploidy evidence is not lost with the row - _infer_locus_ploidy records it
        first, which is what a mixed coding is detected from later.
        """
        vcf_obj = VCF([self._df_with_gts(["0"])], {}, set(), sample="test_sample")
        self.assertEqual(len(vcf_obj.df), 0)

    def test_haploid_zero_is_still_recorded_as_ploidy_evidence(self) -> None:
        """Dropping the row must not drop the fact that the locus was haploid."""
        vcf_obj = VCF([self._df_with_gts(["0"])], {}, set(), sample="test_sample")
        self.assertEqual(len(vcf_obj.df), 0)
        self.assertTrue(any(vcf_obj.haploid_loci.values()))

    def test_rename_chrom(self) -> None:
        """Check rename_chrom removes the 'chr' prefix."""
        vcf_obj = VCF([self.test_df], {}, set(), sample="test_sample")
        self.assertFalse(vcf_obj.df["CHROM"].str.contains("chr").any())

    def test_set_loci(self) -> None:
        """Check set_loci returns a set of chrom:pos identifiers."""
        vcf_obj = VCF([self.test_df], {}, set(), sample="test_sample")
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


    @patch("rbceq2.IO.vcf.variant_in_intervals", return_value=True)
    def test_read_vcf_header_transformation(self, mock_intervals) -> None:
        """Test that a header with 10 columns is transformed."""
        content = (
            "##fileformat=VCFv4.2\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tEXTRA\n"
            "chr1\t100\t.\tA\tT\t.\tPASS\t.\tGT:AD\tdata1\n"
        )
        file_path = self._create_temp_file(content, suffix=".vcf")
        # CHANGED: Added intervals={}
        df = read_vcf(file_path, intervals={})
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
        self.assertEqual(df["SAMPLE"][0], "data1")
        os.remove(file_path)
    
    @patch("rbceq2.IO.vcf.variant_in_intervals", return_value=True)
    def test_read_vcf_non_gz(self, mock_intervals) -> None:
        """Test reading a non‑gzipped VCF with valid header and data."""
        content = (
            "##fileformat=VCFv4.2\n"
            "##meta-info\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\n"
            "chr1\t100\t.\tA\tT\t.\tPASS\t.\tGT:AD\n"
            "chr2\t200\t.\tG\tC\t.\tPASS\t.\tGT:AD\n"
        )
        file_path = self._create_temp_file(content, suffix=".vcf")
        df = read_vcf(file_path, intervals={})
        self.assertEqual(
            df.columns,
            ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"],
        )
        # UPDATED: Expect 'chr1' and 'chr2' (raw input), not '1' and '2'
        self.assertEqual(df["CHROM"][0], "chr1")
        self.assertEqual(df["CHROM"][1], "chr2")
        os.remove(file_path)

    @patch("rbceq2.IO.vcf.variant_in_intervals", return_value=True)
    def test_read_vcf_gzipped(self, mock_intervals) -> None:
        """Test reading a gzipped VCF file."""
        content = (
            "##fileformat=VCFv4.2\n"
            "##gzipped meta\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\n"
            "chr3\t300\t.\tC\tG\t.\tPASS\t.\tGT:AD\n"
        )
        file_path = self._create_temp_gz_file(content, suffix=".vcf.gz")
        df = read_vcf(file_path, intervals={})
        self.assertEqual(
            df.columns,
            ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"],
        )
        # UPDATED: Expect 'chr3', not '3'
        self.assertEqual(df["CHROM"][0], "chr3")
        os.remove(file_path)
   

    def test_read_vcf_no_header(self) -> None:
        """Test that a VCF with no valid header line raises a ValueError."""
        content = (
            "##fileformat=VCFv4.2\n##meta-info\nchr1\t100\t.\tA\tT\t.\tPASS\t.\tGT:AD\n"
        )
        file_path = self._create_temp_file(content, suffix=".vcf")
        with self.assertRaises(VcfMissingHeaderError) as context:
            # CHANGED: Added intervals={}
            _ = read_vcf(file_path, intervals={})

        name = Path(file_path).name
        message = f"VCF header is missing or invalid in file: '{name}'"
        self.assertEqual(str(context.exception), message)
        os.remove(file_path)



"""
Unit tests for the add_lane_variants method of the VCF class.
"""


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
                lane_variants = {"chr1": ["207331122"], "chr9": ["133257521"]}
                # Initialize VCF with the DataFrame wrapped in a list.
                vcf_obj = VCF(
                    [input_df],
                    lane_variants=lane_variants,
                    unique_variants=set(),
                    sample="test_sample",
                )
                final_df = vcf_obj.df.reset_index(drop=True)

                # Expect two rows: one from the original lane and one newly added.
                self.assertEqual(
                    len(final_df),
                    3,
                    "Expected two rows after adding new lane variants.",
                )


if __name__ == "__main__":
    unittest.main()
