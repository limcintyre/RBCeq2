import unittest

import pandas as pd

from rbceq2.IO.PDF_reports import (
    _format_cell_content,
    _get_data_for_sample,
    _normalize_sample_id,
    _prepare_dataframes,
)


class TestReportGenerator(unittest.TestCase):
    def setUp(self):
        """Set up sample data for tests."""
        self.sample_suffix = "_GRCh38_1_22_v4.2.1_benchmark_filtered"
        self.geno_data = {
            f"HG001{self.sample_suffix}": {"ABO": "A/B", "RH": "D+/-"},
            f"HG002{self.sample_suffix}.vcf": {"ABO": "O/O", "RH": "D-/-"},
            f"HG003{self.sample_suffix}": {"ABO": None, "RH": "D+/+"},  # Missing ABO
        }
        self.alpha_data = {
            f"HG001{self.sample_suffix}": {"ABO": "AB", "RH": "D+ | D weak"},
            f"HG002{self.sample_suffix}.vcf": {"ABO": "O", "RH": "D-"},
            # HG003 missing from alpha
        }
        self.num_data = {
            f"HG001{self.sample_suffix}": {"ABO": "ABO:1,2", "RH": "RH:1,-5"},
            f"HG002{self.sample_suffix}.vcf": {"ABO": "ABO:-1,-2", "RH": "RH:-1"},
            f"HG003{self.sample_suffix}": {"ABO": "", "RH": "RH:1"},  # Empty ABO
        }

        self.df_geno = pd.DataFrame.from_dict(self.geno_data, orient="index")
        self.df_alpha = pd.DataFrame.from_dict(self.alpha_data, orient="index")
        self.df_num = pd.DataFrame.from_dict(self.num_data, orient="index")

    def test_normalize_sample_id(self):
        self.assertEqual(_normalize_sample_id(f"HG001{self.sample_suffix}"), "HG001")
        self.assertEqual(
            _normalize_sample_id(f"HG002{self.sample_suffix}.vcf"), "HG002"
        )
        self.assertEqual(_normalize_sample_id("SomeOtherID"), "SomeOtherID")

    def test_format_cell_content(self):
        self.assertEqual(_format_cell_content("A/B", separator=","), "A/B")  # No comma
        self.assertEqual(
            _format_cell_content("A / B", separator=","), "A / B"
        )  # No comma
        self.assertEqual(_format_cell_content("A,B", separator=","), "A<br/>B")
        self.assertEqual(_format_cell_content(" A , B ", separator=","), "A<br/>B")
        self.assertEqual(
            _format_cell_content("D+ | D weak", separator=" | "), "D+<br/>D weak"
        )
        self.assertEqual(_format_cell_content(None), "N/A")
        self.assertEqual(_format_cell_content(float("nan")), "N/A")
        self.assertEqual(_format_cell_content(""), "N/A")
 
    def test_prepare_dataframes(self):
        dfs, ids, id_map = _prepare_dataframes(self.df_geno, self.df_alpha, self.df_num)

        self.assertIsInstance(dfs, dict)
        self.assertIn("genotype", dfs)
        self.assertIn("alpha", dfs)
        self.assertIn("numeric", dfs)
        self.assertIsInstance(dfs["genotype"], pd.DataFrame)
        self.assertIn("SampleID_Normalized", dfs["genotype"].columns)
        self.assertEqual(
            set(dfs["genotype"]["SampleID_Normalized"]), {"HG001", "HG002", "HG003"}
        )

        self.assertEqual(ids, {"HG001", "HG002", "HG003"})

        self.assertIn("HG001", id_map)
        self.assertTrue(
            id_map["HG001"].startswith("HG001")
        )  # Check it maps to an original ID
        self.assertEqual(len(id_map), 3)

    def test_prepare_dataframes_missing_input(self):
        # Test with one df missing
        dfs, ids, id_map = _prepare_dataframes(self.df_geno, None, self.df_num)
        self.assertIsNone(dfs["alpha"])
        self.assertEqual(ids, {"HG001", "HG002", "HG003"})
        self.assertEqual(len(id_map), 3)  # Should still find all IDs

    def test_get_data_for_sample(self):
        processed_dfs, _, _ = _prepare_dataframes(
            self.df_geno, self.df_alpha, self.df_num
        )
        sample_data, keys = _get_data_for_sample("HG001", processed_dfs)
        self.assertEqual(sample_data["genotype"]["ABO"], "A/B")
        self.assertEqual(sample_data["alpha"]["ABO"], "AB")
        self.assertEqual(sample_data["numeric"]["ABO"], "ABO:1,2")
        self.assertEqual(keys, {"ABO", "RH"})

        sample_data_hg3, keys_hg3 = _get_data_for_sample("HG003", processed_dfs)
        self.assertEqual(sample_data_hg3["genotype"]["RH"], "D+/+")
        self.assertIsNone(sample_data_hg3["alpha"])  # HG003 has no alpha data
        self.assertEqual(sample_data_hg3["numeric"]["ABO"], "")
        self.assertEqual(keys_hg3, {"ABO", "RH"})  # Keys from geno and num

