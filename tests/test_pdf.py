import unittest
from pathlib import Path
from unittest.mock import patch

import pandas as pd

from rbceq2.IO import PDF_reports as pdf
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
        for sample_id in (
            "HG001", f"HG001{self.sample_suffix}",
            f"HG002{self.sample_suffix}.vcf", "SomeOtherID", " sample ",
        ):
            with self.subTest(sample_id=sample_id):
                self.assertEqual(_normalize_sample_id(sample_id), sample_id)

    def test_invalid_sample_ids_are_rejected(self):
        for sample_id in ("", None, 1, 1.5):
            with self.subTest(sample_id=sample_id):
                with self.assertRaises(ValueError):
                    _normalize_sample_id(sample_id)
                frame = pd.DataFrame({"ABO": ["A/B"]}, index=[sample_id])
                with self.assertRaises(ValueError):
                    _prepare_dataframes(frame, None, None)

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
        pd.testing.assert_frame_equal(dfs["genotype"], self.df_geno)
        pd.testing.assert_frame_equal(dfs["alpha"], self.df_alpha)
        pd.testing.assert_frame_equal(dfs["numeric"], self.df_num)
        self.assertEqual(ids, set(self.geno_data))
        self.assertEqual(id_map, {sample_id: sample_id for sample_id in ids})

    def test_prepare_dataframes_missing_input(self):
        # Test with one df missing
        dfs, ids, id_map = _prepare_dataframes(self.df_geno, None, self.df_num)
        self.assertIsNone(dfs["alpha"])
        self.assertEqual(ids, set(self.geno_data))
        self.assertEqual(len(id_map), 3)  # Should still find all IDs

    def test_get_data_for_sample(self):
        processed_dfs, _, _ = _prepare_dataframes(
            self.df_geno, self.df_alpha, self.df_num
        )
        sample_data, keys = _get_data_for_sample(
            f"HG001{self.sample_suffix}", processed_dfs
        )
        self.assertEqual(sample_data["genotype"]["ABO"], "A/B")
        self.assertEqual(sample_data["alpha"]["ABO"], "AB")
        self.assertEqual(sample_data["numeric"]["ABO"], "ABO:1,2")
        self.assertEqual(keys, {"ABO", "RH"})

        sample_data_hg3, keys_hg3 = _get_data_for_sample(
            f"HG003{self.sample_suffix}", processed_dfs
        )
        self.assertEqual(sample_data_hg3["genotype"]["RH"], "D+/+")
        self.assertIsNone(sample_data_hg3["alpha"])  # HG003 has no alpha data
        self.assertEqual(sample_data_hg3["numeric"]["ABO"], "")
        self.assertEqual(keys_hg3, {"ABO", "RH"})  # Keys from geno and num


class TestPdfIdentity(unittest.TestCase):
    """Keep sample identity intact in rows, displayed text, and output paths."""

    def test_similar_ids_keep_their_own_rows(self):
        ids = ["sample", "sample.vcf", "sample_GRCh38_1_22_v4.2.1_benchmark_filtered"]
        genotype = pd.DataFrame({"ABO": ["A/A", "B/B", "O/O"]}, index=ids)
        alpha = pd.DataFrame({"ABO": ["O", "A"]}, index=[ids[2], ids[0]])
        numeric = pd.DataFrame(
            {"ABO": ["ABO:2", "ABO:1", "ABO:-1,-2"]},
            index=[ids[1], ids[0], ids[2]],
        )
        processed, found_ids, id_map = _prepare_dataframes(genotype, alpha, numeric)
        self.assertEqual(found_ids, set(ids))
        self.assertEqual(id_map, {sample_id: sample_id for sample_id in ids})
        expected = [
            ("A/A", "A", "ABO:1"),
            ("B/B", None, "ABO:2"),
            ("O/O", "O", "ABO:-1,-2"),
        ]
        for sample_id, (geno, pheno_alpha, pheno_num) in zip(ids, expected):
            with self.subTest(sample_id=sample_id):
                data, keys = _get_data_for_sample(sample_id, processed)
                self.assertEqual(data["genotype"]["ABO"], geno)
                if pheno_alpha is None:
                    self.assertIsNone(data["alpha"])
                else:
                    self.assertEqual(data["alpha"]["ABO"], pheno_alpha)
                self.assertEqual(data["numeric"]["ABO"], pheno_num)
                self.assertEqual(keys, {"ABO"})

    def test_duplicate_indices_are_rejected_in_each_input(self):
        duplicated = pd.DataFrame({"ABO": ["A/A", "B/B"]}, index=["same", "same"])
        for position in range(3):
            with self.subTest(position=position):
                frames = [None, None, None]
                frames[position] = duplicated
                with self.assertRaises(ValueError):
                    _prepare_dataframes(*frames)

    def test_input_frames_and_named_column_are_preserved(self):
        frames = [
            pd.DataFrame({"SampleID_Normalized": [value], "ABO": ["A/B"]},
                         index=["sample.vcf"])
            for value in ("genotype value", "alpha value", "numeric value")
        ]
        originals = [frame.copy(deep=True) for frame in frames]
        processed, _, _ = _prepare_dataframes(*frames)
        data, keys = _get_data_for_sample("sample.vcf", processed)
        self.assertEqual(keys, {"SampleID_Normalized", "ABO"})
        for kind, frame, original in zip(
            ("genotype", "alpha", "numeric"), frames, originals
        ):
            with self.subTest(kind=kind):
                pd.testing.assert_frame_equal(frame, original)
                pd.testing.assert_frame_equal(processed[kind], original)
                self.assertEqual(data[kind]["SampleID_Normalized"], original.iloc[0, 0])

    def test_filenames_are_distinct_bounded_and_safe(self):
        ids = [
            "sample", "sample.vcf", "sample/one", "sample?one", "sample_one",
            "Sample", "SAMPLE", "sample__2", "../sample", "a\\b", "é", "É",
            "<sample>&", "a" * 1000, "a" * 999 + "b",
        ]
        names = [pdf._pdf_filename(sample_id) for sample_id in ids]
        self.assertEqual(len(set(name.casefold() for name in names)), len(ids))
        for sample_id, name in zip(ids, names):
            with self.subTest(sample_id=sample_id):
                self.assertEqual(name, pdf._pdf_filename(sample_id))
                self.assertTrue(name.isascii())
                self.assertLessEqual(len(name.encode("ascii")), 255)
                self.assertEqual(Path(name).name, name)
                self.assertNotIn("/", name)
                self.assertNotIn("\\", name)
                self.assertTrue(name.endswith("_BloodGroupReport.pdf"))
        self.assertEqual(
            list(reversed(names)),
            [pdf._pdf_filename(sample_id) for sample_id in reversed(ids)],
        )
        suffix_like_id = names[0].removesuffix("_BloodGroupReport.pdf")
        self.assertNotIn(pdf._pdf_filename(suffix_like_id).casefold(),
                         {name.casefold() for name in names})

    def test_case_insensitive_filename_collision_fails_before_writes(self):
        frame = pd.DataFrame({"ABO": ["A/A", "B/B"]}, index=["A", "B"])
        with patch.object(pdf, "_pdf_filename", side_effect=["same.pdf", "SAME.pdf"]), \
                patch.object(pdf.os, "makedirs") as mkdir, \
                patch.object(pdf, "_generate_pdf_report_for_sample") as render:
            with self.assertRaises(ValueError):
                pdf.generate_all_reports(frame, None, None, Path("unused"), "test")
        mkdir.assert_not_called()
        render.assert_not_called()

    def test_header_displays_literal_sample_id(self):
        for sample_id in ("<b>sample</b>&one", "sample &lt;two&gt;", 'sample"three'):
            with self.subTest(sample_id=sample_id):
                story = []
                pdf._create_report_header(story, sample_id, pdf._setup_styles())
                paragraphs = [item.getPlainText() for item in story
                              if hasattr(item, "getPlainText")]
                self.assertIn(f"Sample ID: {sample_id}", paragraphs)
                self.assertIn("Not for clinical use", paragraphs)


class TestPdfFailures(unittest.TestCase):
    """Surface report failures while still attempting the remaining samples."""

    def setUp(self):
        self.frame = pd.DataFrame(
            {"YT": ["YT*01/YT*01"] * 3}, index=["A", "B", "C"]
        )

    def test_render_failure_reaches_the_batch_caller(self):
        """Do not swallow a write error from the real single-report builder."""
        processed, _, _ = pdf._prepare_dataframes(self.frame, None, None)
        with patch.object(pdf.BaseDocTemplate, "build", side_effect=OSError("write denied")):
            with self.assertRaisesRegex(OSError, "write denied"):
                pdf._generate_pdf_report_for_sample(
                    "A", "A", processed, pdf._setup_styles(), Path("unused"), "test"
                )

    def test_batch_attempts_every_sample_then_reports_failures(self):
        """Keep successful samples and name every failed sample in the exception."""
        visited = []

        def render(sample, *args):
            visited.append(sample)
            if sample in {"A", "C"}:
                raise OSError(f"write denied for {sample}")

        with patch.object(pdf.os, "makedirs"), patch.object(
            pdf, "_generate_pdf_report_for_sample", side_effect=render
        ):
            with self.assertRaisesRegex(RuntimeError, "2 PDF report") as error:
                pdf.generate_all_reports(self.frame, None, None, Path("unused"), "test")
        self.assertEqual(visited, ["A", "B", "C"])
        self.assertIn("'A'", str(error.exception))
        self.assertIn("'C'", str(error.exception))
        self.assertNotIn("'B'", str(error.exception))

    def test_successful_batch_returns_normally(self):
        """All successful reports remain an ordinary successful call."""
        with patch.object(pdf.os, "makedirs"), patch.object(
            pdf, "_generate_pdf_report_for_sample"
        ) as render:
            self.assertIsNone(
                pdf.generate_all_reports(self.frame, None, None, Path("unused"), "test")
            )
        self.assertEqual([call.args[0] for call in render.call_args_list], ["A", "B", "C"])
