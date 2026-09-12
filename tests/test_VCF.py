#!/usr/bin/env python3
"""
Unit tests for the VCF module.
"""

import gzip
import io
import os
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import pandas as pd
import polars as pl

from rbceq2.IO.vcf import (
    Interval,
    PloidyScan,
    VCF,
    VcfMissingHeaderError,
    check_if_multi_sample_vcf,
    read_vcf,
    split_vcf_to_dfs,
)
from rbceq2.core_logic.data_procesing import get_ref
from rbceq2.core_logic.utils import BeyondLogicError, Zygosity

# Dummy common columns list
COMMON_COLS = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]


class TestCheckIfMultiSampleVcf(unittest.TestCase):
    """Tests for sample-count detection and duplicate cohort header names."""

    @staticmethod
    def _header(samples: list[str]) -> str:
        """Return VCF metadata and a header with the supplied sample names."""
        return "##fileformat=VCFv4.2\n#" + "\t".join(COMMON_COLS + samples) + "\n"

    def test_duplicate_sample_names_raise_header_error(self) -> None:
        """Reject duplicate names through both text-opening branches."""
        content = self._header(["Z", "A", "Z", "A"])
        for suffix, opener in (
            (".vcf", "builtins.open"),
            (".vcf.gz", "rbceq2.IO.vcf.gzip.open"),
        ):
            with self.subTest(suffix=suffix):
                file_path = f"cohort{suffix}"
                with patch(opener, return_value=io.StringIO(content)) as mock_open:
                    with self.assertRaises(VcfMissingHeaderError) as context:
                        check_if_multi_sample_vcf(file_path)
                mock_open.assert_called_once_with(file_path, "rt")
                self.assertEqual(context.exception.filename, file_path)
                self.assertEqual(
                    context.exception.reason, "Duplicate VCF column names: A, Z"
                )

    def test_sample_name_matching_fixed_column_raises_header_error(self) -> None:
        """Keep duplicate detection across all columns of a cohort header."""
        content = self._header(["CHROM", "OTHER"])
        with patch("builtins.open", return_value=io.StringIO(content)):
            with self.assertRaises(VcfMissingHeaderError) as context:
                check_if_multi_sample_vcf("cohort.vcf")
        self.assertEqual(
            context.exception.reason, "Duplicate VCF column names: CHROM"
        )

    def test_unique_multi_sample_header_returns_true(self) -> None:
        """Continue to recognize a cohort with distinct sample names."""
        content = self._header(["FIRST", "SECOND"])
        with patch("builtins.open", return_value=io.StringIO(content)):
            self.assertTrue(check_if_multi_sample_vcf("cohort.vcf"))

    def test_single_sample_header_returns_false(self) -> None:
        """Preserve the ten-column branch, including its validation scope."""
        for sample in ("SINGLE", "CHROM"):
            with self.subTest(sample=sample):
                content = self._header([sample])
                with patch("builtins.open", return_value=io.StringIO(content)):
                    self.assertFalse(check_if_multi_sample_vcf("single.vcf"))

    def test_header_without_sample_column_raises_header_error(self) -> None:
        """Preserve the existing error for headers with fewer than ten columns."""
        content = self._header([])
        with patch("builtins.open", return_value=io.StringIO(content)):
            with self.assertRaises(VcfMissingHeaderError) as context:
                check_if_multi_sample_vcf("missing_sample.vcf")
        self.assertIsNone(context.exception.reason)


class TestConsumedFormat(unittest.TestCase):
    """Check FORMAT where retained rows first supply genotype information."""

    @staticmethod
    def _frame(format_field="GT:DP", sample_field="0/1:30", chrom="chr6", pos="100"):
        """Return one plain substitution with the supplied FORMAT and sample."""
        values = [
            chrom, pos, ".", "A", "G", ".", "PASS", ".", format_field, sample_field
        ]
        return pd.DataFrame([values], columns=COMMON_COLS + ["SAMPLE"])

    @staticmethod
    def _vcf(frame):
        """Construct one sample using a retained locus and no lane synthesis."""
        return VCF(
            frame, {}, {"6:100"}, sample="test_sample", reference_genome="GRCh38"
        )

    def test_retained_misplaced_gt_is_refused_before_ploidy_and_removal(self) -> None:
        """Neither a depth of 30 nor a depth of zero can become a genotype."""
        for sample_field in ("30:0/1", "0:1/1"):
            for input_type in ("pandas", "polars"):
                with self.subTest(sample_field=sample_field, input_type=input_type):
                    frame = self._frame("DP:GT", sample_field)
                    original = frame.copy(deep=True)
                    source = (
                        [frame] if input_type == "pandas" else pl.from_pandas(frame)
                    )
                    with self.assertRaises(BeyondLogicError) as context:
                        self._vcf(source)
                    self.assertEqual(context.exception.raised_by, "VCF/GT_not_first")
                    self.assertIn("6:100", context.exception.context)
                    self.assertIn("FORMAT='DP:GT'", context.exception.context)
                    pd.testing.assert_frame_equal(frame, original)

    def test_retained_gt_absence_is_explicitly_unsupported(self) -> None:
        """GT-less rows at an inference locus cannot become reference calls."""
        for format_field in ("DP", "GTX:DP", "."):
            with self.subTest(format_field=format_field):
                with self.assertRaises(BeyondLogicError) as context:
                    self._vcf([self._frame(format_field, "30")])
                self.assertEqual(context.exception.raised_by, "VCF/missing_GT")

    def test_gt_prefix_does_not_satisfy_the_first_key_requirement(self) -> None:
        """GTX preceding GT must be refused just like DP preceding GT."""
        with self.assertRaises(BeyondLogicError) as context:
            self._vcf([self._frame("GTX:GT", "30:0/1")])
        self.assertEqual(context.exception.raised_by, "VCF/GT_not_first")

    def test_split_refuses_unsupported_format_before_yield(self) -> None:
        """Standalone splitting cannot read depth as GT for its ploidy log."""
        for format_field, expected_reason in (
            ("DP:GT", "VCF/GT_not_first"),
            ("GTX:GT", "VCF/GT_not_first"),
            ("DP", "VCF/missing_GT"),
        ):
            with self.subTest(format_field=format_field):
                frame = self._frame(format_field, "30:0/1")
                frame["SECOND"] = "0:1/1"
                with patch("rbceq2.IO.vcf.logger.info") as log_info:
                    with self.assertRaises(BeyondLogicError) as context:
                        next(split_vcf_to_dfs(frame))
                self.assertEqual(context.exception.raised_by, expected_reason)
                log_info.assert_not_called()

    def test_split_without_rows_or_samples_does_not_validate_format(self) -> None:
        """No row or no sample leaves no genotype for the splitter to consume."""
        frame = self._frame("DP", "30")
        self.assertEqual(list(split_vcf_to_dfs(frame[COMMON_COLS])), [])
        split_frames = list(split_vcf_to_dfs(frame.iloc[:0]))
        self.assertEqual(len(split_frames), 1)
        self.assertTrue(split_frames[0][0].empty)

    def test_direct_get_variants_refuses_misplaced_gt(self) -> None:
        """Direct helper callers receive the same named input refusal."""
        vcf = self._vcf([self._frame()])
        vcf.df.loc[0, "FORMAT"] = "DP:GT"
        vcf.df.loc[0, "SAMPLE"] = "30:0/1"
        with self.assertRaises(BeyondLogicError) as context:
            vcf.get_variants()
        self.assertEqual(context.exception.raised_by, "VCF/GT_not_first")
        self.assertIn("6:100_A_G", context.exception.context)

    def test_gt_only_and_reordered_trailing_fields_keep_their_meaning(self) -> None:
        """The guard preserves GT alone, metrics, phasing and ordinary ploidy."""
        for format_field, sample_field, expected in (
            ("GT", "0/1", {"GT": "0/1"}),
            ("GT:GQ:PS", "0|1:99:150", {"GT": "0|1", "GQ": "99", "PS": "150"}),
            ("GT:PS:GQ", "0|1:150:99", {"GT": "0|1", "PS": "150", "GQ": "99"}),
        ):
            with self.subTest(format_field=format_field):
                vcf = self._vcf([self._frame(format_field, sample_field)])
                self.assertEqual(vcf.variants["6:100_A_G"], expected)
                self.assertEqual(vcf.haploid_chroms, frozenset())
                self.assertEqual(vcf.diploid_loci["6"], frozenset({100}))

    def test_ignored_autosomal_and_unused_phase_neighbor_rows_remain_ignored(self) -> None:
        """Raw chr9 PS candidates need validation only if final selection uses them."""
        frame = pd.concat(
            [
                self._frame(),
                self._frame("DP:GT", "30:0/1", "chr1", "50"),
                self._frame("DP:GT:PS", "0:1/1:100", "chr9", "100"),
                self._frame("GT:PS", "0|1:200", "chr9", "200"),
                self._frame("GT:PS", "0|1:300", "chr9", "300"),
            ],
            ignore_index=True,
        )
        content = "##fileformat=VCFv4.2\n#" + frame.to_csv(sep="\t", index=False)
        scan = PloidyScan("GRCh38")
        with patch("builtins.open", return_value=io.StringIO(content)):
            raw_frame = read_vcf(
                "sample.vcf", {"6": [Interval(1, 1000)], "9": [Interval(1, 1000)]},
                unique_variants={"6:100"}, ploidy_scan=scan,
            )
        self.assertEqual(raw_frame.height, 4)
        vcf = self._vcf(raw_frame)
        self.assertEqual(
            set(vcf.variants), {"6:100_A_G", "9:200_A_G", "9:300_A_G"}
        )
        self.assertEqual(scan.for_sample("SAMPLE"), frozenset())

    def test_discarded_gt_absent_non_par_row_supplies_no_evidence(self) -> None:
        """A discarded depth-only X row has the same result as omitting it."""
        expected = self._vcf([self._frame()])
        for format_field, sample_field in (("DP", "30"), ("GTX:DP", "1:30")):
            with self.subTest(format_field=format_field):
                frame = pd.concat(
                    [
                        self._frame(),
                        self._frame(format_field, sample_field, "chrX", "10000000"),
                    ],
                    ignore_index=True,
                )
                content = "##fileformat=VCFv4.2\n#" + frame.to_csv(sep="\t", index=False)
                scan = PloidyScan("GRCh38")
                with patch("builtins.open", return_value=io.StringIO(content)):
                    raw_frame = read_vcf(
                        "sample.vcf", {"6": [Interval(1, 1000)]},
                        unique_variants={"6:100"}, ploidy_scan=scan,
                    )
                actual = self._vcf(raw_frame)
                self.assertEqual(scan.for_sample("SAMPLE"), frozenset())
                self.assertEqual(actual.variants, expected.variants)
                self.assertEqual(actual.haploid_chroms, expected.haploid_chroms)


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

    def test_synthesised_lane_row_has_no_filter_value(self) -> None:
        """A lane locus with no called row is RBCeq2's own assertion, not a call.

        Its ID, REF, ALT, QUAL, FILTER and INFO were filled with the *column names*, so
        the FILTER column of every such row held the literal string 'FILTER'. Nothing
        read it - only '_ref' tokens carry it and the FILTER lookup skips those - but it
        made the column impossible to survey and record_unused_variants had to explain
        it. '.' is what the specification uses for a field with no value.
        """
        vcf_obj = VCF(
            [self.df_local], {"2": ["9999"]}, set(), sample="test_sample"
        )
        synthesised = vcf_obj.df[vcf_obj.df["variant"] == "2:9999_ref"]
        self.assertEqual(len(synthesised), 1)
        row = synthesised.iloc[0]
        for column in ("ID", "REF", "ALT", "QUAL", "FILTER", "INFO"):
            self.assertEqual(
                row[column], ".", f"{column} should carry no value, not its own name"
            )

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
        gts: list[str], filters: list[str], phase_sets: list[str] | None = None
    ) -> pd.DataFrame:
        """Return a VCF-like DataFrame whose rows are all the same variant.

        The shape a file takes when two callers both report a position and FILTER marks
        which of them is the one that conflicts. One position, one REF/ALT, so every row
        encodes to the same token.

        Args:
            gts (list[str]): One GT per row, in file order.
            filters (list[str]): One FILTER value per row, in file order.
            phase_sets (list[str] | None): One PS per row, in file order. Defaults to
                the same set on every row, which is every case where the rows cannot
                name different ones.

        Returns:
            pd.DataFrame: One row per genotype, all T>C at chr2:2000.
        """
        n = len(gts)
        sets = ["1"] * n if phase_sets is None else phase_sets
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
                "SAMPLE": [f"{gt}:1,1:30:30:{ps}" for gt, ps in zip(gts, sets)],
            }
        )

    def test_get_variants_warns_when_rows_contradict_each_other(self) -> None:
        """Two rows, one token, genotypes that cannot both be right.

        Rows that agree are reconciled before get_variants runs, so what reaches it is
        the case reconciliation declined to arbitrate: the rows disagree about the call
        and nothing at this layer can know which caller to believe. The dict keeps the
        last, so the genotype in use is whichever the file ended with, and this is the
        only place that loss can be reported.
        """
        with patch("rbceq2.IO.vcf.logger") as mock_logger:
            vcf_obj = VCF(
                [
                    self._df_with_two_rows_for_one_variant(
                        ["0|1", "1/1"], ["PASS", "TargetedConflict"]
                    )
                ],
                {},
                set(),
                sample="test_sample",
            )
        self.assertEqual(vcf_obj.variants["2:2000_T_C"]["GT"], "1/1")
        warnings = [call.args[0] for call in mock_logger.warning.call_args_list]
        contradiction = [w for w in warnings if "more than one row" in w]
        self.assertEqual(len(contradiction), 1)
        self.assertIn("2:2000_T_C (0|1 then 1/1)", contradiction[0])
        self.assertIn("test_sample", contradiction[0])

    def test_get_variants_silent_when_the_rows_agree(self) -> None:
        """Rows that agree never reach get_variants - they were reconciled to one."""
        with patch("rbceq2.IO.vcf.logger") as mock_logger:
            vcf_obj = VCF(
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
        self.assertEqual(len(vcf_obj.df), 1)

    def test_reconciliation_keeps_the_phased_row(self) -> None:
        """An unphased genotype does not contradict a phased one - it is less specific.

        Discarding the phase because the unphased row happens to come second is the
        825 heterozygous calls this exists to stop losing.
        """
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
        self.assertEqual(len(vcf_obj.df), 1)
        self.assertEqual(vcf_obj.variants["2:2000_T_C"]["GT"], "0|1")

    def test_reconciliation_joins_every_FILTER_value(self) -> None:
        """A FILTER field already carries several values ';' joined, and the classifier
        already splits on it, so the verdict stops depending on which row is first."""
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
        self.assertEqual(vcf_obj.df["FILTER"].tolist(), ["PASS;TargetedConflict"])

    @staticmethod
    def _phased_rows_naming_sets(
        positions_and_sets: list[tuple[str, str]],
    ) -> pd.DataFrame:
        """Rows whose only job is to put weight behind a phase set.

        Each is a distinct variant at a distinct position, so none of them collides with
        the token under test - they are here to be counted, not reconciled.

        Args:
            positions_and_sets (list[tuple[str, str]]): One (POS, PS) per row.

        Returns:
            pd.DataFrame: One phased T>G row per pair, on chr2.
        """
        n = len(positions_and_sets)
        return pd.DataFrame(
            {
                "CHROM": ["chr2"] * n,
                "POS": [pos for pos, _ in positions_and_sets],
                "ID": ["."] * n,
                "REF": ["T"] * n,
                "ALT": ["G"] * n,
                "QUAL": ["."] * n,
                "FILTER": ["PASS"] * n,
                "INFO": ["."] * n,
                "FORMAT": ["GT:AD:GQ:DP:PS"] * n,
                "SAMPLE": [f"0|1:1,1:30:30:{ps}" for _, ps in positions_and_sets],
            }
        )

    def test_reconciliation_keeps_the_phase_set_the_most_rows_name(self) -> None:
        """The better attested set wins, even though its row is not the first.

        A phase set is spent a comparison at a time, so the set with more variants in it
        is the one that lets the phased filters do more: two variants offer one
        comparison and five offer ten. Here the second row's set is named by three rows
        and the first row's by one, so the second row is kept - which is the whole point
        of counting, since reading the file's order would have taken the first.
        """
        frame = pd.concat(
            [
                self._df_with_two_rows_for_one_variant(
                    ["0|1", "0|1"],
                    ["PASS", "TargetedConflict"],
                    phase_sets=["1500", "1900"],
                ),
                self._phased_rows_naming_sets([("2100", "1900"), ("2200", "1900")]),
            ],
            ignore_index=True,
        )
        vcf_obj = VCF([frame], {}, set(), sample="test_sample")
        self.assertEqual(vcf_obj.variants["2:2000_T_C"]["PS"], "1900")

    def test_reconciliation_prefers_a_row_that_names_a_set_to_one_that_does_not(
        self,
    ) -> None:
        """A phased row naming no set counts zero, so it loses to one that names a set.

        Phase with nothing to measure it against orients the call relative to nothing
        the file declares. The row that names a set is the more specific of the two, in
        the same way a phased genotype is more specific than an unphased one, so the
        same preference applies - and it has to survive the row naming no set coming
        first, which is what this pins.
        """
        frame = self._df_with_two_rows_for_one_variant(
            ["0|1", "0|1"], ["PASS", "TargetedConflict"], phase_sets=["1500", "1900"]
        )
        frame.loc[0, "FORMAT"] = "GT:AD:GQ:DP"
        frame.loc[0, "SAMPLE"] = "0|1:1,1:30:30"
        vcf_obj = VCF([frame], {}, set(), sample="test_sample")
        self.assertEqual(vcf_obj.variants["2:2000_T_C"]["PS"], "1900")

    def test_reconciliation_collapses_rows_that_name_different_phase_sets(self) -> None:
        """Two phase sets are not two calls.

        Both rows say the alternate is on one chromosome and the reference on the other.
        They differ only in which set that orientation is measured against, and two sets
        are by definition not phased relative to one another, so neither row is wrong
        and neither contradicts the other. rows_make_the_same_call compares orders only
        inside one set, so this collapses the way any other agreeing pair does.
        """
        vcf_obj = VCF(
            [
                self._df_with_two_rows_for_one_variant(
                    ["0|1", "0|1"],
                    ["PASS", "TargetedConflict"],
                    phase_sets=["1500", "1900"],
                )
            ],
            {},
            set(),
            sample="test_sample",
        )
        self.assertEqual(len(vcf_obj.df), 1)

    def test_reconciliation_collapses_opposite_orders_in_different_sets(self) -> None:
        """Which way round a set is written cannot decide whether two rows agree.

        The pair above with one set's orientation reversed, which is the same claim: a
        set says how its own members sit relative to each other and nothing about how it
        sits against another set. Both rows still report one alternate copy, so they
        still collapse, and the tie rule still keeps the first phased row.

        This is the shape a run emitting both a targeted caller's calls and the general
        caller's produces - the two rows carrying the orientation from a small local set
        and from the wider one. Reading the genotypes alone, this pair was refused as a
        contradiction while the pair above collapsed, so whether the duplicate survived
        depended on which way round a caller happened to write a set.
        """
        vcf_obj = VCF(
            [
                self._df_with_two_rows_for_one_variant(
                    ["0|1", "1|0"],
                    ["PASS", "TargetedConflict"],
                    phase_sets=["1500", "1900"],
                )
            ],
            {},
            set(),
            sample="test_sample",
        )
        self.assertEqual(len(vcf_obj.df), 1)
        self.assertEqual(vcf_obj.variants["2:2000_T_C"]["PS"], "1500")

    def test_reconciliation_falls_back_to_the_first_phased_row_on_a_tie(self) -> None:
        """Two sets each named by one row are equally attested, so order decides.

        The count cannot separate them and nothing else here can either, so the rule
        runs out and the fallback is what reconciliation did before the count existed.
        It has to be a stated fallback rather than an accident: max keeps the first of
        equal keys, and that is load bearing, not incidental.
        """
        vcf_obj = VCF(
            [
                self._df_with_two_rows_for_one_variant(
                    ["0|1", "0|1"],
                    ["PASS", "TargetedConflict"],
                    phase_sets=["1500", "1900"],
                )
            ],
            {},
            set(),
            sample="test_sample",
        )
        self.assertEqual(vcf_obj.variants["2:2000_T_C"]["PS"], "1500")

    def test_reconciliation_accepts_a_difference_of_ploidy_alone(self) -> None:
        """'1/1/1/1' and '1/1' are the same claim: every copy carries the alternate.

        A gene conversion caller reporting a paralogue pair writes the first and the
        general caller the second. Ploidy is per genotype in a VCF, so differing ploidy
        is not by itself a disagreement.
        """
        vcf_obj = VCF(
            [
                self._df_with_two_rows_for_one_variant(
                    ["1/1/1/1", "1/1"], ["PASS", "TargetedConflict"]
                )
            ],
            {},
            set(),
            sample="test_sample",
        )
        self.assertEqual(len(vcf_obj.df), 1)

    def test_reconciliation_leaves_disagreeing_rows_alone(self) -> None:
        """Deciding which caller to believe is not something this layer can know."""
        vcf_obj = VCF(
            [
                self._df_with_two_rows_for_one_variant(
                    ["0|1", "1/1"], ["PASS", "TargetedConflict"]
                )
            ],
            {},
            set(),
            sample="test_sample",
        )
        self.assertEqual(len(vcf_obj.df), 2)

    def test_reconciliation_will_not_touch_a_structural_row(self) -> None:
        """The structural reader takes its events off these rows, and its size threshold
        is a run time argument this layer never sees."""
        frame = self._df_with_two_rows_for_one_variant(
            ["0/1", "0/1"], ["PASS", "PASS"]
        )
        frame["REF"] = ["N", "N"]
        frame["ALT"] = ["<DEL>", "<DEL>"]
        frame["INFO"] = ["SVTYPE=DEL;END=2500", "SVTYPE=DEL;END=2500"]
        vcf_obj = VCF([frame], {}, set(), sample="test_sample")
        self.assertEqual(len(vcf_obj.df), 2)

    def test_remove_home_ref_drops_haploid_zero(self) -> None:
        """Haploid '0' IS hom ref, wherever it is.

        Reversed when a gene reported consistently at one copy started being read as one
        copy. It used to be kept so it would fail loudly downstream, on the grounds that
        whether the region is single copy was undecided. That turned out to
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



class TestHomRefBoundary(unittest.TestCase):
    """Preserve dosage and no-call evidence before variant interpretation."""

    LOCUS = "7:100893176"
    TOKEN = "7:100893176_G_T"

    @classmethod
    def _vcf(cls, gt, format_field, input_type, *, lanes=False):
        """Construct a sample at the curated YT*02 substitution."""
        sample_field = gt if format_field == "GT" else f"{gt}:30"
        frame = pd.DataFrame(
            [["chr7", "100893176", ".", "G", "T", "50", "PASS", ".",
              format_field, sample_field]],
            columns=COMMON_COLS + ["SAMPLE"],
        )
        source = [frame] if input_type == "pandas" else pl.from_pandas(frame)
        return VCF(
            source,
            {"chr7": ["100893176"]} if lanes else {},
            {cls.LOCUS},
            sample="hom_ref_boundary",
            reference_genome="GRCh38",
        )

    def test_all_reference_calls_are_removed_at_every_supported_ploidy(self):
        """Removing reference rows preserves the earlier locus-copy evidence."""
        for gt in ("0", "0/0", "0|0", "0/0/0", "0/0/0/0", "0|0|0|0"):
            for format_field in ("GT", "GT:DP"):
                for input_type in ("pandas", "polars"):
                    with self.subTest(gt=gt, format=format_field, input=input_type):
                        vcf = self._vcf(gt, format_field, input_type)
                        self.assertTrue(vcf.df.empty)
                        self.assertEqual(vcf.variants, {})
                        self.assertNotIn(self.LOCUS, vcf.loci)
                        if gt == "0":
                            self.assertIn(100893176, vcf.haploid_loci["7"])

    def test_intermediate_dosage_reaches_the_existing_named_refusal(self):
        """Reference prefixes cannot bypass dosage handling in either order."""
        for gt in ("0/0/1", "1/0/0", "0/0/0/1", "1/0/0/0", "0/0/1/1",
                   "0|0|1", "0|0|0|1"):
            for format_field in ("GT", "GT:DP"):
                for input_type in ("pandas", "polars"):
                    with self.subTest(gt=gt, format=format_field, input=input_type):
                        vcf = self._vcf(gt, format_field, input_type)
                        self.assertIn(self.LOCUS, vcf.loci)
                        self.assertIn(self.TOKEN, vcf.variants)
                        metrics = vcf.variants[self.TOKEN]
                        self.assertEqual(metrics["GT"], gt)
                        with self.assertRaises(BeyondLogicError) as caught:
                            get_ref(metrics, self.TOKEN)
                        self.assertEqual(
                            caught.exception.raised_by,
                            "get_ref/dosage_between_the_bounds",
                        )
                        self.assertIn(self.TOKEN, str(caught.exception))

    def test_partial_no_calls_keep_the_locus_without_synthesising_reference(self):
        """A missing allele must reach NO_DATA with its original GT intact."""
        for gt in ("0/0/.", "./0/0", "0|0|.", "./././.", ".|.|.|."):
            for format_field in ("GT", "GT:DP"):
                for input_type in ("pandas", "polars"):
                    with self.subTest(gt=gt, format=format_field, input=input_type):
                        vcf = self._vcf(gt, format_field, input_type, lanes=True)
                        self.assertIn(self.LOCUS, vcf.loci)
                        self.assertEqual(set(vcf.variants), {self.TOKEN})
                        metrics = vcf.variants[self.TOKEN]
                        self.assertEqual(metrics["GT"], gt)
                        self.assertEqual(get_ref(metrics, self.TOKEN), Zygosity.NO_DATA)

    def test_supported_variant_dosage_survives_the_sample_boundary(self):
        """GT-only and trailing metrics preserve the supported dosage paths."""
        for gt, expected, locus_copies in (
            ("1", Zygosity.HEM, 1),
            ("0/1", Zygosity.HET, None),
            ("0|1", Zygosity.HET, None),
            ("1/1", Zygosity.HOM, None),
            ("1/1/1/1", Zygosity.HOM, None),
            ("1|1|1|1", Zygosity.HOM, None),
        ):
            for format_field in ("GT", "GT:DP"):
                for input_type in ("pandas", "polars"):
                    with self.subTest(gt=gt, format=format_field, input=input_type):
                        vcf = self._vcf(gt, format_field, input_type)
                        metrics = vcf.variants[self.TOKEN]
                        expected_metrics = {"GT": gt}
                        if format_field == "GT:DP":
                            expected_metrics["DP"] = "30"
                        self.assertEqual(metrics, expected_metrics)
                        self.assertEqual(
                            get_ref(metrics, self.TOKEN, locus_copies=locus_copies),
                            expected,
                        )


class TestGtValidityBoundary(unittest.TestCase):
    """Only supported GTs naming declared alleles may supply evidence."""

    @staticmethod
    def _vcf(gt, alt="A,C", format_field="GT", input_type="pandas"):
        sample = gt if format_field == "GT" else f"{gt}:30"
        frame = pd.DataFrame({
            "CHROM": ["chr7"], "POS": ["142957921"], "ID": ["."],
            "REF": ["G"], "ALT": [alt], "QUAL": ["50"],
            "FILTER": ["PASS"], "INFO": ["."], "FORMAT": [format_field],
            "SAMPLE": [sample],
        })
        source = [frame] if input_type == "pandas" else pl.from_pandas(frame)
        return VCF(source, {}, {"7:142957921"}, sample="s",
                   reference_genome="GRCh38")

    @staticmethod
    def _scan_row(gt, alt="G", format_field="GT"):
        sample = gt if format_field == "GT" else f"{gt}:30"
        return ["X", "50000000", ".", "A", alt, "50", "PASS", ".",
                format_field, sample]

    def test_malformed_calls_are_refused_before_inference_or_recoding(self):
        for gt in ("garbage", "-1", "1.0", "1//2", "0/0x", "", "0/0/", "1 ", " 1"):
            for fmt in ("GT", "GT:DP"):
                for input_type in ("pandas", "polars"):
                    with self.subTest(gt=gt, format=fmt, input=input_type):
                        with self.assertRaises(BeyondLogicError) as caught:
                            self._vcf(gt, format_field=fmt, input_type=input_type)
                        self.assertEqual(caught.exception.raised_by, "VCF/invalid_GT")
                        self.assertIn("s", caught.exception.context)
                        self.assertIn("142957921", caught.exception.context)

    def test_indices_must_name_an_alt_declared_by_the_row(self):
        for gt, alt in (
            ("2", "A"), ("2/2", "A"), ("./2", "A"), ("0/0/2", "A"),
            ("3/3", "A,C"), ("1/3", "A,C"), ("./3", "A,C"),
            ("1/999999999999999999999", "A,C"), ("1", "."),
        ):
            for fmt in ("GT", "GT:DP"):
                with self.subTest(gt=gt, alt=alt, format=fmt):
                    with self.assertRaises(BeyondLogicError) as caught:
                        self._vcf(gt, alt, fmt)
                    self.assertEqual(
                        caught.exception.raised_by, "VCF/GT_index_out_of_range"
                    )
                    self.assertIn(gt, caught.exception.context)

    def test_unsupported_spellings_are_distinguished_from_malformed_calls(self):
        for gt in ("/0/1", "|1", "01/2"):
            with self.subTest(gt=gt):
                with self.assertRaises(BeyondLogicError) as caught:
                    self._vcf(gt)
                self.assertEqual(
                    caught.exception.raised_by, "VCF/unsupported_GT_encoding"
                )

    def test_multi_digit_indices_and_no_alt_reference_remain_readable(self):
        alts = "A,C,T,GA,GC,GT,GAA,GAC,GAT,GCA"
        vcf = self._vcf("10", alts)
        self.assertEqual(vcf.variants, {"7:142957921_G_GCA": {"GT": "1"}})
        self.assertEqual(self._vcf("0", ".").variants, {})

    def test_direct_multi_alt_extraction_cannot_bypass_range_validation(self):
        vcf = self._vcf("1/2")
        vcf.df.loc[:, "SAMPLE"] = "1/3"
        with self.assertRaises(BeyondLogicError) as caught:
            vcf.get_variants()
        self.assertEqual(caught.exception.raised_by, "VCF/GT_index_out_of_range")

    def test_bad_scan_candidates_warn_and_supply_no_copy_evidence(self):
        for gt, reason in (
            ("garbage", "VCF/invalid_GT"), ("-1", "VCF/invalid_GT"),
            ("1.0", "VCF/invalid_GT"), ("1 ", "VCF/invalid_GT"),
            (" 1", "VCF/invalid_GT"), ("2", "VCF/GT_index_out_of_range"),
            ("01", "VCF/unsupported_GT_encoding"),
        ):
            for fmt in ("GT", "GT:DP"):
                with self.subTest(gt=gt, format=fmt):
                    scan = PloidyScan("GRCh38")
                    with patch("rbceq2.IO.vcf.logger.warning") as warning:
                        scan.observe("X", 50000000, ["s"], self._scan_row(gt, format_field=fmt))
                    self.assertEqual(scan.for_sample("s"), frozenset())
                    warning.assert_called_once()
                    self.assertIn(reason, str(warning.call_args))
                    self.assertIn("50000000", str(warning.call_args))

    def test_valid_evidence_survives_bad_candidates_in_either_order(self):
        for gts in (("garbage", "0"), ("0", "garbage")):
            with self.subTest(gts=gts):
                scan = PloidyScan("GRCh38")
                with patch("rbceq2.IO.vcf.logger.warning") as warning:
                    for gt in gts:
                        scan.observe("X", 50000000, ["s"], self._scan_row(gt))
                self.assertEqual(scan.for_sample("s"), frozenset({"X"}))
                warning.assert_called_once()
                self.assertIn("VCF/invalid_GT", str(warning.call_args))

    def test_scan_gt_only_line_endings_are_not_part_of_the_genotype(self):
        for gt in ("1\n", "1\r\n"):
            with self.subTest(gt=gt):
                scan = PloidyScan("GRCh38")
                with patch("rbceq2.IO.vcf.logger.warning") as warning:
                    scan.observe("X", 50000000, ["s"], self._scan_row(gt))
                self.assertEqual(scan.for_sample("s"), frozenset({"X"}))
                warning.assert_not_called()

    def test_scan_index_bounds_and_samples_are_independent(self):
        scan = PloidyScan("GRCh38")
        row = self._scan_row("2", alt="G,T") + ["3", ".", "./.", "0"]
        with patch("rbceq2.IO.vcf.logger.warning") as warning:
            scan.observe("X", 50000000, ["valid", "bad", "missing", "missing2", "ref"], row)
        self.assertEqual(scan.for_sample("valid"), frozenset({"X"}))
        self.assertEqual(scan.for_sample("ref"), frozenset({"X"}))
        for sample in ("bad", "missing", "missing2"):
            self.assertEqual(scan.for_sample(sample), frozenset())
        warning.assert_called_once()
        self.assertIn("sample=bad", str(warning.call_args))


if __name__ == "__main__":
    unittest.main()
