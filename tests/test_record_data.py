#!/usr/bin/env python3
"""Unit tests for rbceq2.IO.record_data."""

import unittest

from loguru import logger

from rbceq2.IO.record_data import report_no_call_summary


class TestReportNoCallSummary(unittest.TestCase):
    """The run level warning about database loci the caller did not call."""

    def setUp(self) -> None:
        self.messages: list[str] = []
        self.sink_id = logger.add(
            lambda msg: self.messages.append(msg.record["message"]), level="WARNING"
        )

    def tearDown(self) -> None:
        logger.remove(self.sink_id)

    def test_silent_when_there_are_no_no_calls(self) -> None:
        """A clean run must not emit a warning at all - see the public_truth datasets."""
        report_no_call_summary({}, 17)
        self.assertEqual(self.messages, [])

    def test_silent_when_no_samples_were_processed(self) -> None:
        """Guards the denominator; nothing to report and nothing to divide by."""
        report_no_call_summary({("FY", "1:159174683_T_C"): {"s1"}}, 0)
        self.assertEqual(self.messages, [])

    def test_one_line_per_locus_plus_a_header(self) -> None:
        report_no_call_summary(
            {
                ("FY", "1:159174683_T_C"): {"s1", "s2"},
                ("VEL", "1:3691528_A_G"): {"s1"},
            },
            4,
        )
        self.assertEqual(len(self.messages), 3)
        self.assertIn("2 database loci were not called", self.messages[0])
        self.assertIn("no_call_at_defining_variant", self.messages[0])

    def test_loci_are_ranked_worst_first(self) -> None:
        """The dead probe is the point of the warning, so it goes at the top."""
        report_no_call_summary(
            {
                ("VEL", "1:3691528_A_G"): {"s1"},
                ("FY", "1:159174683_T_C"): {"s1", "s2", "s3"},
                ("KN", "1:207760773_C_T"): {"s1", "s2"},
            },
            4,
        )
        self.assertIn("FY", self.messages[1])
        self.assertIn("KN", self.messages[2])
        self.assertIn("VEL", self.messages[3])

    def test_rate_is_reported_against_samples_processed(self) -> None:
        report_no_call_summary({("FY", "1:159174683_T_C"): {"s1", "s2", "s3"}}, 4)
        self.assertIn("no call in 3/4 samples (75.0%)", self.messages[1])

    def test_majority_failure_is_called_out(self) -> None:
        """Above half, the probe is the problem rather than the samples."""
        report_no_call_summary({("FY", "1:159174683_T_C"): {"s1", "s2", "s3"}}, 4)
        self.assertIn("failed in the majority of samples", self.messages[1])

    def test_minority_failure_is_not_called_out(self) -> None:
        report_no_call_summary({("FY", "1:159174683_T_C"): {"s1"}}, 4)
        self.assertNotIn("failed in the majority", self.messages[1])

    def test_exactly_half_is_not_a_majority(self) -> None:
        """Boundary: the marker is for '> 50%', not '>= 50%'."""
        report_no_call_summary({("FY", "1:159174683_T_C"): {"s1", "s2"}}, 4)
        self.assertIn("(50.0%)", self.messages[1])
        self.assertNotIn("failed in the majority", self.messages[1])


if __name__ == "__main__":
    unittest.main()
