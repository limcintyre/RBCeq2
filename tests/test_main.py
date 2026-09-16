"""Test the main module's per-sample worker boundary."""

from pathlib import Path
import unittest
from unittest.mock import patch, sentinel

import pandas as pd

from rbceq2.core_logic.utils import BeyondLogicError
from rbceq2.main import SampleFailure, find_hits_or_failure


class TestFindHitsOrFailure(unittest.TestCase):
    """Preserve results, diagnostic context and process-control exceptions."""

    def test_success_preserves_result_and_forwards_arguments(self):
        source = Path("sample.vcf.gz")
        for expected in (sentinel.result, None):
            with self.subTest(result=expected):
                with patch("rbceq2.main.find_hits", return_value=expected) as find:
                    result = find_hits_or_failure(
                        sentinel.database, source, args=sentinel.args,
                        allele_relationships=sentinel.relationships,
                        excluded=sentinel.excluded, ant_mapping=sentinel.antigens,
                    )
                self.assertIs(result, expected)
                find.assert_called_once_with(
                    sentinel.database, source, args=sentinel.args,
                    allele_relationships=sentinel.relationships,
                    excluded=sentinel.excluded, ant_mapping=sentinel.antigens,
                )

    def test_failures_preserve_sample_identity_reason_and_traceback(self):
        named_error = BeyondLogicError(
            "ambiguous evidence", raised_by="test_worker/named_reason"
        )
        cases = (
            (Path("sample.vcf"), "sample", ValueError("invalid input"),
             "ValueError: invalid input"),
            (Path("sample.vcf.gz"), "sample.vcf", ValueError("invalid input"),
             "ValueError: invalid input"),
            ((pd.DataFrame(), "exact.sample-id"), "exact.sample-id", named_error,
             "BeyondLogicError: [test_worker/named_reason] ambiguous evidence"),
        )
        for source, sample, error, expected_message in cases:
            with self.subTest(sample=sample):
                with patch("rbceq2.main.find_hits", side_effect=error):
                    result = find_hits_or_failure(sentinel.database, source)
                self.assertIsInstance(result, SampleFailure)
                self.assertEqual(result.sample, sample)
                self.assertEqual(result.error, expected_message)
                self.assertIn("Traceback (most recent call last):", result.traceback)
                self.assertIn(expected_message, result.traceback)

    def test_cancellation_and_exit_propagate(self):
        for signal in (KeyboardInterrupt("cancelled"), SystemExit(2)):
            with self.subTest(exception=type(signal).__name__):
                with patch("rbceq2.main.find_hits", side_effect=signal):
                    with self.assertRaises(type(signal)) as caught:
                        find_hits_or_failure(sentinel.database, Path("sample.vcf"))
                self.assertIs(caught.exception, signal)


if __name__ == "__main__":
    unittest.main()
