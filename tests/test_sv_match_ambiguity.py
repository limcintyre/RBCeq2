"""Regression coverage for selecting evidence among equally ranked SV events."""

from dataclasses import replace
from itertools import permutations
import unittest

import pandas as pd

from rbceq2.core_logic.large_variants import SvDef, SvEvent, SvMatcher, SvReader
from rbceq2.core_logic.utils import BeyondLogicError


class TestSvMatchAmbiguity(unittest.TestCase):
    """Keep geometric selection independent of input order and sample metrics."""

    def setUp(self):
        self.db = SvDef("6", 100000, "DEL", 1000, "100000_del_1kb", id="I*test")
        self.matcher = SvMatcher()
        self.event = SvEvent(
            chrom="6", pos=100000, end=101000, svtype="DEL", svlen=-1000,
            alt="<DEL>", id="call1", qual="40", variant="6:100000_A_<DEL>",
            info={}, sample_fmt="GT:PS:DP:GQ", sample_value="0/1:.:12:30",
        )

    def assert_refused_in_both_orders(self, first, second):
        """Assert a named refusal regardless of which tied event arrives first."""
        for events in ((first, second), (second, first)):
            with self.subTest(samples=[event.sample_value for event in events]):
                with self.assertRaises(BeyondLogicError) as caught:
                    self.matcher.match([self.db], events)
                self.assertIn(
                    "SvMatcher.match/ambiguous_equal_best_sv_evidence",
                    str(caught.exception),
                )
                self.assertIn(self.db.raw, str(caught.exception))

    def assert_equivalent_in_both_orders(self, first, second):
        """Allow duplicate evidence to have either equivalent representative."""
        for events in ((first, second), (second, first)):
            with self.subTest(samples=[event.sample_value for event in events]):
                result = self.matcher.match([self.db], events)
                self.assertEqual(len(result), 1)
                self.assertIn(result[0].vcf, events)

    def test_conflicting_genotypes_refuse(self):
        self.assert_refused_in_both_orders(
            self.event, replace(self.event, sample_value="1/1:.:12:30")
        )

    def test_different_ploidy_is_not_collapsed(self):
        self.assert_refused_in_both_orders(
            replace(self.event, sample_value="1/1:.:12:30"),
            replace(self.event, sample_value="1:.:12:30"),
        )

    def test_opposite_orientation_in_same_phase_set_refuses(self):
        self.assert_refused_in_both_orders(
            replace(self.event, sample_value="0|1:100:12:30"),
            replace(self.event, sample_value="1|0:100:12:30"),
        )

    def test_different_phase_linkage_refuses_selection(self):
        # Different sets need not contradict; neither linkage may silently win.
        self.assert_refused_in_both_orders(
            replace(self.event, sample_value="0|1:100:12:30"),
            replace(self.event, sample_value="0|1:200:12:30"),
        )

    def test_filter_difference_refuses(self):
        self.assert_refused_in_both_orders(
            replace(self.event, filter_value="PASS"),
            replace(self.event, filter_value="LowQual"),
        )

    def test_metrics_and_format_order_are_not_evidence_conflicts(self):
        self.assert_equivalent_in_both_orders(
            self.event,
            replace(
                self.event, id="call2", qual="99", sample_fmt="GT:GQ:DP:PS",
                sample_value="0/1:99:100:.", info={"UNRELATED": "value"},
            ),
        )

    def test_unphased_allele_order_is_equivalent(self):
        self.assert_equivalent_in_both_orders(
            self.event, replace(self.event, sample_value="1/0:.:12:30")
        )

    def test_unphased_phase_set_is_ignored(self):
        self.assert_equivalent_in_both_orders(
            replace(self.event, sample_value="0/1:100:12:30"),
            replace(self.event, sample_value="0/1:200:12:30"),
        )

    def test_filter_order_and_repetition_are_equivalent(self):
        self.assert_equivalent_in_both_orders(
            replace(self.event, filter_value="LowQual;cnvLength"),
            replace(self.event, filter_value="cnvLength;LowQual;LowQual"),
        )

    def test_distinct_equally_ranked_physical_events_refuse(self):
        left = replace(self.event, pos=99990, end=100990, variant="left")
        right = replace(self.event, pos=100010, end=101010, variant="right")
        self.assertEqual(
            self.matcher._score_unclamped(self.db, left),
            self.matcher._score_unclamped(self.db, right),
        )
        self.assert_refused_in_both_orders(left, right)

    def test_worse_conflict_is_superseded_by_unique_best(self):
        worse = replace(self.event, pos=100010, end=101010, variant="worse")
        conflict = replace(worse, sample_value="1/1:.:12:30")
        for events in permutations((worse, conflict, self.event)):
            with self.subTest(order=[event.sample_value for event in events]):
                result = self.matcher.match([self.db], events)
                self.assertEqual(len(result), 1)
                self.assertIs(result[0].vcf, self.event)

    def test_public_zero_score_does_not_erase_geometric_winner(self):
        worse = replace(
            self.event, pos=100010, end=101010, variant="worse",
            sample_value="1/1:.:12:30",
        )
        self.assertEqual(self.matcher.score(self.db, self.event)[0], 0)
        self.assertEqual(self.matcher.score(self.db, worse)[0], 0)
        for events in ((worse, self.event), (self.event, worse)):
            result = self.matcher.match([self.db], events)
            self.assertEqual(len(result), 1)
            self.assertIs(result[0].vcf, self.event)


class TestSvReaderFilterEvidence(unittest.TestCase):
    """Retain FILTER evidence while supporting existing frames without FILTER."""

    def frame(self):
        return pd.DataFrame([{
            "CHROM": "6", "POS": 100000, "ID": "call1", "REF": "A",
            "ALT": "<DEL>", "QUAL": "40",
            "INFO": "SVTYPE=DEL;SVLEN=-1000;END=101000",
            "FORMAT": "GT:PS", "SAMPLE": "0|1:100",
            "variant": "6:100000_A_<DEL>",
        }])

    def test_filter_is_propagated(self):
        frame = self.frame()
        frame["FILTER"] = "LowQual;cnvLength"
        events = list(SvReader(frame).events())
        self.assertEqual(len(events), 1)
        self.assertEqual(events[0].filter_value, "LowQual;cnvLength")

    def test_absent_filter_defaults_to_dot(self):
        events = list(SvReader(self.frame()).events())
        self.assertEqual(len(events), 1)
        self.assertEqual(events[0].filter_value, ".")


if __name__ == "__main__":
    unittest.main()
