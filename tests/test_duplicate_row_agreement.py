"""Duplicate rows agree when they make the same claim, whatever order they write it in.

`rows_make_the_same_call` decides whether the rows a file gives more than once for one
token can be collapsed. Where any row is unphased there is no order being claimed on
both sides, and these cases pin down that the answer does not depend on which side of
the bar the alternate happens to sit.

The rest of these cases are about which orders are comparable at all. An order is
relative to its phase set, so two phased rows contradict each other only inside one
set; across two sets they are no more comparable than a phased row and an unphased one.
Rows naming no set are compared with each other, since the specification puts phased
genotypes carrying no PS in one common set, but not with a row that names one.
"""

import unittest

from rbceq2.IO.vcf import rows_make_the_same_call


class TestDuplicateRowAgreement(unittest.TestCase):
    def test_phased_and_unphased_agree_in_either_order(self):
        """These two say the same thing and used to answer differently."""
        self.assertTrue(rows_make_the_same_call(["0|1", "0/1"], ["100", ""]))
        self.assertTrue(rows_make_the_same_call(["1|0", "0/1"], ["100", ""]))
        self.assertTrue(rows_make_the_same_call(["0/1", "1|0"], ["", "100"]))

    def test_two_phased_rows_in_one_set_in_opposite_orders_disagree(self):
        """Inside one set that is a real contradiction and nothing here resolves it."""
        self.assertFalse(rows_make_the_same_call(["0|1", "1|0"], ["100", "100"]))

    def test_phased_rows_naming_no_set_are_compared_with_each_other(self):
        """The specification puts phased genotypes with no PS in one common set."""
        self.assertFalse(rows_make_the_same_call(["0|1", "1|0"], ["", ""]))

    def test_two_phased_rows_in_different_sets_are_not_comparable(self):
        """The defect: a set's orientation is arbitrary, so this was a coin toss.

        Both rows report one alternate copy. Their bars are relative to different
        sets, so neither confirms nor contradicts the other, and the copies decide.
        """
        for orders in (["0|1", "1|0"], ["1|0", "0|1"], ["0|1", "0|1"]):
            with self.subTest(orders=orders):
                self.assertTrue(rows_make_the_same_call(orders, ["100", "200"]))

    def test_an_unnamed_set_is_not_compared_with_a_named_one(self):
        self.assertTrue(rows_make_the_same_call(["0|1", "1|0"], ["100", ""]))
        self.assertTrue(rows_make_the_same_call(["0|1", "1|0"], ["", "100"]))

    def test_a_contradiction_inside_one_set_survives_a_third_row_elsewhere(self):
        """One set disagreeing with itself is not excused by a row from another."""
        self.assertFalse(
            rows_make_the_same_call(["0|1", "1|0", "0|1"], ["100", "100", "200"])
        )

    def test_agreeing_sets_with_a_third_row_elsewhere_still_reconcile(self):
        self.assertTrue(
            rows_make_the_same_call(["0|1", "0|1", "1|0"], ["100", "100", "200"])
        )

    def test_identical_calls_agree(self):
        for gts in (["0|1", "0|1"], ["1|0", "1|0"], ["0/1", "0/1"], ["1/1", "1|1"]):
            with self.subTest(gts=gts):
                self.assertTrue(rows_make_the_same_call(gts, ["100", "200"]))

    def test_differing_dosage_is_a_disagreement(self):
        """Different sets do not make different copy numbers agree."""
        self.assertFalse(rows_make_the_same_call(["0|1", "1/1"], ["100", "200"]))
        self.assertFalse(rows_make_the_same_call(["1|0", "0/0"], ["100", "200"]))
        self.assertFalse(rows_make_the_same_call(["0|1", "1|1"], ["100", "200"]))

    def test_a_no_call_is_never_reconciled(self):
        for sets in (["", ""], ["100", "200"]):
            with self.subTest(phase_sets=sets):
                self.assertFalse(rows_make_the_same_call(["0/1", "./."], sets))
                self.assertFalse(rows_make_the_same_call(["1|0", "."], sets))

    def test_all_reference_and_all_alternate_still_agree_across_ploidy(self):
        self.assertTrue(rows_make_the_same_call(["1/1/1/1", "1/1"], ["", ""]))
        self.assertTrue(rows_make_the_same_call(["0/0/0", "0|0"], ["", "100"]))

    def test_multi_copy_orders_agree_when_a_row_is_unphased(self):
        """Same two copies, one row naming an order and one not."""
        self.assertTrue(rows_make_the_same_call(["1|0", "0/1", "0/1"], ["", "", ""]))
        self.assertFalse(rows_make_the_same_call(["1|0", "0|1", "0/1"], ["", "", ""]))
        self.assertTrue(
            rows_make_the_same_call(["1|0", "0|1", "0/1"], ["100", "200", ""])
        )


if __name__ == "__main__":
    unittest.main()
