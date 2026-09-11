"""Duplicate rows agree when they make the same claim, whatever order they write it in.

`rows_make_the_same_call` decides whether the rows a file gives more than once for one
token can be collapsed. It is given the genotypes and nothing else, so it cannot tell a
contradiction inside one phase set from two orders that are relative to different sets.
Where any row is unphased there is no order being claimed on both sides, and these
cases pin down that the answer no longer depends on which side of the bar the alternate
happens to sit.
"""

import unittest

from rbceq2.IO.vcf import rows_make_the_same_call


class TestDuplicateRowAgreement(unittest.TestCase):
    def test_phased_and_unphased_agree_in_either_order(self):
        """The defect: these two say the same thing and answered differently."""
        self.assertTrue(rows_make_the_same_call(["0|1", "0/1"]))
        self.assertTrue(rows_make_the_same_call(["1|0", "0/1"]))
        self.assertTrue(rows_make_the_same_call(["0/1", "1|0"]))

    def test_two_phased_rows_in_opposite_orders_are_still_a_disagreement(self):
        """Within one set that is a real contradiction, and the sets are not visible."""
        self.assertFalse(rows_make_the_same_call(["0|1", "1|0"]))

    def test_identical_calls_agree(self):
        for gts in (["0|1", "0|1"], ["1|0", "1|0"], ["0/1", "0/1"], ["1/1", "1|1"]):
            with self.subTest(gts=gts):
                self.assertTrue(rows_make_the_same_call(gts))

    def test_differing_dosage_is_a_disagreement(self):
        self.assertFalse(rows_make_the_same_call(["0|1", "1/1"]))
        self.assertFalse(rows_make_the_same_call(["1|0", "0/0"]))

    def test_a_no_call_is_never_reconciled(self):
        self.assertFalse(rows_make_the_same_call(["0/1", "./."]))
        self.assertFalse(rows_make_the_same_call(["1|0", "."]))

    def test_all_reference_and_all_alternate_still_agree_across_ploidy(self):
        self.assertTrue(rows_make_the_same_call(["1/1/1/1", "1/1"]))
        self.assertTrue(rows_make_the_same_call(["0/0/0", "0|0"]))

    def test_multi_copy_orders_agree_when_a_row_is_unphased(self):
        """Same two copies, one row naming an order and one not."""
        self.assertTrue(rows_make_the_same_call(["1|0", "0/1", "0/1"]))
        self.assertFalse(rows_make_the_same_call(["1|0", "0|1", "0/1"]))


if __name__ == "__main__":
    unittest.main()
