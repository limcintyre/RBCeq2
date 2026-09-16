"""Phase-set provenance: a block's span is not evidence of membership.

`assign_ref_phase_set` fills in the phase set of a reference token, and of any token
whose PS field the input omitted. Three sources are legitimate, in order: the token's
own row, since a lane reference row is a copy of the measured row and carries its PS;
the row of an alternate measured at the same position, which describes the same locus;
and, failing both, nothing. What it may not do is answer from a position merely lying
inside another block's span, which is what these cases pin down.
"""

import unittest

from rbceq2.core_logic.alleles import BloodGroup
from rbceq2.core_logic.data_procesing import add_phasing
from rbceq2.core_logic.utils import Zygosity
from rbceq2.filters.shared_filter_functionality import carries_phase

LOW = "18:45730438_G_A"
MID_REF = "18:45736484_ref"
MID_ALT = "18:45736484_A_G"
HIGH = "18:45748385_C_T"
# One block, declared by the two outer tokens, whose span covers the middle position.
SPANNING_BLOCK = {"18": {100: (45730438, 45748385)}}


class TestPhaseSetProvenance(unittest.TestCase):
    def run_add_phasing(self, pool, metrics, phase_sets=SPANNING_BLOCK):
        bg = BloodGroup(type="JK", sample="provenance", alleles={},
                        variant_pool=dict(pool))
        add_phasing({"JK": bg}, phased=True, variant_metrics=metrics,
                    phase_sets=phase_sets)
        return bg

    def test_containment_alone_does_not_name_a_block(self):
        """A reference token inside a block's span borrows nothing from it."""
        pool = {LOW: Zygosity.HET, MID_REF: Zygosity.HET, HIGH: Zygosity.HET}
        metrics = {LOW: {"GT": "0|1", "PS": "100"},
                   MID_REF: {"GT": "1|0"},
                   HIGH: {"GT": "0|1", "PS": "100"}}
        bg = self.run_add_phasing(pool, metrics)
        self.assertEqual(bg.variant_pool_phase_set[MID_REF], "unknown")
        self.assertEqual(bg.variant_pool_phase_set[LOW], "100")
        self.assertEqual(bg.variant_pool_phase_set[HIGH], "100")

    def test_same_position_partner_is_still_borrowed(self):
        """The partner borrow is a statement about one locus and is unaffected."""
        pool = {LOW: Zygosity.HET, MID_REF: Zygosity.HET, MID_ALT: Zygosity.HET}
        metrics = {LOW: {"GT": "0|1", "PS": "100"},
                   MID_REF: {"GT": "1|0"},
                   MID_ALT: {"GT": "0|1", "PS": "250"}}
        bg = self.run_add_phasing(pool, metrics)
        self.assertEqual(bg.variant_pool_phase_set[MID_REF], "250")

    def test_partner_without_a_phase_set_still_declines(self):
        """A partner that names no set cannot lend one, and the span may not fill in."""
        pool = {LOW: Zygosity.HET, MID_REF: Zygosity.HET, MID_ALT: Zygosity.HET}
        metrics = {LOW: {"GT": "0|1", "PS": "100"},
                   MID_REF: {"GT": "1|0"},
                   MID_ALT: {"GT": "0|1"}}
        bg = self.run_add_phasing(pool, metrics)
        self.assertEqual(bg.variant_pool_phase_set[MID_REF], "unknown")

    def test_homozygous_reference_token_is_unchanged(self):
        """The homozygous short circuit runs before any of this."""
        pool = {LOW: Zygosity.HET, MID_REF: Zygosity.HOM}
        metrics = {LOW: {"GT": "0|1", "PS": "100"}, MID_REF: {"GT": "1/1"}}
        bg = self.run_add_phasing(pool, metrics)
        self.assertEqual(bg.variant_pool_phase_set[MID_REF], ".")

    def test_neither_spelling_of_no_set_is_answered_with_a_block(self):
        """A reported "." stays "."; a field the row omitted becomes "unknown".

        The two are not the same string and should not be: one is what the file said,
        the other is that the file said nothing. What matters is that neither is
        answered with a block identifier, and that `carries_phase` refuses both, so
        neither can supply heterozygous phase evidence.
        """
        pool = {LOW: Zygosity.HET, MID_REF: Zygosity.HET, HIGH: Zygosity.HET}
        common = {LOW: {"GT": "0|1", "PS": "100"}, HIGH: {"GT": "0|1", "PS": "100"}}
        absent = self.run_add_phasing(pool, {**common, MID_REF: {"GT": "1|0"}})
        dotted = self.run_add_phasing(
            pool, {**common, MID_REF: {"GT": "1|0", "PS": "."}})
        self.assertEqual(absent.variant_pool_phase_set[MID_REF], "unknown")
        self.assertEqual(dotted.variant_pool_phase_set[MID_REF], ".")
        for bg in (absent, dotted):
            self.assertFalse(carries_phase(bg.variant_pool_phase_set[MID_REF]))

    def test_a_reference_token_keeps_the_set_its_own_row_reported(self):
        """A lane reference row copies the measured row's PS, which is not an inference.

        Its alternate half is often recorded as unused, so there is no partner in the
        pool to borrow from; before this the span of a nearby block answered instead,
        and happened to give the same number.
        """
        pool = {LOW: Zygosity.HET, MID_REF: Zygosity.HET, HIGH: Zygosity.HET}
        metrics = {LOW: {"GT": "0|1", "PS": "100"},
                   MID_REF: {"GT": "0|1", "PS": "250"},
                   HIGH: {"GT": "0|1", "PS": "100"}}
        bg = self.run_add_phasing(pool, metrics)
        self.assertEqual(bg.variant_pool_phase_set[MID_REF], "250")

    def test_non_reference_token_with_no_ps_also_declines(self):
        """The other route in: a non-reference token whose PS field the input omits.

        No reference partner shares its position here, so before the correction the
        span of block 100 answered for it.
        """
        pool = {LOW: Zygosity.HET, MID_ALT: Zygosity.HET, HIGH: Zygosity.HET}
        metrics = {LOW: {"GT": "0|1", "PS": "100"},
                   MID_ALT: {"GT": "1|0"},
                   HIGH: {"GT": "0|1", "PS": "100"}}
        bg = self.run_add_phasing(pool, metrics)
        self.assertEqual(bg.variant_pool_phase_set[MID_ALT], "unknown")

    def test_a_reported_set_on_a_non_reference_token_is_untouched(self):
        """A set the input actually reported stays exactly as reported."""
        pool = {LOW: Zygosity.HET, MID_ALT: Zygosity.HET}
        metrics = {LOW: {"GT": "0|1", "PS": "100"},
                   MID_ALT: {"GT": "1|0", "PS": "777"}}
        bg = self.run_add_phasing(pool, metrics)
        self.assertEqual(bg.variant_pool_phase_set[MID_ALT], "777")


if __name__ == "__main__":
    unittest.main()
