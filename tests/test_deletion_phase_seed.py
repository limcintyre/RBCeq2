"""The deletion phase repair may only seed from, and write to, the block it declares.

`modify_phase_of_large_indel` gives an unphased deletion an orientation copied from an
overlapping phased token, and records it against the one phase set the blood group's
pool declares. The token it copies from therefore has to be in that phase set, and has
to be heterozygous, or the orientation written is relative to something else while
being indistinguishable afterwards from one the input reported.

The same applies to what it writes. A reference token is oriented here only as the
opposite side of a deletion it overlaps, so only a deletion in that same phase set can
place it, and a reference token overlapping no deletion is not this routine's to
relabel at all.
"""

import unittest

from rbceq2.core_logic.alleles import BloodGroup
from rbceq2.core_logic.data_procesing import modify_phase_of_large_indel
from rbceq2.core_logic.utils import Zygosity

DELETION = "19:44812161_del_27kb"   # LU*02N.06, spans 44812161-44839161
INSIDE = "19:44812284_G_A"          # LU*02.-05, inside the deletion
OUTSIDE = "19:44811263_G_A"         # LU*02.-28, outside it, declares block 100


class TestDeletionPhaseSeed(unittest.TestCase):
    @staticmethod
    def build(seed_phase, seed_ps, seed_zygosity=Zygosity.HET, del_phase="0/1"):
        """One unphased deletion, one candidate seed inside it, one token in block 100."""
        return BloodGroup(
            type="LU", sample="seed", alleles={},
            variant_pool={DELETION: Zygosity.HET, INSIDE: seed_zygosity,
                          OUTSIDE: Zygosity.HET},
            variant_pool_phase={DELETION: del_phase, INSIDE: seed_phase,
                                OUTSIDE: "0|1"},
            variant_pool_phase_set={DELETION: ".", INSIDE: seed_ps,
                                    OUTSIDE: "100"},
        )

    def repair(self, bg):
        modify_phase_of_large_indel({"LU": bg}, phased=True)
        return bg.variant_pool_phase[DELETION], bg.variant_pool_phase_set[DELETION]

    def test_seed_in_the_written_block_still_places_the_deletion(self):
        """The case the routine exists for is unchanged."""
        phase, phase_set = self.repair(self.build("0|1", "100"))
        self.assertEqual(phase, "1|0")
        self.assertEqual(phase_set, "100")

    def test_seed_naming_no_phase_set_is_refused(self):
        """A bar whose own phase set is "." is relative to nothing the file declares."""
        phase, phase_set = self.repair(self.build("0|1", "."))
        self.assertEqual(phase, "0/1")
        self.assertEqual(phase_set, ".")

    def test_seed_with_an_unknown_phase_set_is_refused(self):
        """"unknown" is excluded from the phase set gate, so it cannot seed either."""
        phase, phase_set = self.repair(self.build("0|1", "unknown"))
        self.assertEqual(phase, "0/1")
        self.assertEqual(phase_set, ".")

    def test_flipping_a_seed_outside_the_block_no_longer_moves_the_deletion(self):
        """The defect: a token in no block decided which side the deletion took."""
        forward = self.repair(self.build("0|1", "."))
        reverse = self.repair(self.build("1|0", "."))
        self.assertEqual(forward, reverse)

    def test_a_homozygous_seed_cannot_place_the_deletion(self):
        """A token on both copies does not locate anything on one of them."""
        phase, phase_set = self.repair(
            self.build("1|1", "100", seed_zygosity=Zygosity.HOM))
        self.assertEqual(phase, "0/1")
        self.assertEqual(phase_set, ".")


INSIDE_REF = "19:44812188_ref"       # the LU lane locus, inside the deletion
OUTSIDE_REF = "19:44811264_ref"      # before the deletion starts, overlapping none
FAR_DELETION = "19:44900000_del_1kb"  # a second deletion, well clear of the first
FAR_REF = "19:44900500_ref"          # inside that second deletion


class TestDeletionReferenceAssignment(unittest.TestCase):
    """What the repair may write on a reference token, and on which."""

    @staticmethod
    def build(outside_ref_phase="0|1", outside_ref_ps=".", far=False):
        """One repairable deletion in block 100, with reference tokens around it.

        INSIDE_REF sits inside that deletion, OUTSIDE_REF before it overlaps nothing.
        With far=True a second deletion the caller phased itself, naming no set, holds
        a reference token of its own.
        """
        pool = {DELETION: Zygosity.HET, INSIDE: Zygosity.HET,
                INSIDE_REF: Zygosity.HET, OUTSIDE_REF: Zygosity.HET}
        phase = {DELETION: "0/1", INSIDE: "0|1", INSIDE_REF: "1|0",
                 OUTSIDE_REF: outside_ref_phase}
        phase_sets = {DELETION: ".", INSIDE: "100", INSIDE_REF: ".",
                      OUTSIDE_REF: outside_ref_ps}
        if far:
            pool[FAR_DELETION] = Zygosity.HET
            pool[FAR_REF] = Zygosity.HET
            phase[FAR_DELETION] = "1|0"
            phase[FAR_REF] = "1|0"
            phase_sets[FAR_DELETION] = "."
            phase_sets[FAR_REF] = "."
        return BloodGroup(
            type="LU", sample="reference_assignment", alleles={},
            variant_pool=pool, variant_pool_phase=phase,
            variant_pool_phase_set=phase_sets,
        )

    @staticmethod
    def repair(bg):
        modify_phase_of_large_indel({"LU": bg}, phased=True)
        return bg

    def test_a_reference_inside_the_placed_deletion_takes_the_other_side(self):
        """The case the reference assignment exists for is unchanged."""
        bg = self.repair(self.build())
        self.assertEqual(bg.variant_pool_phase[DELETION], "1|0")
        self.assertEqual(bg.variant_pool_phase[INSIDE_REF], "0|1")
        self.assertEqual(bg.variant_pool_phase_set[INSIDE_REF], "100")

    def test_a_reference_overlapping_no_deletion_keeps_its_own_phase_set(self):
        """The defect: this token was written into a block nothing linked it to."""
        for phase_set in (".", "unknown"):
            with self.subTest(phase_set=phase_set):
                bg = self.repair(self.build(outside_ref_ps=phase_set))
                self.assertEqual(bg.variant_pool_phase_set[OUTSIDE_REF], phase_set)
                self.assertEqual(bg.variant_pool_phase[OUTSIDE_REF], "0|1")

    def test_a_reference_already_in_the_written_block_is_left_where_it_is(self):
        """Nothing to correct where the file put the token in that block itself."""
        bg = self.repair(self.build(outside_ref_ps="100"))
        self.assertEqual(bg.variant_pool_phase_set[OUTSIDE_REF], "100")
        self.assertEqual(bg.variant_pool_phase[OUTSIDE_REF], "0|1")

    def test_flipping_a_reference_that_names_no_set_changes_nothing(self):
        """Which side of an undeclared bar a token sits on is not evidence.

        The shape found in a real three row input: both orientations are the same
        claim, and the stamped token used to carry them into the only block in the
        pool, where a phase filter spent them on opposite candidates.
        """
        forward = self.repair(self.build(outside_ref_phase="0|1"))
        reverse = self.repair(self.build(outside_ref_phase="1|0"))
        for bg in (forward, reverse):
            self.assertEqual(bg.variant_pool_phase_set[OUTSIDE_REF], ".")
        self.assertEqual(
            {token: forward.variant_pool_phase_set[token]
             for token in forward.variant_pool_phase_set},
            {token: reverse.variant_pool_phase_set[token]
             for token in reverse.variant_pool_phase_set},
        )

    def test_a_deletion_naming_no_set_orients_no_reference(self):
        """Its own bar is relative to nothing declared, so it places nothing."""
        bg = self.repair(self.build(far=True))
        self.assertEqual(bg.variant_pool_phase[FAR_REF], "1|0")
        self.assertEqual(bg.variant_pool_phase_set[FAR_REF], ".")
        # the repairable deletion in block 100 is placed exactly as before
        self.assertEqual(bg.variant_pool_phase[DELETION], "1|0")
        self.assertEqual(bg.variant_pool_phase[INSIDE_REF], "0|1")


if __name__ == "__main__":
    unittest.main()
