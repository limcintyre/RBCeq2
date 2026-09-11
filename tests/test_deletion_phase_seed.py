"""The deletion phase repair may only seed from a token in the block it writes.

`modify_phase_of_large_indel` gives an unphased deletion an orientation copied from an
overlapping phased token, and records it against the one phase set the blood group's
pool declares. The token it copies from therefore has to be in that phase set, and has
to be heterozygous, or the orientation written is relative to something else while
being indistinguishable afterwards from one the input reported.
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


if __name__ == "__main__":
    unittest.main()
