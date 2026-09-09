"""Regressions for phase evidence used to remove alleles and KN CO pairs."""

import unittest

from rbceq2.core_logic.alleles import Allele, BloodGroup, Pair
from rbceq2.core_logic.constants import AlleleState, SYNTHESISED_HOM_REF_GT
from rbceq2.core_logic.utils import Zygosity
from rbceq2.filters.knops import remove_unphased_co
from rbceq2.filters.phased import remove_unphased
from rbceq2.filters.shared_filter_functionality import identify_unphased


MISSING = object()


class TestAllelePhaseEvidence(unittest.TestCase):
    """Exercise actual removal decisions without depending on evidence helpers."""

    FIRST = "1:101_A_G"
    SECOND = "1:202_C_T"
    CONTROL = "1:303_T_C"

    @staticmethod
    def _allele(name, variants):
        """Build a generic TEST allele with no biological phenotype claim."""
        return Allele(
            genotype=name,
            phenotype=".",
            genotype_alt=".",
            phenotype_alt=".",
            defining_variants=frozenset(variants),
            null=False,
            weight_geno=1000,
            sub_type="TEST*01",
        )

    def _group(
        self,
        phases=("0|1", "1|0"),
        phase_sets=("100", "100"),
        zygosities=(Zygosity.HET, Zygosity.HET),
        second=None,
        extra=(),
    ):
        """Build independent allele and pair states over the same evidence.

        The KN type enables the CO consumer; the TEST alleles are generic
        fixtures rather than additions to the curated KN definitions.
        """
        second = second or self.SECOND
        tokens = [self.FIRST, second]
        observations = list(zip(tokens, phases, phase_sets, zygosities))
        observations.extend(extra)
        target = self._allele("TEST*01.01", [row[0] for row in observations])
        control = self._allele("TEST*01.02", [self.CONTROL])
        target_pair = Pair(target, control)
        control_pair = Pair(control, control)
        bg = BloodGroup(
            type="KN",
            sample="phase_evidence",
            alleles={
                AlleleState.FILT: [target, control],
                AlleleState.NORMAL: [],
                AlleleState.CO: [target_pair, control_pair],
            },
            variant_pool={self.CONTROL: Zygosity.HOM},
            variant_pool_phase={self.CONTROL: "1/1"},
            variant_pool_phase_set={self.CONTROL: "."},
        )
        for token, phase, phase_set, zygosity in observations:
            if zygosity is not MISSING:
                bg.variant_pool[token] = zygosity
            if phase is not MISSING:
                bg.variant_pool_phase[token] = phase
            if phase_set is not MISSING:
                bg.variant_pool_phase_set[token] = phase_set
        return bg, target, control, target_pair, control_pair

    def _assert_outcome(self, excluded, **observations):
        """Check the shared decision and each consumer with fresh real objects."""
        for consumer in ("identify", "normal", "co"):
            with self.subTest(consumer=consumer):
                bg, target, control, target_pair, control_pair = self._group(
                    **observations
                )
                original_evidence = (
                    dict(bg.variant_pool),
                    dict(bg.variant_pool_phase),
                    dict(bg.variant_pool_phase_set),
                )
                if consumer == "identify":
                    self.assertEqual(
                        identify_unphased(bg, [target, control]),
                        [target] if excluded else [],
                    )
                    self.assertEqual(dict(bg.filtered_out), {})
                elif consumer == "normal":
                    remove_unphased({"KN": bg}, phased=True)
                    self.assertEqual(
                        bg.alleles[AlleleState.FILT],
                        [control] if excluded else [target, control],
                    )
                    self.assertEqual(
                        dict(bg.filtered_out),
                        {"remove_unphased": [target]} if excluded else {},
                    )
                    self.assertEqual(
                        bg.alleles[AlleleState.CO], [target_pair, control_pair]
                    )
                else:
                    remove_unphased_co({"KN": bg}, phased=True)
                    self.assertEqual(
                        bg.alleles[AlleleState.CO],
                        [control_pair] if excluded else [target_pair, control_pair],
                    )
                    self.assertEqual(
                        dict(bg.filtered_out),
                        {"remove_unphased_co": [target_pair]} if excluded else {},
                    )
                    self.assertEqual(bg.alleles[AlleleState.FILT], [target, control])
                self.assertEqual(
                    (
                        bg.variant_pool,
                        bg.variant_pool_phase,
                        bg.variant_pool_phase_set,
                    ),
                    original_evidence,
                )

    def test_known_shared_block_cis_retains_allele_and_pairs(self):
        for phase_set in ("100", "0"):
            for phase in ("0|1", "1|0"):
                with self.subTest(phase_set=phase_set, phase=phase):
                    self._assert_outcome(
                        False,
                        phases=(phase, phase),
                        phase_sets=(phase_set, phase_set),
                    )

    def test_known_shared_block_trans_removes_with_existing_reasons(self):
        for phase_set in ("100", "0"):
            for phases in (("0|1", "1|0"), ("1|0", "0|1")):
                with self.subTest(phase_set=phase_set, phases=phases):
                    self._assert_outcome(
                        True, phases=phases, phase_sets=(phase_set, phase_set)
                    )

    def test_independent_block_flips_cannot_remove_an_allele(self):
        for first_phase in ("0|1", "1|0"):
            for second_phase in ("0|1", "1|0"):
                with self.subTest(first=first_phase, second=second_phase):
                    self._assert_outcome(
                        False,
                        phases=(first_phase, second_phase),
                        phase_sets=("100", "200"),
                    )

    def test_unavailable_phase_sets_cannot_prove_trans(self):
        for unavailable in (".", "unknown", "", None, MISSING):
            for phase_sets in (
                (unavailable, unavailable),
                ("100", unavailable),
                (unavailable, "100"),
            ):
                with self.subTest(phase_sets=phase_sets):
                    self._assert_outcome(False, phase_sets=phase_sets)

    def test_unusable_orientations_cannot_prove_trans_with_a_known_ps(self):
        # Include HOM and haploid GT spellings beside an HET pool entry: raw
        # strings alone must not manufacture a heterozygous orientation.
        for unavailable in (
            "unknown",
            ".",
            "",
            SYNTHESISED_HOM_REF_GT,
            "0/1",
            "1/0",
            "./.",
            "1/1",
            "1|1",
            "0|0",
            "1",
            "0",
            ".|1",
            "1|.",
            None,
            MISSING,
        ):
            with self.subTest(orientation=unavailable):
                self._assert_outcome(False, phases=("0|1", unavailable))

    def test_non_het_tokens_cannot_supply_a_conflicting_side(self):
        for zygosity in (
            Zygosity.HOM,
            Zygosity.HEM,
            Zygosity.NO_DATA,
            Zygosity.NO_COPIES,
            MISSING,
        ):
            for zygosities in (
                (Zygosity.HET, zygosity),
                (zygosity, Zygosity.HET),
            ):
                with self.subTest(zygosities=zygosities):
                    self._assert_outcome(False, zygosities=zygosities)

    def test_same_ps_label_on_different_chromosomes_is_not_a_shared_block(self):
        self._assert_outcome(False, second="2:202_C_T")

    def test_chromosome_prefix_does_not_hide_a_real_shared_block_conflict(self):
        self._assert_outcome(True, second="chr1:202_C_T")

    def test_unresolved_third_token_does_not_hide_a_known_trans_conflict(self):
        for phase_set in (".", "unknown", "200", MISSING):
            with self.subTest(third_phase_set=phase_set):
                self._assert_outcome(
                    True,
                    extra=(("1:404_G_T", "0|1", phase_set, Zygosity.HET),),
                )

    def test_consumers_keep_all_candidates_when_phasing_is_disabled(self):
        bg, target, control, target_pair, control_pair = self._group()
        self.assertEqual(identify_unphased(bg, [target, control]), [target])
        remove_unphased({"KN": bg}, phased=False)
        remove_unphased_co({"KN": bg}, phased=False)
        self.assertEqual(bg.alleles[AlleleState.FILT], [target, control])
        self.assertEqual(bg.alleles[AlleleState.CO], [target_pair, control_pair])
        self.assertEqual(dict(bg.filtered_out), {})


if __name__ == "__main__":
    unittest.main()
