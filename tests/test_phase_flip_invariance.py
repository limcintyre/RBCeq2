"""Reversing a phase block must not change any phase consumer's decision.

A phase set states how the heterozygous genotypes inside it are oriented relative to
each other. It says nothing about how that block sits against another block, so
reversing every genotype in one block - '0|1' to '1|0' - changes no claim the input
makes. A decision that moves when that is done was read from something the input
never said.

Each test builds a state spanning at least two flip units, runs one consumer on it,
then runs the consumer again on every combination of units reversed and requires the
same surviving alleles and pairs, the same single-slot genotypes and the same named
exclusions. A unit is a (chromosome, phase set) block. Barred genotypes whose phase
set is '.', 'unknown' or absent form one further unit per chromosome: they carry no
orientation evidence, so reversing them on their own must not matter either.
Reversing every unit together is one of the combinations, which covers the shared
block case.

Each test also checks that its consumer acted on the in-block evidence, so a fixture
on which nothing happens cannot pass as invariance.

Consumers, in pipeline order: remove_unphased, filter_pairs_by_phase,
narrow_second_slot_candidates_by_phase, cant_name_second_slot_cuz_hom_ref_impossible
after cant_be_hom_ref_due_to_HET_SNP, cant_name_second_slot_cuz_ref_not_phased after
remove_unphased and ref_not_phased, remove_unphased_co, impossible_alleles_phased,
filter_if_all_HET_vars_on_same_side_and_phased, the three relationship filters,
rm_ref_if_2x_HET_phased and low_weight_hom.
"""

import csv
import unittest
from collections import defaultdict
from itertools import combinations
from pathlib import Path
from unittest.mock import patch

from rbceq2.core_logic.alleles import Allele, BloodGroup, Pair
from rbceq2.core_logic.constants import UNDETERMINED_SLOT, AlleleState
from rbceq2.core_logic.utils import Zygosity
from rbceq2.filters.knops import remove_unphased_co
from rbceq2.filters.phased import (
    cant_be_hom_ref_due_to_HET_SNP,
    cant_name_second_slot_cuz_hom_ref_impossible,
    cant_name_second_slot_cuz_ref_not_phased,
    filter_if_all_HET_vars_on_same_side_and_phased,
    filter_on_in_relationship_if_all_HOM_and_phased,
    filter_on_in_relationship_if_HET_vars_on_dif_side_and_phased,
    filter_on_in_relationship_when_HOM_cant_be_on_one_side,
    filter_pairs_by_phase,
    impossible_alleles_phased,
    low_weight_hom,
    narrow_second_slot_candidates_by_phase,
    ref_not_phased,
    remove_unphased,
    rm_ref_if_2x_HET_phased,
)

DB = Path(__file__).resolve().parents[1] / "src/rbceq2/resources/db.tsv"
NO_SET = "<no phase set>"
HET, HOM = Zygosity.HET, Zygosity.HOM


def curated(names: set[str]) -> dict[str, Allele]:
    """Read allele definitions from the curated database, GRCh38 tokens."""
    found: dict[str, Allele] = {}
    with DB.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            name = row["Genotype"]
            if name not in names:
                continue
            if name in found:
                raise ValueError(f"Duplicate curated fixture allele: {name}")
            chrom = row["Chrom"].removeprefix("chr")
            found[name] = Allele(
                genotype=name,
                genotype_alt=row["Genotype_alt"],
                defining_variants=frozenset(
                    f"{chrom}:{token.strip()}" for token in row["GRCh38"].split(",")
                ),
                phenotype=row["Phenotype_change"] or row["Phenotype"],
                phenotype_alt=row["Phenotype_alt_change"] or row["Phenotype_alt"],
                null=False,
                weight_geno=int(row["Weight_of_genotype"] or 1000),
                reference=row["Reference_genotype"] == "Yes",
                sub_type=row["Sub_type"],
            )
    missing = names - set(found)
    if missing:
        raise ValueError(f"Missing curated fixture alleles: {sorted(missing)}")
    return found


def synthetic(name, variants, *, weight=1000, reference=False, sub_type="TEST*01"):
    """A definition built for the test, with no biological claim beyond its tokens."""
    return Allele(
        genotype=name,
        phenotype=".",
        genotype_alt=name,
        phenotype_alt=".",
        defining_variants=frozenset(variants),
        null=False,
        weight_geno=weight,
        reference=reference,
        sub_type=sub_type,
    )


def group(bg_type, sample, alleles, observations, filtered_out=None):
    """A BloodGroup whose pools hold one (zygosity, GT, phase set) per token.

    A phase set of None leaves the token out of the phase set pool.
    """
    return BloodGroup(
        type=bg_type,
        sample=sample,
        alleles=alleles,
        variant_pool={token: obs[0] for token, obs in observations.items()},
        variant_pool_phase={token: obs[1] for token, obs in observations.items()},
        variant_pool_phase_set={
            token: obs[2] for token, obs in observations.items() if obs[2] is not None
        },
        filtered_out=defaultdict(list) if filtered_out is None else filtered_out,
    )


def flip_unit(bg: BloodGroup, variant: str):
    """The unit a token's orientation belongs to, or None if it has no orientation."""
    gt = bg.variant_pool_phase.get(variant)
    if not isinstance(gt, str) or "/" in gt or "|" not in gt:
        return None
    sides = gt.split("|")
    if len(set(sides)) == 1:
        return None  # '1|1' and '0|0' are on both haplotypes; reversing is a no-op
    chrom = variant.partition(":")[0].removeprefix("chr")
    phase_set = bg.variant_pool_phase_set.get(variant)
    if phase_set in (None, "", ".", "unknown"):
        phase_set = NO_SET
    return chrom, phase_set


def reverse(bg: BloodGroup, units: set) -> None:
    """Reverse the haplotype order of every genotype in the chosen units."""
    for variant, gt in list(bg.variant_pool_phase.items()):
        if flip_unit(bg, variant) in units:
            bg.variant_pool_phase[variant] = "|".join(reversed(gt.split("|")))


def label(item) -> str:
    if isinstance(item, Pair):
        return f"{item.allele1.genotype}/{item.allele2.genotype}"
    if isinstance(item, Allele):
        return item.genotype
    return repr(item)


def decision(bg: BloodGroup) -> dict:
    """Everything a phase consumer decides: survivors, single slots, named trail."""
    return {
        "alleles": {
            str(state): None if items is None else [label(item) for item in items]
            for state, items in bg.alleles.items()
        },
        "filtered_out": {
            reason: [label(item) for item in items]
            for reason, items in bg.filtered_out.items()
            if items
        },
        "single_slot_genotypes": list(bg.single_slot_genotypes),
    }


def run(*stages):
    """Apply stages in order, each as (function, keyword arguments)."""

    def apply(bg: BloodGroup) -> None:
        for stage, kwargs in stages:
            stage({bg.type: bg}, **kwargs)

    return apply


PHASED = {"phased": True}


class _FlipInvariance(unittest.TestCase):
    def assert_flip_invariant(self, build, apply) -> dict:
        """Decide on the state as built and on every combination of reversed units.

        Returns the decision, for the caller to check that the consumer acted.
        """
        reference = build()
        units = sorted(
            {
                unit
                for unit in (flip_unit(reference, v) for v in reference.variant_pool_phase)
                if unit is not None
            }
        )
        self.assertGreaterEqual(len(units), 2, "one unit cannot show independence")
        apply(reference)
        expected = decision(reference)
        for size in range(1, len(units) + 1):
            for chosen in combinations(units, size):
                with self.subTest(reversed_units=chosen):
                    bg = build()
                    reverse(bg, set(chosen))
                    apply(bg)
                    self.assertEqual(decision(bg), expected)
        return expected


class TestAlleleRemoval(_FlipInvariance):
    """remove_unphased and remove_unphased_co: opposite sides inside one block only."""

    FIRST, SECOND = "1:101_A_G", "1:102_C_T"  # block 100, opposite sides
    OTHER_BLOCK = "1:201_G_A"  # block 200
    NO_SET_TOKEN = "1:301_T_C"  # barred, no phase set

    def build(self):
        split = synthetic("TEST*01.01", (self.FIRST, self.SECOND))
        spread = synthetic("TEST*01.02", (self.FIRST, self.OTHER_BLOCK, self.NO_SET_TOKEN))
        straddling = synthetic("TEST*01.03", (self.SECOND, self.OTHER_BLOCK))
        # KN, so that the CO consumer has a state to work on.
        return group(
            "KN",
            "flip_allele_removal",
            {
                AlleleState.FILT: [split, spread, straddling],
                AlleleState.NORMAL: [],
                AlleleState.CO: [Pair(split, spread), Pair(spread, straddling)],
            },
            {
                self.FIRST: (HET, "0|1", "100"),
                self.SECOND: (HET, "1|0", "100"),
                self.OTHER_BLOCK: (HET, "0|1", "200"),
                self.NO_SET_TOKEN: (HET, "1|0", "."),
            },
        )

    def test_remove_unphased(self):
        expected = self.assert_flip_invariant(self.build, run((remove_unphased, PHASED)))
        self.assertEqual(expected["filtered_out"], {"remove_unphased": ["TEST*01.01"]})

    def test_remove_unphased_co(self):
        expected = self.assert_flip_invariant(self.build, run((remove_unphased_co, PHASED)))
        removed = self.build().alleles[AlleleState.CO][0]
        self.assertEqual(expected["filtered_out"], {"remove_unphased_co": [label(removed)]})


class TestPairFilters(_FlipInvariance):
    """Pair stages that compare two alleles' sides: only inside one shared block."""

    LEFT = "1:100_A_G"  # '0|1', block 1000
    RIGHT = "1:200_C_T"  # '1|0', block 1000: opposite LEFT
    SAME = "1:250_G_A"  # '0|1', block 1000: beside LEFT
    OTHER_BLOCK = "1:400_T_C"  # block 2000
    NO_SET_TOKEN = "1:500_G_T"  # barred, no phase set
    REFERENCE = synthetic("TEST*01", (), reference=True)

    def build(self, bg_type="TEST", co=False):
        left = synthetic("TEST*01.01", (self.LEFT,), weight=1)
        right = synthetic("TEST*01.02", (self.RIGHT,), weight=2)
        same = synthetic("TEST*01.03", (self.SAME,), weight=4)
        other = synthetic("TEST*01.04", (self.OTHER_BLOCK,), weight=8)
        unset = synthetic("TEST*01.05", (self.NO_SET_TOKEN,), weight=16)
        ref = self.REFERENCE
        pairs = [
            Pair(left, right), Pair(right, same), Pair(left, same),
            Pair(left, other), Pair(right, other), Pair(left, unset),
            Pair(ref, left), Pair(ref, right), Pair(ref, other),
        ]
        return group(
            bg_type,
            "flip_pair_filters",
            {AlleleState.NORMAL: list(pairs), AlleleState.CO: list(pairs) if co else None},
            {
                self.LEFT: (HET, "0|1", "1000"),
                self.RIGHT: (HET, "1|0", "1000"),
                self.SAME: (HET, "0|1", "1000"),
                self.OTHER_BLOCK: (HET, "0|1", "2000"),
                self.NO_SET_TOKEN: (HET, "1|0", None),
            },
        )

    def pair(self, first, second):
        return label(next(
            p for p in self.build().alleles[AlleleState.NORMAL]
            if {p.allele1.genotype, p.allele2.genotype} == {first, second}
        ))

    def test_filter_pairs_by_phase(self):
        expected = self.assert_flip_invariant(
            self.build,
            run((filter_pairs_by_phase,
                 {"phased": True, "reference_alleles": {"TEST": self.REFERENCE}})),
        )
        self.assertEqual(
            expected["filtered_out"],
            {"filter_pairs_by_phase": [self.pair("TEST*01.01", "TEST*01.03")]},
        )

    def test_same_side_in_both_states(self):
        reason = "filter_if_all_HET_vars_on_same_side_and_phased"
        expected = self.assert_flip_invariant(
            lambda: self.build("KN", co=True),
            run((filter_if_all_HET_vars_on_same_side_and_phased, PHASED)),
        )
        self.assertEqual(set(expected["filtered_out"]), {reason})
        self.assertEqual(
            set(expected["filtered_out"][reason]), {self.pair("TEST*01.01", "TEST*01.03")}
        )

    def test_rm_ref_if_2x_HET_phased(self):
        expected = self.assert_flip_invariant(
            self.build, run((rm_ref_if_2x_HET_phased, PHASED))
        )
        self.assertEqual(
            set(expected["filtered_out"]["rm_ref_if_2x_HET_phased"]),
            {self.pair("TEST*01", name) for name in ("TEST*01.01", "TEST*01.02", "TEST*01.04")},
        )

    def test_low_weight_hom(self):
        expected = self.assert_flip_invariant(self.build, run((low_weight_hom, PHASED)))
        self.assertEqual(
            expected["alleles"][str(AlleleState.NORMAL)],
            [self.pair("TEST*01.01", "TEST*01.02")],
        )


class TestRelationshipFilters(_FlipInvariance):
    """LU: two children trans in block 1000, a third child alone in block 2000."""

    HOM_REF = "19:44812188_ref"
    CHILD_04 = "19:44813269_G_A"
    CHILD_05 = "19:44812284_G_A"
    CHILD_19 = "19:44819487_A_G"
    NAMES = {"LU*02", "LU*02.-04.1", "LU*02.-05", "LU*02.19"}

    @classmethod
    def setUpClass(cls):
        cls.lu = curated(cls.NAMES)

    def build(self, hom_pair=False):
        ref, c04, c05, c19 = (self.lu[n] for n in ("LU*02", "LU*02.-04.1", "LU*02.-05", "LU*02.19"))
        pairs = [
            Pair(ref, c04), Pair(ref, c05), Pair(ref, c19),
            Pair(c04, c05), Pair(c04, c19), Pair(c05, c19),
        ]
        if hom_pair:
            pairs.insert(0, Pair(ref, ref))
        return group(
            "LU",
            "flip_relationships",
            {AlleleState.RAW: [ref, c04, c05, c19], AlleleState.NORMAL: pairs,
             AlleleState.CO: None},
            {
                self.HOM_REF: (HOM, "1/1", "."),
                self.CHILD_04: (HET, "1|0", "1000"),
                self.CHILD_05: (HET, "0|1", "1000"),
                self.CHILD_19: (HET, "0|1", "2000"),
            },
        )

    def in_block_parent_pairs(self):
        return {label(Pair(self.lu["LU*02"], self.lu[n])) for n in ("LU*02.-04.1", "LU*02.-05")}

    def test_het_vars_on_different_sides(self):
        reason = "filter_on_in_relationship_if_HET_vars_on_dif_side_and_phased"
        expected = self.assert_flip_invariant(
            self.build,
            run((filter_on_in_relationship_if_HET_vars_on_dif_side_and_phased, PHASED)),
        )
        self.assertLessEqual(self.in_block_parent_pairs(), set(expected["filtered_out"][reason]))

    def test_hom_cannot_be_on_one_side(self):
        reason = "filter_on_in_relationship_when_HOM_cant_be_on_one_side"
        expected = self.assert_flip_invariant(
            self.build,
            run((filter_on_in_relationship_when_HOM_cant_be_on_one_side, PHASED)),
        )
        self.assertLessEqual(self.in_block_parent_pairs(), set(expected["filtered_out"][reason]))

    def test_all_hom_and_phased(self):
        reason = "filter_on_in_relationship_if_all_HOM_and_phased"
        expected = self.assert_flip_invariant(
            lambda: self.build(hom_pair=True),
            run((filter_on_in_relationship_if_all_HOM_and_phased, PHASED)),
        )
        self.assertEqual(
            expected["filtered_out"], {reason: [label(Pair(self.lu["LU*02"], self.lu["LU*02"]))]}
        )


class TestImpossibleAlleles(_FlipInvariance):
    """KN: a reference containment in block 1000 and a genuine subset in block 2000."""

    BACKGROUND = "1:207609571_A_T"
    REF_TOKEN = "1:207609424_ref"
    OTHER_TOKEN = "1:207609424_G_A"
    CHILD_TOKEN = "1:207587428_C_T"
    SUPERSET_ONLY = "1:207609544_A_G"
    SHARED_SUBSET = "1:207609586_A_G"
    NAMES = {"KN*01", "KN*01.-05", "KN*02", "KN*01.10", "KN*01.07"}

    @classmethod
    def setUpClass(cls):
        cls.kn = curated(cls.NAMES)

    def pairs(self):
        other = self.kn["KN*02"]
        return [Pair(self.kn[n], other) for n in ("KN*01", "KN*01.-05", "KN*01.10", "KN*01.07")]

    def build(self, co=False):
        return group(
            "KN",
            "flip_impossible_alleles",
            {AlleleState.NORMAL: self.pairs(), AlleleState.CO: self.pairs() if co else None},
            {
                self.BACKGROUND: (HOM, "1/1", "."),
                self.REF_TOKEN: (HET, "0|1", "1000"),
                self.OTHER_TOKEN: (HET, "1|0", "1000"),
                self.CHILD_TOKEN: (HET, "0|1", "1000"),
                self.SUPERSET_ONLY: (HET, "0|1", "2000"),
                self.SHARED_SUBSET: (HET, "0|1", "2000"),
            },
        )

    def check(self, co):
        expected = self.assert_flip_invariant(
            lambda: self.build(co=co), run((impossible_alleles_phased, PHASED))
        )
        removed = {label(p) for p in self.pairs()[0::2]}  # KN*01 and KN*01.10
        self.assertEqual(
            set(expected["filtered_out"]["filter_impossible_alleles_phased"]), removed
        )

    def test_normal_state(self):
        self.check(co=False)

    def test_co_state(self):
        self.check(co=True)


class TestSingleSlotNarrowing(_FlipInvariance):
    """Candidates for one named slot, narrowed only by forced in-block extensions."""

    RHCE = {
        "RHCE*01.01", "RHCE*01.02.01", "RHCE*01.20.02.01",
        "RHCE*01.20.02.02", "RHCE*01.20.04.01", "RHCE*01.20.04.02",
    }
    RHCE_INDEPENDENT = "1:25385759_G_A"
    RHCE_LINKED = {"1:25390806_A_G", "1:25390817_G_C", "1:25408711_ref", "1:25420682_G_A"}
    RHCE_HOM = {"1:25390874_ref", "1:25420739_ref"}
    REASON = "narrow_second_slot_candidates_by_phase"

    @classmethod
    def setUpClass(cls):
        cls.rhce = curated(cls.RHCE)

    def build_rhce(self):
        alleles = [self.rhce[name] for name in sorted(self.RHCE)]
        observations = {token: (HOM, "1/1", ".") for token in self.RHCE_HOM}
        observations.update({token: (HET, "0|1", "25233074") for token in self.RHCE_LINKED})
        observations[self.RHCE_INDEPENDENT] = (HET, "0|1", "25385759")
        bg = group(
            "RHCE", "flip_rhce_narrowing",
            {AlleleState.RAW: alleles, AlleleState.NORMAL: [], AlleleState.CO: None},
            observations,
        )
        bg.single_slot_genotypes = [f"{a.genotype}/{UNDETERMINED_SLOT}" for a in alleles]
        return bg

    def test_rhce_independent_block_keeps_two_candidates(self):
        expected = self.assert_flip_invariant(
            self.build_rhce, run((narrow_second_slot_candidates_by_phase, PHASED))
        )
        self.assertEqual(
            set(expected["single_slot_genotypes"]),
            {f"{n}/{UNDETERMINED_SLOT}" for n in ("RHCE*01.20.02.02", "RHCE*01.20.04.02")},
        )

    ANCHOR, OTHER_BLOCK, EXTRA = "1:100_A_G", "1:200_C_T", "1:150_G_T"

    def build_forced(self):
        subset = synthetic("TEST*01.01", (self.ANCHOR, self.OTHER_BLOCK))
        superset = synthetic("TEST*01.02", (self.ANCHOR, self.OTHER_BLOCK, self.EXTRA))
        bg = group(
            "TEST", "flip_forced_extension",
            {AlleleState.RAW: [subset, superset], AlleleState.NORMAL: []},
            {
                self.ANCHOR: (HET, "0|1", "1000"),
                self.OTHER_BLOCK: (HET, "0|1", "2000"),
                self.EXTRA: (HET, "0|1", "1000"),
            },
        )
        bg.single_slot_genotypes = [f"{a.genotype}/{UNDETERMINED_SLOT}" for a in (subset, superset)]
        return bg

    def test_forced_in_block_extension_prunes_the_subset(self):
        expected = self.assert_flip_invariant(
            self.build_forced, run((narrow_second_slot_candidates_by_phase, PHASED))
        )
        self.assertEqual(expected["single_slot_genotypes"], [f"TEST*01.02/{UNDETERMINED_SLOT}"])
        self.assertEqual(expected["filtered_out"], {self.REASON: ["TEST*01.01"]})


class TestSecondSlotRescues(_FlipInvariance):
    """Naming one slot needs every HET definition of that allele placed in one block."""

    REF_874, ALT_874 = "1:25390874_ref", "1:25390874_C_G"
    REF_711, ALT_711 = "1:25408711_ref", "1:25408711_G_A"
    REF_739, ALT_739 = "1:25420739_ref", "1:25420739_G_C"
    UNRELATED = "1:25385759_G_A"  # in the pool, in no definition, in its own block

    def reference(self):
        return synthetic("RHCE*01", (self.REF_874, self.REF_711, self.ALT_739),
                         reference=True, sub_type="RHCE*01")

    def build_hom_ref(self, alt_739_block):
        ref = self.reference()
        return group(
            "RHCE", "flip_hom_ref_rescue",
            {AlleleState.NORMAL: [Pair(ref, ref)], AlleleState.CO: None},
            {
                self.REF_874: (HOM, "1/1", "."),
                self.ALT_711: (HET, "0|1", "25211850"),
                self.REF_711: (HET, "1|0", "25211850"),
                self.ALT_739: (HET, "1|0", alt_739_block),
                self.UNRELATED: (HET, "0|1", "25385759"),
            },
        )

    def hom_ref_chain(self):
        return run((cant_be_hom_ref_due_to_HET_SNP, PHASED),
                   (cant_name_second_slot_cuz_hom_ref_impossible, PHASED))

    def test_reference_named_when_its_het_definitions_share_a_block(self):
        expected = self.assert_flip_invariant(
            lambda: self.build_hom_ref("25211850"), self.hom_ref_chain()
        )
        self.assertEqual(expected["single_slot_genotypes"], [f"RHCE*01/{UNDETERMINED_SLOT}"])

    def test_reference_not_named_across_independent_blocks(self):
        expected = self.assert_flip_invariant(
            lambda: self.build_hom_ref("2000"), self.hom_ref_chain()
        )
        self.assertEqual(expected["single_slot_genotypes"], [])
        self.assertIn("cant_be_hom_ref_due_to_HET_SNP", expected["filtered_out"])

    def build_ref_not_phased(self, block_739):
        ref = self.reference()
        partner = synthetic("RHCE*03", (self.ALT_874, self.REF_711, self.ALT_739),
                            sub_type="RHCE*01")
        return group(
            "RHCE", "flip_partner_rescue",
            # NORMAL holds the pair process_genetic_data re-adds the reference into.
            {AlleleState.FILT: [ref, partner], AlleleState.NORMAL: [Pair(ref, partner)],
             AlleleState.CO: None},
            {
                self.REF_874: (HET, "0|1", "25214110"),
                self.ALT_874: (HET, "1|0", "25214110"),
                self.REF_711: (HET, "1|0", "25214110"),
                self.ALT_711: (HET, "0|1", "25214110"),
                self.REF_739: (HET, "0|1", block_739),
                self.ALT_739: (HET, "1|0", block_739),
                self.UNRELATED: (HET, "0|1", "25385759"),
            },
        )

    def partner_chain(self):
        """The three stages that reach the partner rescue, in pipeline order.

        ref_not_phased empties the blood group here, which is the shape the rescue
        exists for, so its 'all pairs removed' warning is expected rather than
        informative and is kept out of the suite's output.
        """

        def apply(bg):
            with patch("rbceq2.core_logic.alleles.logger"):
                run((remove_unphased, PHASED), (ref_not_phased, PHASED),
                    (cant_name_second_slot_cuz_ref_not_phased, PHASED))(bg)

        return apply

    def test_partner_named_when_its_het_definitions_share_a_block(self):
        expected = self.assert_flip_invariant(
            lambda: self.build_ref_not_phased("25214110"), self.partner_chain()
        )
        self.assertEqual(expected["single_slot_genotypes"], [f"RHCE*03/{UNDETERMINED_SLOT}"])
        self.assertEqual(expected["filtered_out"]["remove_unphased"], ["RHCE*01"])

    def test_partner_not_named_across_independent_blocks(self):
        expected = self.assert_flip_invariant(
            lambda: self.build_ref_not_phased("2000"), self.partner_chain()
        )
        self.assertEqual(expected["single_slot_genotypes"], [])
        self.assertEqual(expected["filtered_out"]["remove_unphased"], ["RHCE*01"])


if __name__ == "__main__":
    unittest.main()
