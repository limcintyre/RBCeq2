"""Regressions for observed allele coverage in ABO and FUT3 candidate pairs."""

import csv
from itertools import combinations
from pathlib import Path
import unittest

from rbceq2.core_logic.alleles import Allele, BloodGroup, Pair
from rbceq2.core_logic.constants import AlleleState
from rbceq2.core_logic.utils import Zygosity
from rbceq2.filters.geno import ensure_HET_SNP_used
from rbceq2.filters.phased import (
    filter_on_in_relationship_when_HOM_cant_be_on_one_side,
)


class _CoverageFixtures(unittest.TestCase):
    """Load complete definitions from the snapshot database, without altering them."""

    @classmethod
    def setUpClass(cls):
        wanted = {
            "ABO*O.01.02", "ABO*O.01.13", "ABO*O.01.39", "ABO*O.01.83",
            "FUT3*01.17.01", "FUT3*01N.17.03", "FUT3*01N.17.05", "FUT3*01N.17.19",
            "FUT3*01N.01.01", "FUT3*01N.01.02", "FUT3*01N.01.08",
            "FUT3*01N.01.10", "FUT3*01N.01.11", "FUT3*01N.01.12",
        }
        cls.curated = {}
        db = Path(__file__).resolve().parents[1] / "src/rbceq2/resources/db.tsv"
        with db.open(newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                name = row["Genotype"]
                if name not in wanted:
                    continue
                if name in cls.curated:
                    raise ValueError(f"Duplicate curated coverage fixture: {name}")
                chrom = row["Chrom"].removeprefix("chr")
                cls.curated[name] = Allele(
                    genotype=name, genotype_alt=row["Genotype_alt"],
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
        if set(cls.curated) != wanted:
            raise ValueError("Missing curated ABO/FUT3 coverage fixture")

    def _pair(self, first, second):
        return Pair(self.curated[first], self.curated[second])

    def _assert_unchanged(self, bg, stage, **kwargs):
        before = list(bg.alleles[AlleleState.NORMAL])

        stage({bg.type: bg}, **kwargs)

        self.assertEqual(bg.alleles[AlleleState.NORMAL], before)
        self.assertEqual(dict(bg.filtered_out), {})


class TestABOHeterozygousCoverage(_CoverageFixtures):
    """NA18520 has one HET token and the full O.01.02 HOM background."""

    HET = "9:133259833_G_A"
    ADDITIONAL_HOM = "9:133256085_A_T"
    REASON = "ensure_HET_SNP_used"

    def _group(self):
        names = ("ABO*O.01.02", "ABO*O.01.13", "ABO*O.01.39", "ABO*O.01.83")
        alleles = [self.curated[name] for name in names]
        pool = {variant: Zygosity.HOM for variant in alleles[0].defining_variants}
        pool[self.HET] = Zygosity.HET
        return BloodGroup(
            type="ABO", sample="NA18520_coverage",
            alleles={AlleleState.RAW: alleles,
                     AlleleState.NORMAL: [Pair(first, second) for first, second in combinations(alleles, 2)],
                     AlleleState.CO: None},
            variant_pool=pool,
        )

    def test_one_het_with_hom_extensions_keeps_only_o0102_counterparts(self):
        bg = self._group()
        before = list(bg.alleles[AlleleState.NORMAL])
        self.assertEqual(bg.variant_pool_phase, {})
        self.assertEqual(bg.variant_pool_phase_set, {})
        self.assertEqual(sum(z == Zygosity.HET for z in bg.variant_pool.values()), 1)
        replacement = self.curated["ABO*O.01.02"]
        for smaller in bg.alleles[AlleleState.RAW][1:]:
            # These are multi-token extensions; the old exact +one-HET check is insufficient.
            self.assertGreater(len(replacement.defining_variants - smaller.defining_variants), 1)

        ensure_HET_SNP_used({"ABO": bg})

        self.assertEqual(bg.alleles[AlleleState.NORMAL], before[:3])
        self.assertEqual(set(bg.filtered_out), {self.REASON})
        self.assertCountEqual(bg.filtered_out[self.REASON], before[3:])

    def test_additional_definition_must_be_observed_hom(self):
        for zygosity in (Zygosity.HEM, Zygosity.NO_DATA, None):
            with self.subTest(zygosity=zygosity):
                bg = self._group()
                first, second = self.curated["ABO*O.01.39"], self.curated["ABO*O.01.83"]
                replacement = self.curated["ABO*O.01.02"]
                self.assertNotIn(self.ADDITIONAL_HOM, first.defining_variants | second.defining_variants)
                self.assertIn(self.ADDITIONAL_HOM, replacement.defining_variants)
                bg.alleles[AlleleState.RAW] = [first, second, replacement]
                bg.alleles[AlleleState.NORMAL] = [
                    Pair(first, second), Pair(replacement, first), Pair(replacement, second)
                ]
                if zygosity is None:
                    bg.variant_pool.pop(self.ADDITIONAL_HOM)
                else:
                    bg.variant_pool[self.ADDITIONAL_HOM] = zygosity

                self._assert_unchanged(bg, ensure_HET_SNP_used)

    def test_missing_replacement_does_not_invent_an_allele(self):
        bg = self._group()
        replacement = self.curated["ABO*O.01.02"]
        bg.alleles[AlleleState.RAW].remove(replacement)
        bg.alleles[AlleleState.NORMAL] = [
            pair for pair in bg.alleles[AlleleState.NORMAL] if replacement not in pair
        ]

        self._assert_unchanged(bg, ensure_HET_SNP_used)

    def test_without_a_het_this_stage_does_not_apply_the_new_rule(self):
        bg = self._group()
        bg.variant_pool[self.HET] = Zygosity.HOM

        self._assert_unchanged(bg, ensure_HET_SNP_used)


class TestFUT3ConditionalCoverage(_CoverageFixtures):
    """Known partner tokens can force a replacement despite its unphased additions."""

    REASON = "filter_on_in_relationship_when_HOM_cant_be_on_one_side"
    STAGE = staticmethod(filter_on_in_relationship_when_HOM_cant_be_on_one_side)

    def _na21101(self):
        names = ("FUT3*01.17.01", "FUT3*01N.17.03", "FUT3*01N.17.05", "FUT3*01N.17.19")
        hom, left, unphased, right = [self.curated[name] for name in names]
        return BloodGroup(
            type="FUT3", sample="NA21101_conditional_coverage",
            alleles={AlleleState.RAW: [hom, left, unphased, right],
                     AlleleState.NORMAL: [Pair(left, unphased), Pair(left, right),
                                          Pair(hom, left), Pair(hom, unphased), Pair(hom, right)],
                     AlleleState.CO: None},
            variant_pool={"19:5844781_A_C": Zygosity.HOM, "19:5843773_A_T": Zygosity.HET,
                          "19:5843779_C_T": Zygosity.HET, "19:5844332_C_T": Zygosity.HET},
            variant_pool_phase={"19:5844781_A_C": "1/1", "19:5843773_A_T": "0|1",
                                "19:5843779_C_T": "1|0", "19:5844332_C_T": "0/1"},
            variant_pool_phase_set={"19:5844781_A_C": ".", "19:5843773_A_T": "5843773",
                                    "19:5843779_C_T": "5843773", "19:5844332_C_T": "unknown"},
        )

    def _hg03388(self):
        names = [f"FUT3*01N.01.{suffix}" for suffix in ("01", "02", "08", "10", "11", "12")]
        hom, d, cd, ad, bc, b = [self.curated[name] for name in names]
        return BloodGroup(
            type="FUT3", sample="HG03388_conditional_coverage",
            alleles={AlleleState.RAW: [hom, d, cd, ad, bc, b],
                     AlleleState.NORMAL: [Pair(hom, d), Pair(hom, cd), Pair(hom, ad),
                                          Pair(hom, bc), Pair(hom, b), Pair(d, bc),
                                          Pair(d, b), Pair(cd, b), Pair(ad, bc), Pair(ad, b)],
                     AlleleState.CO: None},
            variant_pool={"19:5844356_C_T": Zygosity.HOM, "19:5844827_C_T": Zygosity.HOM,
                          "19:5843866_G_A": Zygosity.HET, "19:5843872_C_G": Zygosity.HET,
                          "19:5844032_C_T": Zygosity.HET, "19:5844173_C_T": Zygosity.HET},
            variant_pool_phase={"19:5844356_C_T": "1/1", "19:5844827_C_T": "1/1",
                                "19:5843866_G_A": "0|1", "19:5843872_C_G": "1|0",
                                "19:5844032_C_T": "0/1", "19:5844173_C_T": "0/1"},
            variant_pool_phase_set={"19:5844356_C_T": ".", "19:5844827_C_T": ".",
                                    "19:5843866_G_A": "5843866", "19:5843872_C_G": "5843866",
                                    "19:5844032_C_T": "unknown", "19:5844173_C_T": "unknown"},
        )

    def test_na21101_rejects_only_the_forced_underdescribed_pair(self):
        bg = self._na21101()
        before = list(bg.alleles[AlleleState.NORMAL])
        bad = self._pair("FUT3*01.17.01", "FUT3*01N.17.19")
        self.assertEqual(bg.variant_pool_phase["19:5844332_C_T"], "0/1")
        self.assertIn("19:5844332_C_T", self.curated["FUT3*01N.17.19"].defining_variants)

        self.STAGE({"FUT3": bg}, phased=True)

        self.assertEqual(bg.alleles[AlleleState.NORMAL], [pair for pair in before if pair != bad])
        self.assertEqual(set(bg.filtered_out), {self.REASON})
        self.assertEqual(bg.filtered_out[self.REASON], [bad])
        self.assertIn(self._pair("FUT3*01.17.01", "FUT3*01N.17.03"), bg.alleles[AlleleState.NORMAL])

    def test_hg03388_retains_the_other_two_new_pairs(self):
        bg = self._hg03388()
        before = list(bg.alleles[AlleleState.NORMAL])
        bad = self._pair("FUT3*01N.01.01", "FUT3*01N.01.10")
        self.assertEqual(bg.variant_pool_phase["19:5844173_C_T"], "0/1")
        self.assertIn("19:5844173_C_T", self.curated["FUT3*01N.01.10"].defining_variants)

        self.STAGE({"FUT3": bg}, phased=True)

        self.assertEqual(bg.alleles[AlleleState.NORMAL], [pair for pair in before if pair != bad])
        self.assertEqual(set(bg.filtered_out), {self.REASON})
        self.assertEqual(bg.filtered_out[self.REASON], [bad])
        for partner in ("FUT3*01N.01.11", "FUT3*01N.01.12"):
            self.assertIn(self._pair("FUT3*01N.01.01", partner), bg.alleles[AlleleState.NORMAL])

    def test_extra_must_be_proved_opposite_in_a_known_shared_block(self):
        for field, value in (
            ("variant_pool_phase_set", "unknown"),
            ("variant_pool_phase_set", "9000000"),
            ("variant_pool_phase", "1|0"),
        ):
            with self.subTest(field=field, value=value):
                bg = self._na21101()
                getattr(bg, field)["19:5843773_A_T"] = value

                self._assert_unchanged(bg, self.STAGE, phased=True)

    def test_additional_hem_or_no_data_token_cannot_force_replacement(self):
        for zygosity in (Zygosity.HEM, Zygosity.NO_DATA):
            with self.subTest(zygosity=zygosity):
                bg = self._na21101()
                bg.variant_pool["19:5843773_A_T"] = zygosity

                self._assert_unchanged(bg, self.STAGE, phased=True)

    def test_missing_forced_replacement_cannot_be_replaced_by_an_unphased_extension(self):
        bg = self._hg03388()
        replacement = self.curated["FUT3*01N.01.12"]
        bg.alleles[AlleleState.RAW].remove(replacement)
        bg.alleles[AlleleState.NORMAL] = [
            pair for pair in bg.alleles[AlleleState.NORMAL] if replacement not in pair
        ]
        self.assertIn(self.curated["FUT3*01N.01.11"], bg.alleles[AlleleState.RAW])

        self._assert_unchanged(bg, self.STAGE, phased=True)

    def test_unphased_flag_leaves_both_cases_unchanged(self):
        for factory in (self._na21101, self._hg03388):
            with self.subTest(case=factory.__name__):
                self._assert_unchanged(factory(), self.STAGE, phased=False)


if __name__ == "__main__":
    unittest.main()
