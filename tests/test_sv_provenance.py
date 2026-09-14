"""Keep a matched structural event's evidence together through inference."""

import unittest
from types import SimpleNamespace

import pandas as pd

from rbceq2.core_logic.alleles import BloodGroup
from rbceq2.core_logic.constants import AlleleState
from rbceq2.core_logic.data_procesing import (
    filter_value_for,
    filter_values_for,
    record_unused_variants,
    variant_was_discarded,
)
from rbceq2.core_logic.utils import Zygosity, sub_alleles_relationships
from rbceq2.db.db import Db, build_antigen_map_for_checks, prepare_db
from rbceq2.main import find_hits, parse_args


class TestSelectedSvPipeline(unittest.TestCase):
    """Use curated GCNT2 definitions and competing synthetic source rows."""

    @classmethod
    def setUpClass(cls):
        cls.db = Db(ref="GRCh38", df=prepare_db())
        kn = {"KN": [a for a in cls.db.make_alleles() if a.blood_group == "KN"]}
        relationships, bg_type = sub_alleles_relationships(kn, "KN")
        cls.relationships = {bg_type: relationships}
        cls.ant_mapping = build_antigen_map_for_checks(cls.db.df)
        row = cls.db.df.loc[cls.db.df["Genotype"] == "GCNT2*01N.06"].iloc[0]
        cls.token = "6:" + row["GRCh38"].strip()
        cls.position = int(row["GRCh38"].strip().split("_")[0])

    def run_rows(self, selected_filter, reverse=False, no_filter=False,
                 missing_ps=False, selected_gt="0|1"):
        """Run both source rows, preserving their independent evidence."""
        lower_filter = "LowQual" if selected_filter == "PASS" else "PASS"
        rows = []
        for length, quality, fmt, sample in (
            (41000, selected_filter, "GT" if missing_ps else "GT:PS",
             selected_gt if missing_ps else f"{selected_gt}:10"),
            (42900 if missing_ps else 41900, lower_filter, "GT:PS",
             "1|1:77" if missing_ps else "1|1:20"),
        ):
            rows.append([
                "chr6", str(self.position), ".", "N", "<DEL>", "50",
                quality,
                f"SVTYPE=DEL;END={self.position + length};SVLEN=-{length}",
                fmt, sample,
            ])
        frame = pd.DataFrame(
            rows[::-1] if reverse else rows,
            columns=["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER",
                     "INFO", "FORMAT", "SAMPLE"],
        )
        flags = ["--reference_genome", "GRCh38", "--HPAs", "--phased"]
        if no_filter:
            flags.append("--no_filter")
        return find_hits(
            self.db, (frame, "selected_sv_source"), args=parse_args(flags),
            allele_relationships=self.relationships,
            excluded=["RHD", "RHCE"], ant_mapping=self.ant_mapping,
        )

    def assert_selected_heterozygote(self, result, phase_set="10"):
        """Check the result and the source evidence that justified it."""
        _, genotypes, numeric, alpha, groups, _ = result
        bg = groups["GCNT2"]
        expected = ["GCNT2*01", "GCNT2*01N.06"]
        self.assertCountEqual(genotypes["GCNT2"].split("/"), expected)
        self.assertEqual(numeric["GCNT2"], "I:1")
        self.assertEqual(alpha["GCNT2"], "I+")
        pairs = bg.alleles[AlleleState.NORMAL]
        self.assertEqual(len(pairs), 1)
        self.assertCountEqual(pairs[0].genotypes, expected)
        self.assertEqual(bg.variant_pool[self.token], Zygosity.HET)
        self.assertEqual(bg.variant_pool_phase[self.token], "0|1")
        self.assertEqual(bg.variant_pool_phase_set[self.token], phase_set)

    def test_selected_pass_survives_lower_ranked_lowqual_in_both_orders(self):
        for reverse in (False, True):
            with self.subTest(reverse=reverse):
                self.assert_selected_heterozygote(self.run_rows("PASS", reverse))

    def test_selected_lowqual_is_excluded_despite_lower_ranked_pass(self):
        for reverse in (False, True):
            with self.subTest(reverse=reverse):
                _, genotypes, _, _, groups, _ = self.run_rows("LowQual", reverse)
                bg = groups["GCNT2"]
                self.assertEqual(genotypes["GCNT2"], "GCNT2*01/GCNT2*01")
                self.assertNotIn(self.token, bg.variant_pool)
                self.assertIn(
                    "GCNT2*01N.06",
                    [a.genotype for a in bg.filtered_out["FILTER_not_PASS"]],
                )

    def test_selected_homozygote_retains_null_phenotypes(self):
        """The selected passing homozygote must retain both null phenotypes."""
        for reverse in (False, True):
            with self.subTest(reverse=reverse):
                _, genotype, numeric, alpha, _, _ = self.run_rows(
                    "PASS", reverse, selected_gt="1|1"
                )
                self.assertEqual(genotype["GCNT2"], "GCNT2*01N.06/GCNT2*01N.06")
                self.assertEqual(numeric["GCNT2"], "I:-1")
                self.assertEqual(alpha["GCNT2"], "I-")

    def test_no_filter_keeps_selected_genotype_and_phase(self):
        for selected_filter in ("PASS", "LowQual"):
            for reverse in (False, True):
                with self.subTest(selected_filter=selected_filter, reverse=reverse):
                    self.assert_selected_heterozygote(
                        self.run_rows(selected_filter, reverse, no_filter=True)
                    )

    def test_missing_selected_ps_is_not_borrowed_from_another_source_row(self):
        for reverse in (False, True):
            with self.subTest(reverse=reverse):
                self.assert_selected_heterozygote(
                    self.run_rows("PASS", reverse, missing_ps=True), phase_set="."
                )


class TestSelectedSvFilterEvidence(unittest.TestCase):
    """Selected DB tokens override ambiguous rows without changing other lookups."""

    def setUp(self):
        self.token = "6:10585989_del_41kb"
        self.small = "6:10586805_C_G"
        self.df = pd.DataFrame({
            "variant": [self.token, self.token, self.small],
            "FILTER": ["LowQual", "PASS", "PASS"],
        })

    def test_selected_filter_overrides_first_row_and_preserves_fallback(self):
        selected = {self.token: "PASS"}
        self.assertEqual(filter_values_for(
            self.token, self.df, selected_sv_filters=selected), ["PASS"])
        self.assertEqual(filter_value_for(
            self.token, self.df, selected_sv_filters=selected), "PASS")
        self.assertEqual(filter_values_for(self.token, self.df), ["LowQual", "PASS"])
        self.assertEqual(filter_values_for(
            self.small, self.df, selected_sv_filters=selected), ["PASS"])
        self.assertIsNone(filter_value_for(
            "6:1_A_G", self.df, selected_sv_filters=selected))

    def test_discarded_lookup_uses_selected_filter(self):
        for value, expected in (("PASS", False), ("LowQual", True)):
            with self.subTest(value=value):
                self.assertEqual(variant_was_discarded(
                    self.token, self.df, selected_sv_filters={self.token: value}
                ), expected)

    def test_unused_audit_keeps_selected_gt_ps_and_filter_together(self):
        bg = BloodGroup(type="GCNT2", alleles={}, sample="source_audit")
        vcf = SimpleNamespace(variants={self.token: {"GT": "0|1", "PS": "10"}})
        record_unused_variants(
            {"GCNT2": bg}, vcf=vcf, df=self.df,
            loci_by_type={"GCNT2": {"6": frozenset({10585989})}},
            selected_sv_filters={self.token: "PASS"},
        )
        self.assertEqual(bg.unused_pool[self.token], Zygosity.HET)
        self.assertEqual(bg.unused_pool_phase[self.token], "0|1")
        self.assertEqual(bg.unused_pool_phase_set[self.token], "10")
        self.assertEqual(bg.unused_pool_filters[self.token], "PASS")


if __name__ == "__main__":
    unittest.main()
