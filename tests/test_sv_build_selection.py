"""Keep structural matching within the selected genome build."""

import unittest

import pandas as pd

from rbceq2.core_logic.constants import AlleleState
from rbceq2.core_logic.utils import sub_alleles_relationships
from rbceq2.db.db import Db, build_antigen_map_for_checks, prepare_db
from rbceq2.main import find_hits, parse_args


class TestSvBuildSelection(unittest.TestCase):
    """A wrong-build exact match must not displace an accepted local match.

    These synthetic deletions exercise the existing matching thresholds against
    curated GCNT2 definitions; they do not establish breakpoint equivalence.
    """

    @classmethod
    def setUpClass(cls):
        """Load the curated database and pipeline relationships for both builds."""
        cls.contexts = {}
        for build in ("GRCh37", "GRCh38"):
            db = Db(ref=build, df=prepare_db())
            kn = {"KN": [a for a in db.make_alleles() if a.blood_group == "KN"]}
            relationships, bg_type = sub_alleles_relationships(kn, "KN")
            cls.contexts[build] = (
                db, {bg_type: relationships}, build_antigen_map_for_checks(db.df)
            )

    def assert_selected_build_result(self, build, input_build):
        """Check matching, allele selection, genotype and both phenotypes."""
        db, relationships, ant_mapping = self.contexts[build]
        row = db.df.loc[db.df["Genotype"] == "GCNT2*01N.06"].iloc[0]
        selected_token = row[build].strip()
        position = int(row[input_build].strip().split("_")[0])
        frame = pd.DataFrame(
            [["chr6", str(position), ".", "N", "<DEL>", "50", "PASS",
              f"SVTYPE=DEL;END={position + 41000};SVLEN=-41000", "GT", "1/1"]],
            columns=["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER",
                     "INFO", "FORMAT", "SAMPLE"],
        )
        result = find_hits(
            db, (frame, f"{build}_at_{input_build}_position"),
            args=parse_args(["--reference_genome", build, "--HPAs"]),
            allele_relationships=relationships, excluded=["RHD", "RHCE"],
            ant_mapping=ant_mapping,
        )
        _, genotypes, numeric, alphanumeric, groups, var_map = result
        expected = "GCNT2*01N.06/GCNT2*01N.06"
        with self.subTest(output="selected-build token"):
            self.assertIn(f"6:{selected_token}", var_map)
            other = "GRCh38" if build == "GRCh37" else "GRCh37"
            self.assertNotIn(f"6:{row[other].strip()}", var_map)
        with self.subTest(output="internal pairs"):
            self.assertCountEqual(
                ["/".join(p.genotypes)
                 for p in groups["GCNT2"].alleles[AlleleState.NORMAL] or []],
                [expected],
            )
        with self.subTest(output="genotype"):
            self.assertEqual(genotypes["GCNT2"], expected)
        with self.subTest(output="numeric phenotype"):
            self.assertEqual(numeric["GCNT2"], "I:-1")
        with self.subTest(output="alphanumeric phenotype"):
            self.assertEqual(alphanumeric["GCNT2"], "I-")

    def test_grch37_exact_position(self):
        self.assert_selected_build_result("GRCh37", "GRCh37")

    def test_grch37_nearby_position_matching_other_build_exactly(self):
        self.assert_selected_build_result("GRCh37", "GRCh38")

    def test_grch38_exact_position(self):
        self.assert_selected_build_result("GRCh38", "GRCh38")

    def test_grch38_nearby_position_matching_other_build_exactly(self):
        self.assert_selected_build_result("GRCh38", "GRCh37")


if __name__ == "__main__":
    unittest.main()
