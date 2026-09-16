"""Infer Lewis phenotypes from curated FUT3 and secretor alleles."""

import unittest

import pandas as pd

from rbceq2.core_logic.constants import AlleleState
from rbceq2.core_logic.utils import sub_alleles_relationships
from rbceq2.db.db import Db, build_antigen_map_for_checks, prepare_db
from rbceq2.main import find_hits, parse_args


class TestFut3Secretor(unittest.TestCase):
    """Keep active Lewis independent of H expression through inference."""

    @classmethod
    def setUpClass(cls):
        cls.contexts = {}
        for build in ("GRCh37", "GRCh38"):
            db = Db(ref=build, df=prepare_db())
            kn = {"KN": [a for a in db.make_alleles() if a.blood_group == "KN"]}
            relationships, bg_type = sub_alleles_relationships(kn, "KN")
            cls.contexts[build] = (
                db, {bg_type: relationships}, build_antigen_map_for_checks(db.df),
            )

    def assert_curated_result(self, build, h_allele, se_allele, h_status,
                              se_status, le_status, active=True):
        """Run the curated alleles without changing their database definitions."""
        db, relationships, mapping = self.contexts[build]

        def definition(genotype):
            return db.df.loc[db.df["Genotype"] == genotype, build].iloc[0]

        def row(position, ref, alt, gt):
            return ["chr19", position, ".", ref, alt, "50", "PASS", ".", "GT", gt]

        # These reverse-strand substitutions restore transcript reference bases.
        # Homozygous ALT prevents the implicit genome-reference lane tokens.
        # With 0/0, the same loci instead supply FUT3*01N.03.02's two tokens.
        lane_gt = "1/1" if active else "0/0"
        rows = [
            row(definition("FUT3*01.04").split("_")[0], "A", "G", lane_gt),
            row(definition("FUT3*01N.03.01").split("_")[0], "G", "A", lane_gt),
        ]
        for allele in (h_allele, se_allele):
            if allele in ("FUT1*01", "FUT2*01"):
                continue
            for token in definition(allele).split(","):
                position, ref, alt = token.split("_")
                rows.append(row(position, ref, alt, "1/1"))
        rows.sort(key=lambda values: int(values[1]))
        frame = pd.DataFrame(rows, columns=[
            "CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO",
            "FORMAT", "SAMPLE",
        ])
        result = find_hits(
            db, (frame, f"Lewis_{build}_{h_status}_{se_status}_{active}"),
            args=parse_args(["--reference_genome", build, "--HPAs"]),
            allele_relationships=relationships, excluded=["RHD", "RHCE"],
            ant_mapping=mapping,
        )
        _, genotypes, numeric, alpha, groups, _ = result
        le_allele = "FUT3*01.01" if active else "FUT3*01N.03.02"
        expected_h = f"{h_allele}/{h_allele}"
        expected_se = f"{se_allele}/{se_allele}"
        expected_le = f"{le_allele}/{le_allele}"
        for group, expected in (("FUT1", expected_h), ("FUT2", expected_se),
                                ("FUT3", expected_le)):
            with self.subTest(group=group, output="internal pairs"):
                self.assertCountEqual(
                    ["/".join(pair.genotypes)
                     for pair in groups[group].alleles[AlleleState.NORMAL] or []],
                    [expected],
                )
        self.assertEqual(genotypes["FUT3"], expected_le)
        for group in ("FUT1", "FUT2"):
            self.assertCountEqual(genotypes[group].split(","),
                                  [expected_h, expected_se])
        self.assertEqual(alpha["FUT3"], le_status)
        self.assertEqual(numeric["FUT3"], "")
        self.assertEqual(alpha["FUT2"], se_status)
        self.assertEqual(alpha["FUT1"], f"{h_status},{se_status}")

    def test_curated_h_and_secretor_combinations(self):
        cases = (
            ("FUT1*01", "FUT2*01", "H+", "Se+", "Le(a-b+)"),
            ("FUT1*01N.01", "FUT2*01", "H-", "Se+", "Le(a-b+)"),
            ("FUT1*01W.01", "FUT2*01", "H+w", "Se+", "Le(a-b+)"),
            ("FUT1*01", "FUT2*01W.02.01", "H+", "Se+w", "Le(a+b+)"),
            ("FUT1*01", "FUT2*01N.01", "H+", "Se-", "Le(a+b-)"),
        )
        for build in self.contexts:
            for case in cases:
                with self.subTest(build=build, h_allele=case[0], se_allele=case[1]):
                    self.assert_curated_result(build, *case)

    def test_curated_null_lewis_with_weak_secretor(self):
        for build in self.contexts:
            with self.subTest(build=build):
                self.assert_curated_result(
                    build, "FUT1*01", "FUT2*01W.02.01", "H+", "Se+w",
                    "Le(a-b-)", active=False,
                )
