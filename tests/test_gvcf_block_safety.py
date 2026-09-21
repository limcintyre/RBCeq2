"""Guard block projection against fabricated alleles and phase evidence."""

import unittest

import pandas as pd

from rbceq2.IO.vcf import VCF, is_single_copy
from rbceq2.core_logic.constants import PAR
from rbceq2.core_logic.utils import BeyondLogicError
from rbceq2.db.db import Db, prepare_db


class TestGvcfBlockSafety(unittest.TestCase):
    """Check bounded block handling through the existing sample VCF API."""

    @classmethod
    def setUpClass(cls):
        """Load the curated YTLI locus without changing its definitions."""
        cls.db = Db(ref="GRCh38", df=prepare_db())
        allele = cls.db.df.loc[
            cls.db.df["Genotype"] == "YT*01.-04"
        ].iloc[0]
        token = next(t.strip() for t in allele["GRCh38"].split(",")
                     if not t.strip().endswith("_ref"))
        pos, cls.ref, cls.alt = token.split("_")
        cls.pos = int(pos)
        cls.locus = f"7:{pos}"

    def block(self, gt="./.", info=None, alt="<NON_REF>", fmt="GT"):
        """Return a synthetic block beginning outside the target locus."""
        if info is None:
            info = f"END={self.pos + 3}"
        return ["7", str(self.pos - 3), ".", "A", alt, ".", ".",
                info, fmt, gt]

    def read(self, rows):
        """Construct a production VCF from independent test inputs."""
        frame = pd.DataFrame(rows, columns=[
            "CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER",
            "INFO", "FORMAT", "SAMPLE",
        ])
        return VCF([frame], self.db.lane_variants, self.db.unique_variants,
                   "block_safety", "GRCh38")

    def test_invalid_end_is_named(self):
        """Malformed intervals cannot silently act as point rows."""
        for info in ("END=bad", "END=", f"END={self.pos - 4}",
                     f"END={self.pos},other",
                     f"END={self.pos};END={self.pos}"):
            with self.subTest(info=info):
                with self.assertRaisesRegex(BeyondLogicError,
                                            "VCF/invalid_gVCF_END"):
                    self.read([self.block(info=info)])

    def test_called_placeholder_and_partial_gt_are_named(self):
        """Unspecified or incomplete calls are not block reference."""
        for gt in ("1", "0/1", "1/1", "0/."):
            with self.subTest(gt=gt):
                with self.assertRaisesRegex(BeyondLogicError,
                                            "VCF/unsupported_gVCF_block_GT"):
                    self.read([self.block(gt=gt)])

    def test_explicit_row_overlap_refuses_in_both_orders(self):
        """Overlaps must not make row order authoritative."""
        explicit = ["7", str(self.pos), ".", self.ref, self.alt,
                    "50", "PASS", ".", "GT", "0/1"]
        block = self.block()
        for rows in ([block, explicit], [explicit, block]):
            with self.subTest(reverse=rows[0] is explicit):
                with self.assertRaisesRegex(
                    BeyondLogicError, "VCF/ambiguous_gVCF_block_overlap"
                ):
                    self.read(rows)

    def test_conflicting_block_ploidy_refuses_in_both_orders(self):
        """Opposing copy observations cannot be resolved by row ordering."""
        one, two = self.block("0"), self.block("0/0")
        for rows in ([one, two], [two, one]):
            with self.subTest(reverse=rows[0] is two):
                with self.assertRaisesRegex(
                    BeyondLogicError, "VCF/ambiguous_gVCF_block_overlap"
                ):
                    self.read(rows)

    def test_identical_blocks_do_not_duplicate_target_evidence(self):
        """Identical blocks have the same target observations."""
        row = self.block()
        once, twice = self.read([row]), self.read([row, list(row)])
        self.assertEqual(once.variants, twice.variants)
        target = once.df.loc[once.df["POS"].astype(int) == self.pos]
        self.assertEqual(len(target), 1)
        self.assertNotEqual(target.iloc[0]["REF"], row[3])
        self.assertEqual(target.iloc[0]["SAMPLE"], "./.")
        sources = once.gvcf_block_sources
        self.assertEqual(sources[0]["INFO"], row[7])
        self.assertEqual(sources[0]["SAMPLE"], row[9])

    def test_block_ps_does_not_create_variant_phase_span(self):
        """Block PS must not assert haplotype links at new positions."""
        for gt in ("0|0", ".|."):
            with self.subTest(gt=gt):
                vcf = self.read([self.block(gt=f"{gt}:123", fmt="GT:PS")])
                self.assertEqual(vcf.phase_sets, {})

    def test_par_boundary_preserves_coordinate_specific_copy_evidence(self):
        """A missing block cannot supply chromosome or locus copy evidence."""
        for build in ("GRCh37", "GRCh38"):
            boundary = PAR[build]["X"][0][1]
            for gt in ("0", "0/0", ".", "./."):
                with self.subTest(build=build, gt=gt):
                    frame = pd.DataFrame([[
                        "X", str(boundary), ".", "N", "<NON_REF>", ".", ".",
                        f"END={boundary + 1}", "GT", gt,
                    ]], columns=["CHROM", "POS", "ID", "REF", "ALT", "QUAL",
                                 "FILTER", "INFO", "FORMAT", "SAMPLE"])
                    vcf = VCF([frame], {},
                              {f"X:{boundary}", f"X:{boundary + 1}"},
                              "boundary", build)
                    self.assertFalse(is_single_copy("X", boundary, build))
                    self.assertTrue(is_single_copy("X", boundary + 1, build))
                    self.assertEqual(vcf.haploid_chroms,
                                     frozenset({"X"}) if gt == "0"
                                     else frozenset())
                    if "." in gt:
                        self.assertEqual(vcf.haploid_loci, {})
                        self.assertEqual(vcf.diploid_loci, {})
                    else:
                        observed = (vcf.haploid_loci if gt == "0"
                                    else vcf.diploid_loci)
                        self.assertEqual(observed["X"],
                                         frozenset({boundary, boundary + 1}))

    def test_other_alt_forms_are_not_reference_blocks(self):
        """END cannot turn structural or mixed ALT rows into blocks."""
        for alt in ("<DEL>", "*", "T,<NON_REF>"):
            with self.subTest(alt=alt):
                vcf = self.read([self.block(gt="0/0", alt=alt)])
                self.assertNotIn(self.pos, vcf.diploid_loci.get("7", ()))


if __name__ == "__main__":
    unittest.main()
