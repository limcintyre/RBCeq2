"""Preserve LU reference-background inference without inventing measured calls.

The curated LU*02.19 compound definition contains the LU lane reference token
and a separate alternate modifier. Fully missing lane evidence follows the
project's reference-background convention; missing modifier evidence must never
supply that alternate modifier. These cases exercise both genome builds through
the production plain and gzip VCF reader.
"""

import gzip
import tempfile
import unittest
from pathlib import Path

from rbceq2.IO.vcf import VCF, build_intervals, read_vcf
from rbceq2.core_logic.constants import AlleleState
from rbceq2.core_logic.utils import Zygosity, sub_alleles_relationships
from rbceq2.db.db import Db, build_antigen_map_for_checks, prepare_db
from rbceq2.main import find_hits, parse_args


class TestLuReferenceFallback(unittest.TestCase):
    """Distinguish reference fallback from measured lane evidence."""

    @classmethod
    def setUpClass(cls):
        """Load the curated compound alleles and create owned input storage."""
        directory = tempfile.TemporaryDirectory(prefix="rbceq2_lu_reference_")
        cls.addClassCleanup(directory.cleanup)
        cls.directory = Path(directory.name)
        cls.contexts = {}
        cls.results = {}
        for build in ("GRCh37", "GRCh38"):
            db = Db(ref=build, df=prepare_db())
            rows = db.df.set_index("Genotype")
            lane = rows.loc["LU*01", build].strip()
            lu01 = {token.strip() for token in
                    rows.loc["LU*01.19", build].split(",")}
            lu02 = {token.strip() for token in
                    rows.loc["LU*02.19", build].split(",")}
            modifier, = lu01 & lu02
            # Check fixture relationships against the database, not literals.
            if lu01 != {lane, modifier} or lu02 != {
                f"{lane.split('_')[0]}_ref", modifier
            }:
                raise AssertionError("Curated LU compound definitions changed")
            kn = {"KN": [a for a in db.make_alleles()
                         if a.blood_group == "KN"]}
            relationships, group = sub_alleles_relationships(kn, "KN")
            cls.contexts[build] = {
                "db": db, "relationships": {group: relationships},
                "mapping": build_antigen_map_for_checks(db.df),
                "chrom": rows.loc["LU*01", "Chrom"].removeprefix("chr"),
                "lane": lane.split("_"), "modifier": modifier.split("_"),
                "rows": rows,
            }

    def write_input(self, build, compressed, lane_kind, lane_gt, modifier_gt):
        """Write nonoverlapping LU rows using database-derived coordinates."""
        context = self.contexts[build]
        chrom = context["chrom"]
        position, ref, alt = context["modifier"]
        rows = [[chrom, position, ".", ref, alt, "50", "PASS", ".",
                 "GT", modifier_gt]]
        if lane_kind != "omitted":
            position, ref, alt = context["lane"]
            info = "."
            if lane_kind in ("point", "interior"):
                target = int(position)
                start, end = ((target, target) if lane_kind == "point"
                              else (target - 3, target + 3))
                self.assertFalse(start <= int(context["modifier"][0]) <= end)
                position = str(start)
                ref = ref if start == target else "N"
                alt = "<NON_REF>"
                info = f"END={end}"
            rows.append([chrom, position, ".", ref, alt, "50", "PASS",
                         info, "GT", lane_gt])
        content = (
            "##fileformat=VCFv4.2\n"
            '##INFO=<ID=END,Number=1,Type=Integer,Description="End">\n'
            '##ALT=<ID=NON_REF,Description="Unspecified alternate allele">\n'
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
            + "".join("\t".join(row) + "\n"
                      for row in sorted(rows, key=lambda row: int(row[1])))
        )
        name = f"{build}_{lane_kind}_{lane_gt}_{modifier_gt}".replace("/", "_")
        path = self.directory / f"{name}.vcf{'.gz' if compressed else ''}"
        opener = gzip.open if compressed else open
        with opener(path, "wt") as handle:
            handle.write(content)
        return path

    def infer(self, build, compressed, lane_kind="ordinary", lane_gt="0/0",
              modifier_gt="1/1"):
        """Return calls, pairs, and evidence from the file path API."""
        key = (build, compressed, lane_kind, lane_gt, modifier_gt)
        if key not in self.results:
            context = self.contexts[build]
            path = self.write_input(*key)
            original = path.read_bytes()
            _, genotype, numeric, alpha, groups, _ = find_hits(
                context["db"], path,
                args=parse_args(["--reference_genome", build, "--HPAs",
                                 "--debug"]),
                allele_relationships=context["relationships"],
                excluded=["RHD", "RHCE"], ant_mapping=context["mapping"],
            )
            self.assertEqual(path.read_bytes(), original)
            group = groups.get("LU")
            pairs = [] if group is None else [
                frozenset(pair.genotypes)
                for pair in group.alleles[AlleleState.NORMAL] or []
            ]
            evidence = None if group is None else {
                **group.variant_pool, **group.unused_pool
            }
            calls = (genotype.get("LU"), numeric.get("LU"), alpha.get("LU"))
            self.results[key] = calls, pairs, evidence, path
        return self.results[key]

    def assert_call(self, build, result, alleles, numeric=None, alpha=None):
        """Check reported alleles, pairs, and phenotype token sets."""
        calls, pairs, _, _ = result
        self.assertCountEqual(calls[0].split("/"), alleles)
        self.assertEqual(pairs, [frozenset(alleles)])
        if numeric is None:
            row = self.contexts[build]["rows"].loc[alleles[0]]
            numeric = row["Phenotype_change"]
            # Output combines each antigen pair into one display token.
            alpha = {"LU*01.19": "Lu(a+b-),Au(a-b+)",
                     "LU*02.19": "Lu(a-b+),Au(a-b+)"}[alleles[0]]
        self.assertEqual(set(calls[1].split(",")), set(numeric.split(",")))
        self.assertEqual(set(calls[2].split(",")), set(alpha.split(",")))

    def assert_missing_evidence(self, build, result, target):
        """Retain missing state without measured chromosome copies."""
        context = self.contexts[build]
        _, _, evidence, path = result
        chrom = context["chrom"]
        position = int(context[target][0])
        prefix = f"{chrom}:{position}_"
        db = context["db"]
        frame = read_vcf(str(path), build_intervals(db.df, build),
                         db.unique_variants)
        vcf = VCF(frame, db.lane_variants, db.unique_variants,
                  path.stem, build)
        rows = vcf.df.loc[vcf.df["POS"].astype(int) == position]
        self.assertFalse(rows.empty)
        source_gts = {value.split(":")[0] for value in rows["SAMPLE"]}
        self.assertIn("./.", source_gts)
        self.assertFalse(source_gts & {"0/0", "0|0"})
        self.assertNotIn(position, vcf.haploid_loci.get(chrom, ()))
        self.assertNotIn(position, vcf.diploid_loci.get(chrom, ()))
        if evidence is not None:
            self.assertTrue(any(
                token.startswith(prefix) and state == Zygosity.NO_DATA
                for token, state in evidence.items()
            ))

    def check_missing_lane(self, kind):
        """Keep the called modifier on reference background."""
        for build in self.contexts:
            for compressed in (False, True):
                with self.subTest(build=build, compressed=compressed):
                    result = self.infer(build, compressed, kind, "./.")
                    with self.subTest(check="missing evidence"):
                        self.assert_missing_evidence(build, result, "lane")
                    with self.subTest(check="reported result"):
                        self.assert_call(
                            build, result, ["LU*02.19", "LU*02.19"]
                        )
                        context = self.contexts[build]
                        modifier = "_".join(context["modifier"])
                        self.assertIsNotNone(result[2])
                        self.assertEqual(
                            result[2].get(f'{context["chrom"]}:{modifier}'),
                            Zygosity.HOM,
                        )
                        control = self.infer(build, compressed)
                        self.assertEqual(result[0], control[0])

    def test_omitted_and_reference_lane_controls(self):
        """Omitted and reference lanes support the LU*02.19 pair."""
        for build in self.contexts:
            for compressed in (False, True):
                for kind in ("omitted", "ordinary"):
                    with self.subTest(build=build, compressed=compressed,
                                      kind=kind):
                        result = self.infer(build, compressed, kind)
                        self.assert_call(
                            build, result, ["LU*02.19", "LU*02.19"]
                        )

    def test_missing_point_lane(self):
        """A one-base missing gVCF block permits reference fallback."""
        self.check_missing_lane("point")

    def test_missing_interior_lane(self):
        """A missing block containing the lane permits the same fallback."""
        self.check_missing_lane("interior")

    def test_missing_concrete_alt_lane(self):
        """A missing ordinary lane row follows the same convention."""
        self.check_missing_lane("ordinary")

    def test_called_heterozygous_lane(self):
        """A heterozygous lane requires one compound allele per background."""
        for build in self.contexts:
            for compressed in (False, True):
                with self.subTest(build=build, compressed=compressed):
                    result = self.infer(build, compressed, lane_gt="0/1")
                    self.assert_call(
                        build, result, ["LU*01.19", "LU*02.19"],
                        "LU:1,2,-18,19", "Lu(a+b+),Au(a-b+)"
                    )

    def test_called_homozygous_lane(self):
        """A called alternate lane cannot fall back to LU*02."""
        for build in self.contexts:
            for compressed in (False, True):
                with self.subTest(build=build, compressed=compressed):
                    result = self.infer(build, compressed, lane_gt="1/1")
                    self.assert_call(build, result, ["LU*01.19", "LU*01.19"])

    def test_missing_modifier_does_not_supply_alternate(self):
        """Unknown modifier evidence cannot construct either LU*.19 allele."""
        for build in self.contexts:
            for compressed in (False, True):
                with self.subTest(build=build, compressed=compressed):
                    result = self.infer(build, compressed, modifier_gt="./.")
                    self.assert_missing_evidence(build, result, "modifier")
                    self.assertNotIn("LU*02.19", result[0][0].split("/"))
                    self.assertNotIn("LU*01.19", result[0][0].split("/"))
                    self.assertTrue(all(not pair & {"LU*01.19", "LU*02.19"}
                                        for pair in result[1]))


if __name__ == "__main__":
    unittest.main()
