"""Preserve measured and missing YTLI evidence across inclusive gVCF blocks.

These synthetic inputs specify interval semantics, not a particular caller's
missing-call output. The separate YT*02 heterozygote supplies the called pair;
missing YTLI evidence alone cannot rule out its reference partner.
"""

import gzip
import tempfile
import unittest
from pathlib import Path

from rbceq2.IO.vcf import VCF, Interval, build_intervals, read_vcf
from rbceq2.core_logic.constants import AlleleState
from rbceq2.core_logic.utils import Zygosity, sub_alleles_relationships
from rbceq2.db.db import Db, build_antigen_map_for_checks, prepare_db
from rbceq2.main import find_hits, parse_args


class TestGvcfBlocks(unittest.TestCase):
    """Retain covered database loci through the production file reader."""

    @classmethod
    def setUpClass(cls):
        """Load curated coordinates and prepare owned temporary inputs."""
        directory = tempfile.TemporaryDirectory(prefix="rbceq2_gvcf_blocks_")
        cls.addClassCleanup(directory.cleanup)
        cls.directory = Path(directory.name)
        cls.contexts = {}
        cls.baselines = {}
        for build in ("GRCh37", "GRCh38"):
            db = Db(ref=build, df=prepare_db())
            kn = {
                "KN": [a for a in db.make_alleles() if a.blood_group == "KN"]
            }
            relationships, group = sub_alleles_relationships(kn, "KN")
            lane = db.df.loc[db.df["Genotype"] == "YT*02"].iloc[0]
            ytli = db.df.loc[db.df["Genotype"] == "YT*01.-04"].iloc[0]
            target = next(token.strip() for token in ytli[build].split(",")
                          if not token.strip().endswith("_ref"))
            cls.contexts[build] = {
                "db": db, "relationships": {group: relationships},
                "ant_mapping": build_antigen_map_for_checks(db.df),
                "chrom": lane["Chrom"].removeprefix("chr"),
                "lane": lane[build].strip().split("_"),
                "target": target.split("_"),
            }

    def write_input(self, build, compressed, name, gt=None, offsets=None):
        """Write sorted, nonoverlapping rows with one called YT lane SNP."""
        context = self.contexts[build]
        chrom = context["chrom"]
        position, ref, alt = context["lane"]
        rows = [(int(position), f"{chrom}\t{position}\t.\t{ref}\t{alt}"
                 "\t50\tPASS\t.\tGT\t0/1\n")]
        if gt is not None:
            target, target_ref, _ = context["target"]
            start, end = (int(target) + offset for offset in offsets)
            self.assertTrue(int(position) < start or int(position) > end)
            # N is an unspecified anchor for these synthetic interval rows.
            block_ref = target_ref if start == int(target) else "N"
            rows.append((start, f"{chrom}\t{start}\t.\t{block_ref}"
                         f"\t<NON_REF>\t50\tPASS\tEND={end}\tGT\t{gt}\n"))
        content = (
            "##fileformat=VCFv4.2\n"
            '##INFO=<ID=END,Number=1,Type=Integer,'
            'Description="End position">\n'
            '##ALT=<ID=NON_REF,Description="Unspecified alternate allele">\n'
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
            + "".join(row for _, row in sorted(rows))
        )
        path = self.directory / (
            f"{build}_{name}.g.vcf{'.gz' if compressed else ''}"
        )
        opener = gzip.open if compressed else open
        with opener(path, "wt") as handle:
            handle.write(content)
        return path

    def read_input(self, build, path, intervals=None):
        """Apply both production reader filters before constructing the VCF."""
        db = self.contexts[build]["db"]
        if intervals is None:
            intervals = build_intervals(db.df, build)
        frame = read_vcf(str(path), intervals, db.unique_variants)
        return VCF(
            frame, db.lane_variants, db.unique_variants, path.stem, build
        )

    def infer(self, build, path):
        """Return reported YT calls, candidate pairs, and retained evidence."""
        context = self.contexts[build]
        _, genotype, numeric, alphanumeric, groups, _ = find_hits(
            context["db"], path,
            args=parse_args([
                "--reference_genome", build, "--HPAs", "--debug"
            ]),
            allele_relationships=context["relationships"],
            excluded=["RHD", "RHCE"], ant_mapping=context["ant_mapping"],
        )
        group = groups["YT"]
        pairs = ["/".join(pair.genotypes)
                 for pair in group.alleles[AlleleState.NORMAL] or []]
        evidence = {**group.variant_pool, **group.unused_pool}
        return (
            (genotype["YT"], numeric["YT"], alphanumeric["YT"]),
            pairs, evidence
        )

    def check_missing(self, offsets, covered, name):
        """Check calls and missing evidence without invented copy counts."""
        for build, context in self.contexts.items():
            chrom = context["chrom"]
            target = int(context["target"][0])
            prefix = f"{chrom}:{target}_"
            for compressed in (False, True):
                with self.subTest(build=build, compressed=compressed):
                    key = (build, compressed)
                    if key not in self.baselines:
                        baseline_path = self.write_input(
                            build, compressed, "ordinary"
                        )
                        self.baselines[key] = self.infer(build, baseline_path)
                    baseline, _, baseline_evidence = self.baselines[key]
                    self.assertFalse(any(
                        token.startswith(prefix) and state == Zygosity.NO_DATA
                        for token, state in baseline_evidence.items()
                    ))
                    path = self.write_input(
                        build, compressed, name, "./.", offsets
                    )
                    calls, pairs, evidence = self.infer(build, path)
                    self.assertEqual(calls[0], "YT*01/YT*02")
                    self.assertCountEqual(pairs, ["YT*01/YT*02"])
                    self.assertEqual(calls, baseline)
                    vcf = self.read_input(build, path)
                    self.assertNotIn(target, vcf.haploid_loci.get(chrom, ()))
                    self.assertNotIn(target, vcf.diploid_loci.get(chrom, ()))
                    missing = any(
                        token.startswith(prefix) and state == Zygosity.NO_DATA
                        for token, state in evidence.items()
                    )
                    self.assertEqual(
                        missing, covered,
                        f"Missing-call evidence at {chrom}:{target}"
                    )

    def test_no_call_point_preserves_evidence_and_calls(self):
        """A one-base block supplies a passing boundary control."""
        self.check_missing((0, 0), True, "missing_point")

    def test_no_call_block_covers_interior(self):
        """An interior database locus retains missing-call evidence."""
        self.check_missing((-3, 3), True, "missing_interior")

    def test_no_call_block_covers_inclusive_end(self):
        """END includes the final covered database base."""
        self.check_missing((-3, 0), True, "missing_end")

    def test_no_call_block_does_not_cover_adjacent_loci(self):
        """Blocks outside YTLI supply no evidence there."""
        self.check_missing((-3, -1), False, "missing_before")
        self.check_missing((1, 3), False, "missing_after")

    def test_reference_point_preserves_observed_ploidy(self):
        """A reference call at the locus supplies a passing ploidy control."""
        self.check_reference_ploidy((0, 0), "reference_point")

    def test_reference_block_preserves_observed_ploidy(self):
        """Retain observed diploid and haploid GTs at covered database loci."""
        for offsets, name in (((-3, 3), "reference_interior"),
                              ((-3, 0), "reference_end")):
            self.check_reference_ploidy(offsets, name)

    def check_reference_ploidy(self, offsets, name):
        """Inspect raw locus observations without interpreting dosage."""
        for build, context in self.contexts.items():
            chrom = context["chrom"]
            target = int(context["target"][0])
            for compressed in (False, True):
                for gt in ("0/0", "0"):
                    with self.subTest(build=build, compressed=compressed,
                                      gt=gt, offsets=offsets):
                        path = self.write_input(
                            build, compressed, f"{name}_{len(gt)}", gt, offsets
                        )
                        vcf = self.read_input(build, path)
                        observed = (target in vcf.haploid_loci.get(chrom, ()),
                                    target in vcf.diploid_loci.get(chrom, ()))
                        self.assertEqual(observed, (gt == "0", gt == "0/0"))

    def test_block_crossing_reader_interval_is_retained(self):
        """Retain a block whose END reaches the reader interval."""
        for build, context in self.contexts.items():
            target = int(context["target"][0])
            lane = int(context["lane"][0])
            chrom = context["chrom"]
            intervals = {
                chrom: [Interval(lane, lane), Interval(target, target)]
            }
            for compressed in (False, True):
                with self.subTest(build=build, compressed=compressed):
                    path = self.write_input(
                        build, compressed, "crossing", "0/0", (-3, 0)
                    )
                    # Isolate the coarse filter from the fine locus filter.
                    frame = read_vcf(str(path), intervals)
                    self.assertGreater(frame.height, 1)


if __name__ == "__main__":
    unittest.main()
