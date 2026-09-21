"""Preserve YT inference when an uncalled gVCF ALT placeholder is present.

These minimal synthetic variant rows are not caller-produced fixtures. They
exercise known called alleles and guard against synthesising reference evidence
from a selected placeholder or partial no-call. Reference blocks are out of
scope.
"""

import gzip
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from rbceq2.IO.vcf import VCF

from rbceq2.core_logic.constants import AlleleState
from rbceq2.core_logic.utils import sub_alleles_relationships
from rbceq2.db.db import Db, build_antigen_map_for_checks, prepare_db
from rbceq2.main import find_hits, parse_args


class TestGvcfUncalledPlaceholder(unittest.TestCase):
    """Keep the curated YT reference partner with unselected placeholders."""

    @classmethod
    def setUpClass(cls):
        """Load both builds and register cleanup for owned inputs."""
        directory = tempfile.TemporaryDirectory(prefix="rbceq2_gvcf_")
        cls.addClassCleanup(directory.cleanup)
        cls.directory = Path(directory.name)
        cls.contexts = {}
        for build in ("GRCh37", "GRCh38"):
            db = Db(ref=build, df=prepare_db())
            kn = {
                "KN": [a for a in db.make_alleles() if a.blood_group == "KN"]
            }
            relationships, group = sub_alleles_relationships(kn, "KN")
            cls.contexts[build] = (
                db, {group: relationships}, build_antigen_map_for_checks(db.df)
            )

    def infer(self, build, compressed, gt, placeholder=None,
              placeholder_first=False, extra_no_call=False):
        """Read plain or gzip variants through the production input path."""
        db, relationships, ant_mapping = self.contexts[build]
        row = db.df.loc[db.df["Genotype"] == "YT*02"].iloc[0]
        position, ref, alt = row[build].strip().split("_")
        chrom = row["Chrom"]
        if placeholder is not None:
            if placeholder_first:
                alt = f"{placeholder},{alt}"
                gt = gt.replace("1", "2")
            else:
                alt = f"{alt},{placeholder}"
        content = (
            "##fileformat=VCFv4.2\n"
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
            f"{chrom}\t{position}\t.\t{ref}\t{alt}\t50\tPASS\t.\tGT\t{gt}\n"
        )
        if extra_no_call:
            other = db.df.loc[db.df["Genotype"] == "YT*01.-04"].iloc[0]
            other_token = next(
                token.strip() for token in other[build].split(",")
                if not token.strip().endswith("_ref")
            )
            other_position, other_ref, _ = other_token.split("_")
            content += (
                f"{chrom}\t{other_position}\t.\t{other_ref}\t<NON_REF>"
                "\t50\tPASS"
                "\t.\tGT\t./.\n"
            )
        suffix = ".g.vcf" if placeholder else ".vcf"
        kind = {
            None: "ordinary", "<NON_REF>": "non_ref", "<*>": "star"
        }[placeholder]
        name = f"{build}_{kind}_{gt.replace('/', '_')}_{placeholder_first}"
        if extra_no_call:
            name += "_extra_no_call"
        path = self.directory / (name + suffix + (".gz" if compressed else ""))
        opener = gzip.open if compressed else open
        with opener(path, "wt") as handle:
            handle.write(content)
        result = find_hits(
            db, path, args=parse_args(["--reference_genome", build, "--HPAs"]),
            allele_relationships=relationships, excluded=["RHD", "RHCE"],
            ant_mapping=ant_mapping,
        )
        _, genotypes, numeric, alphanumeric, groups, _ = result
        pairs = ["/".join(pair.genotypes)
                 for pair in groups["YT"].alleles[AlleleState.NORMAL] or []]
        return genotypes["YT"], numeric["YT"], alphanumeric["YT"], pairs

    def check_calls(self, gt, placeholder=None):
        """Assert explicit ordinary results and placeholder equivalence."""
        expected_genotype = "YT*01/YT*02" if gt == "0/1" else "YT*02/YT*02"
        for build in self.contexts:
            for compressed in (False, True):
                orders = (False, True) if placeholder else (False,)
                for placeholder_first in orders:
                    with self.subTest(build=build, compressed=compressed,
                                      placeholder_first=placeholder_first):
                        actual = self.infer(
                            build, compressed, gt, placeholder,
                            placeholder_first
                        )
                        self.assertEqual(actual[0], expected_genotype)
                        self.assertCountEqual(actual[3], [expected_genotype])
                        if placeholder:
                            ordinary = self.infer(build, compressed, gt)
                            self.assertEqual(actual[:3], ordinary[:3])

    def test_ordinary_heterozygote_has_reference_partner(self):
        self.check_calls("0/1")

    def test_ordinary_homozygote_has_two_variant_copies(self):
        self.check_calls("1/1")

    def test_non_ref_homozygote_matches_ordinary(self):
        self.check_calls("1/1", "<NON_REF>")

    def test_star_homozygote_matches_ordinary(self):
        self.check_calls("1/1", "<*>")

    def test_non_ref_heterozygote_keeps_reference_partner(self):
        self.check_calls("0/1", "<NON_REF>")

    def test_star_heterozygote_keeps_reference_partner(self):
        self.check_calls("0/1", "<*>")

    def lane_vcf(self, alt_template, gt, build="GRCh38"):
        """Encode YT evidence while checking raw ALT and GT preservation."""
        db, _, _ = self.contexts[build]
        row = db.df.loc[db.df["Genotype"] == "YT*02"].iloc[0]
        position, ref, alt = row[build].strip().split("_")
        alt_field = alt_template.format(alt=alt)
        frame = pd.DataFrame(
            [[row["Chrom"], position, ".", ref, alt_field, "50", "PASS", ".",
              "GT:PS", f"{gt}:123"]],
            columns=["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER",
                     "INFO", "FORMAT", "SAMPLE"],
        )
        vcf = VCF(
            [frame], db.lane_variants, db.unique_variants, "sample", build
        )
        original = vcf.df.loc[vcf.df["ALT"] == alt_field]
        self.assertEqual(len(original), 1)
        self.assertEqual(original.iloc[0]["SAMPLE"], f"{gt}:123")
        locus = f"{row['Chrom'].removeprefix('chr')}:{position}"
        return vcf, locus, f"{locus}_{ref}_{alt}"

    def test_phased_placeholder_keeps_complementary_reference_gt_and_ps(self):
        """Keep reference evidence on the opposite haplotype in the same PS."""
        for build in self.contexts:
            for placeholder in ("<NON_REF>", "<*>"):
                for first in (False, True):
                    template = (placeholder + ",{alt}" if first
                                else "{alt}," + placeholder)
                    index = "2" if first else "1"
                    for gt, called_gt, ref_gt in (
                        (f"0|{index}", "0|1", "1|0"),
                        (f"{index}|0", "1|0", "0|1"),
                    ):
                        with self.subTest(build=build, placeholder=placeholder,
                                          first=first, gt=gt):
                            vcf, locus, token = self.lane_vcf(
                                template, gt, build
                            )
                            self.assertEqual(vcf.variants[token],
                                             {"GT": called_gt, "PS": "123"})
                            self.assertEqual(vcf.variants[f"{locus}_ref"],
                                             {"GT": ref_gt, "PS": "123"})

    def test_no_call_elsewhere_preserves_called_reference_partner(self):
        """A placeholder no-call at YTLI cannot erase the observed YT pair."""
        for build in self.contexts:
            with self.subTest(build=build):
                actual = self.infer(build, False, "0/1", "<NON_REF>",
                                    extra_no_call=True)
                self.assertEqual(actual[0], "YT*01/YT*02")
                self.assertCountEqual(actual[3], ["YT*01/YT*02"])
                self.assertEqual(
                    actual[:3], self.infer(build, False, "0/1")[:3]
                )

    def test_selected_placeholder_does_not_synthesise_observed_reference(self):
        """A selected unspecified allele cannot supply the lane SNP."""
        for placeholder in ("<NON_REF>", "<*>"):
            for first in (False, True):
                template = (placeholder + ",{alt}" if first
                            else "{alt}," + placeholder)
                placeholder_index, real_index = (
                    ("1", "2") if first else ("2", "1")
                )
                for gt in (f"0/{placeholder_index}",
                           f"{real_index}/{placeholder_index}"):
                    with self.subTest(
                        placeholder=placeholder, first=first, gt=gt
                    ):
                        vcf, locus, _ = self.lane_vcf(template, gt)
                        self.assertNotIn(f"{locus}_ref", vcf.variants)

    def test_partial_no_call_does_not_synthesise_heterozygous_reference(self):
        """An unknown chromosome copy cannot supply observed reference."""
        for template, index in (("{alt},<NON_REF>", "1"),
                                ("<NON_REF>,{alt}", "2")):
            for gt in (f"./{index}", f"{index}|.", "0/."):
                with self.subTest(template=template, gt=gt):
                    vcf, locus, _ = self.lane_vcf(template, gt)
                    self.assertNotIn(f"{locus}_ref", vcf.variants)

    def test_ordinary_multi_alt_mapping_is_unchanged(self):
        """Keep multi-ALT token dosage without adding a lane partner."""
        for gt, expected_real, expected_other in (
            ("0/1", "0/1", None),
            ("1/2", "1/0", "0/1"),
            ("2/2", None, "1/1"),
        ):
            with self.subTest(gt=gt):
                vcf, locus, token = self.lane_vcf("{alt},C", gt)
                self.assertNotIn(f"{locus}_ref", vcf.variants)
                other_token = token.rsplit("_", 1)[0] + "_C"
                for variant, expected_gt in ((token, expected_real),
                                             (other_token, expected_other)):
                    if expected_gt is None:
                        self.assertNotIn(variant, vcf.variants)
                    else:
                        self.assertEqual(vcf.variants[variant],
                                         {"GT": expected_gt, "PS": "123"})


if __name__ == "__main__":
    unittest.main()
