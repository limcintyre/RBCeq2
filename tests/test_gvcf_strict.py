"""Specify optional strict reference fallback without fabricating called data.

Both builds use curated allele definitions and run through the production file
reader with owned temporary input files.
"""

import contextlib
import io
import tempfile
import unittest
from pathlib import Path

from rbceq2.IO.vcf import VCF, build_intervals, read_vcf
from rbceq2.core_logic.constants import AlleleState
from rbceq2.core_logic.data_procesing import get_ref
from rbceq2.core_logic.utils import Zygosity, sub_alleles_relationships
from rbceq2.db.db import Db, build_antigen_map_for_checks, prepare_db
from rbceq2.main import find_hits, parse_args


REASON = "reference_fallback_blocked_by_no_data"


class TestGvcfStrict(unittest.TestCase):
    """Exercise reference fallback and independent called alternate
    evidence.
    """

    @classmethod
    def setUpClass(cls):
        """Derive fixtures from the database and own generated input
        storage.
        """
        directory = tempfile.TemporaryDirectory(prefix="rbceq2_gvcf_strict_")
        cls.addClassCleanup(directory.cleanup)
        cls.directory = Path(directory.name)
        cls.contexts = {}
        cls.results = {}
        for build in ("GRCh37", "GRCh38"):
            db = Db(ref=build, df=prepare_db())
            rows = db.df.set_index("Genotype")

            def tokens(name):
                return {value.strip() for value in
                        rows.loc[name, build].split(",")}
            crom, = tokens("CROM*01")
            lu, = tokens("LU*02")
            kn, = {token for token in tokens("KN*01")
                   if token.endswith("_ref")}
            kn_shared, = tokens("KN*01") - {kn}
            kn_alt, = tokens("KN*01.-05") - {kn_shared}
            modifier, = tokens("LU*02.19") - {lu}
            unrelated, = tokens("CROM*01.-05") - {crom}
            reference_bases = {}
            for reference in (crom, lu, kn):
                position = reference.split("_")[0]
                bases = {
                    parts[1] for definition in rows[build].dropna()
                    for token in definition.split(",")
                    if len(parts := token.strip().split("_")) == 3
                    and parts[0] == position
                    and parts[1] in {"A", "C", "G", "T"}
                    and parts[2] in {"A", "C", "G", "T"}
                }
                if len(bases) != 1:
                    raise AssertionError("Reference locus has ambiguous bases")
                reference_bases[position] = bases.pop()
            kn_alleles = {"KN": [a for a in db.make_alleles()
                                 if a.blood_group == "KN"]}
            relations, group = sub_alleles_relationships(kn_alleles, "KN")
            cls.contexts[build] = {
                "db": db, "mapping": build_antigen_map_for_checks(db.df),
                "relations": {group: relations},
                "reference_bases": reference_bases,
                "targets": [("1", crom), ("19", lu), ("1", kn)],
                "kn_calls": [("1", kn_alt, "0/1"),
                             ("1", kn_shared, "1/1")],
                "modifier": ("19", modifier),
                "unrelated": ("1", unrelated),
            }

    def input_path(self, build, kind, gt, modifier=False, unrelated=False):
        """Create independent rows with optional missing blocks or Au calls."""
        context = self.contexts[build]
        rows = []
        for chrom, token, call in context["kn_calls"]:
            pos, ref, alt = token.split("_")
            rows.append([chrom, pos, ".", ref, alt, "50", "PASS",
                         ".", "GT", call])
        if modifier:
            chrom, token = context["modifier"]
            pos, ref, alt = token.split("_")
            rows.append([chrom, pos, ".", ref, alt, "50", "PASS",
                         ".", "GT", "1/1"])
        if kind != "absent":
            for chrom, token in context["targets"]:
                pos = int(token.split("_")[0])
                start = pos - 3 if kind == "block" else pos
                end = pos + 3 if kind == "block" else pos
                ref = ("N" if kind == "block" else
                       context["reference_bases"][str(pos)])
                rows.append([chrom, str(start), ".", ref, "<NON_REF>",
                             "50", "PASS", f"END={end}", "GT", gt])
        if unrelated:
            chrom, token = context["unrelated"]
            pos, ref, alt = token.split("_")
            rows.append([chrom, pos, ".", ref, alt, "50", "PASS",
                         ".", "GT", "./."])
        path = self.directory / (
            f"{build}_{kind}_{gt.replace('/', '_')}_{modifier}_{unrelated}.vcf"
        )
        if not path.exists():
            path.write_text(
                "##fileformat=VCFv4.2\n"
                '##INFO=<ID=END,Number=1,Type=Integer,Description="End">\n'
                '##ALT=<ID=NON_REF,'
                'Description="Unspecified alternate allele">\n'
                '##FORMAT=<ID=GT,Number=1,Type=String,'
                'Description="Genotype">\n'
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
                "\tFORMAT\tSAMPLE\n"
                + "".join("\t".join(row) + "\n" for row in
                          sorted(rows, key=lambda r: (int(r[0]), int(r[1]))))
            )
        return path

    def infer(self, build, kind, strict, gt="./.", modifier=False,
              unrelated=False, no_filter=False):
        """Use the production pipeline with the proposed option set
        explicitly.
        """
        key = (build, kind, strict, gt, modifier, unrelated, no_filter)
        if key not in self.results:
            context = self.contexts[build]
            path = self.input_path(build, kind, gt, modifier, unrelated)
            argv = ["--reference_genome", build, "--HPAs", "--debug"]
            if no_filter:
                argv.append("--no_filter")
            args = parse_args(argv)
            # Keep behavior failures distinct from the separate CLI contract.
            args.GVCF_strict = strict
            original = path.read_bytes()
            result = find_hits(
                context["db"], path, args=args,
                allele_relationships=context["relations"],
                excluded=["RHD", "RHCE"], ant_mapping=context["mapping"],
            )
            self.assertEqual(path.read_bytes(), original)
            self.results[key] = result
        return self.results[key]

    def assert_calls(self, result, expected):
        """Check output genotype multiplicity, not display ordering."""
        for group, alleles in expected.items():
            with self.subTest(group=group):
                self.assertCountEqual(
                    result[1].get(group, "").split("/"), alleles
                )

    def assert_blocked(self, result, groups=("CROM", "LU", "KN")):
        """Require named exclusions and prevent later pair restoration."""
        for name in groups:
            with self.subTest(group=name):
                with self.subTest(check="phenotypes"):
                    self.assertFalse(result[2].get(name))
                    self.assertFalse(result[3].get(name))
                self.assertIn(name, result[4])
                group = result[4][name]
                removed = group.filtered_out.get(REASON, [])
                with self.subTest(check="named exclusion"):
                    self.assertTrue(removed, f"{name} has no named exclusion")
                with self.subTest(check="sample identity"):
                    self.assertEqual(group.sample, result[0])
                final = group.alleles[AlleleState.NORMAL] or []
                removed_pairs = [a for a in removed if hasattr(a, "genotypes")]
                for pair in removed_pairs:
                    self.assertNotIn(pair, final)
                    self.assertNotIn(sorted(pair.genotypes),
                                     [sorted(p.genotypes) for p in final])

    def test_cli_contract(self):
        """The flag is opt-in and parsed with the documented capitalization."""
        with self.subTest(mode="default"):
            self.assertIs(getattr(parse_args(["--reference_genome", "GRCh38"]),
                                  "GVCF_strict", None), False)
        with self.subTest(mode="strict"):
            with contextlib.redirect_stderr(io.StringIO()):
                try:
                    args = parse_args(["--reference_genome", "GRCh38",
                                       "--GVCF_strict"])
                except SystemExit as error:
                    self.fail(f"--GVCF_strict was rejected: {error}")
            self.assertIs(args.GVCF_strict, True)

    def test_no_call_tokens_never_become_called_reference(self):
        """Raw missing observations cannot acquire synthetic reference
        calls.
        """
        for build, context in self.contexts.items():
            for kind in ("point", "block"):
                with self.subTest(build=build, kind=kind):
                    path = self.input_path(build, kind, "./.")
                    db = context["db"]
                    frames = read_vcf(str(path), build_intervals(db.df, build),
                                      db.unique_variants)
                    vcf = VCF(frames, db.lane_variants, db.unique_variants,
                              path.stem, build)
                    for chrom, token in context["targets"]:
                        pos = int(token.split("_")[0])
                        locus = f"{chrom}:{pos}"
                        observed = {
                            key: metrics["GT"]
                            for key, metrics in vcf.variants.items()
                            if key.startswith(locus + "_")
                        }
                        with self.subTest(locus=locus):
                            self.assertTrue(observed)
                            self.assertTrue(
                                any(gt in (".", "./.", ".|.")
                                    for gt in observed.values()), observed
                            )
                            for key, gt in observed.items():
                                with self.subTest(token=key):
                                    self.assertEqual(get_ref({"GT": gt}, key),
                                                     Zygosity.NO_DATA)
                            self.assertNotIn(
                                pos, vcf.haploid_loci.get(chrom, ())
                            )
                            self.assertNotIn(
                                pos, vcf.diploid_loci.get(chrom, ())
                            )

    def test_missing_reference_default_fallback(self):
        """Default mode may report reference alleles despite missing
        evidence.
        """
        for build in self.contexts:
            for kind in ("point", "block"):
                with self.subTest(build=build, kind=kind):
                    self.assert_calls(self.infer(build, kind, False), {
                        "CROM": ["CROM*01", "CROM*01"],
                        "LU": ["LU*02", "LU*02"],
                        "KN": ["KN*01", "KN*01.-05"],
                    })

    def test_missing_reference_strict(self):
        """Strict mode keeps the independently evidenced KN alternate."""
        for build in self.contexts:
            for kind in ("point", "block"):
                with self.subTest(build=build, kind=kind):
                    result = self.infer(build, kind, True)
                    self.assert_calls(result, {
                        "CROM": ["Undetermined", "Undetermined"],
                        "LU": ["Undetermined", "Undetermined"],
                        "KN": ["KN*01.-05", "Undetermined"],
                    })
                    self.assert_blocked(result)

    def test_called_reference_and_absent_row_controls(self):
        """Strict mode leaves called reference and legacy absent rows
        unchanged.
        """
        for build in self.contexts:
            for kind in ("point", "block", "absent"):
                with self.subTest(build=build, kind=kind):
                    default = self.infer(build, kind, False, gt="0/0")
                    strict = self.infer(build, kind, True, gt="0/0")
                    self.assert_calls(strict, {
                        "CROM": ["CROM*01", "CROM*01"],
                        "LU": ["LU*02", "LU*02"],
                        "KN": ["KN*01", "KN*01.-05"],
                    })
                    self.assertEqual(default[1:4], strict[1:4])

    def test_unrelated_missing_row_does_not_block_reference(self):
        """An unknown CROM modifier does not erase its called reference
        locus.
        """
        for build in self.contexts:
            with self.subTest(build=build):
                result = self.infer(
                    build, "point", True, gt="0/0", unrelated=True
                )
                self.assert_calls(result, {"CROM": ["CROM*01", "CROM*01"]})

    def test_called_au_modifier_needs_reference_background(self):
        """Au evidence cannot confirm LU*02.19 when its background is
        unknown.
        """
        for build in self.contexts:
            for kind in ("point", "block"):
                with self.subTest(build=build, kind=kind):
                    result = self.infer(build, kind, True, modifier=True)
                    self.assert_calls(
                        result, {"LU": ["Undetermined", "Undetermined"]}
                    )
                    self.assert_blocked(result, ("LU",))

    def test_no_filter_does_not_disable_strict_evidence_rule(self):
        """Disabling ordinary filters must not license missing reference
        calls.
        """
        result = self.infer("GRCh38", "block", True, no_filter=True)
        self.assert_calls(result, {"CROM": ["Undetermined", "Undetermined"],
                                   "LU": ["Undetermined", "Undetermined"]})
        self.assert_blocked(result, ("CROM", "LU"))


if __name__ == "__main__":
    unittest.main()
