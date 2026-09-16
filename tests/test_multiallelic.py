#!/usr/bin/env python3
"""Unit tests for multi-allelic locus support.

A row with two alternates carries one genotype describing both. The encoder already
emitted one token per alternate; what was missing was a genotype per token, so both
tokens of a '1/2' were handed '1/2' and get_ref refused the row as multi-allelic.

Two things are being pinned here and they fail differently:

    zygosity   '1/2' has to reach get_ref as '1/0' and '0/1', so each alternate is HET
    phase      the phased filters compare phase strings literally, so '1|2' has to
               become '1|0' and '0|1' or two different alternates at one position read
               as being in phase with each other - a wrong answer rather than a raise
"""

import unittest

import pandas as pd

from rbceq2.core_logic.data_procesing import get_ref
from rbceq2.core_logic.utils import BeyondLogicError, Zygosity
from rbceq2.IO.vcf import VCF, recode_gt_for_alt_index


COMMON = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]


def rows_df(rows: list[tuple[str, str, str, str, str]]) -> pd.DataFrame:
    """A VCF frame. rows are (chrom, pos, REF, ALT, SAMPLE)."""
    return pd.DataFrame(
        {
            "CHROM": [r[0] for r in rows],
            "POS": [r[1] for r in rows],
            "ID": ["."] * len(rows),
            "REF": [r[2] for r in rows],
            "ALT": [r[3] for r in rows],
            "QUAL": ["."] * len(rows),
            "FILTER": ["PASS"] * len(rows),
            "INFO": ["."] * len(rows),
            "FORMAT": ["GT"] * len(rows),
            "SAMPLE": [r[4] for r in rows],
        }
    )


def a_vcf(rows: list[tuple[str, str, str, str, str]]) -> VCF:
    return VCF(
        [rows_df(rows)],
        lane_variants={},
        unique_variants=set(),
        sample="s",
        reference_genome="GRCh38",
    )


class TestRecodeGtForAltIndex(unittest.TestCase):
    """The whole change is this function. The truth table is the specification."""

    def test_het_for_two_different_alternates(self) -> None:
        """The common multi-allelic genotype: one copy of each."""
        self.assertEqual(recode_gt_for_alt_index("1/2", 1), "1/0")
        self.assertEqual(recode_gt_for_alt_index("1/2", 2), "0/1")

    def test_hom_for_the_second_alternate(self) -> None:
        """'2/2' is two copies of the second and none of the first."""
        self.assertEqual(recode_gt_for_alt_index("2/2", 1), "0/0")
        self.assertEqual(recode_gt_for_alt_index("2/2", 2), "1/1")

    def test_het_with_reference(self) -> None:
        self.assertEqual(recode_gt_for_alt_index("0/2", 1), "0/0")
        self.assertEqual(recode_gt_for_alt_index("0/2", 2), "0/1")

    def test_the_separator_is_preserved(self) -> None:
        """A phased genotype has to stay phased or the phase is silently discarded."""
        self.assertEqual(recode_gt_for_alt_index("1|2", 1), "1|0")
        self.assertEqual(recode_gt_for_alt_index("1|2", 2), "0|1")

    def test_a_half_call_keeps_its_no_call(self) -> None:
        """'.' is the absence of a measurement and is not recoded to anything."""
        self.assertEqual(recode_gt_for_alt_index("./2", 2), "./1")
        self.assertEqual(recode_gt_for_alt_index("2/.", 1), "0/.")

    def test_haploid(self) -> None:
        self.assertEqual(recode_gt_for_alt_index("2", 2), "1")
        self.assertEqual(recode_gt_for_alt_index("2", 1), "0")

    def test_an_index_beyond_the_alternates_reads_as_absent(self) -> None:
        """Third alternate named, two tokens - each says 'not mine', which is true."""
        self.assertEqual(recode_gt_for_alt_index("1/3", 2), "0/0")

    def test_an_empty_genotype_is_left_alone(self) -> None:
        self.assertEqual(recode_gt_for_alt_index("", 1), "")


class TestGetVariantsSplitsTheGenotype(unittest.TestCase):
    """End to end through the VCF layer, which is where the recode is applied."""

    def test_each_alternate_gets_its_own_genotype(self) -> None:
        vcf = a_vcf([("chr1", "25272548", "A", "G,T", "1/2")])
        self.assertEqual(vcf.variants["1:25272548_A_G"]["GT"], "1/0")
        self.assertEqual(vcf.variants["1:25272548_A_T"]["GT"], "0/1")

    def test_an_alternate_the_sample_does_not_carry_is_absent(self) -> None:
        """Zero copies is encoded as absence, the same as a hom ref row."""
        vcf = a_vcf([("chr1", "25272548", "A", "G,T", "2/2")])
        self.assertNotIn("1:25272548_A_G", vcf.variants)
        self.assertEqual(vcf.variants["1:25272548_A_T"]["GT"], "1/1")

    def test_the_metrics_dicts_are_not_shared(self) -> None:
        """They were one object before. Nothing depended on it; nothing should."""
        vcf = a_vcf([("chr1", "25272548", "A", "G,T", "1/2")])
        self.assertIsNot(
            vcf.variants["1:25272548_A_G"], vcf.variants["1:25272548_A_T"]
        )

    def test_other_format_fields_survive_the_split(self) -> None:
        df = rows_df([("chr1", "25272548", "A", "G,T", "1|2:30:159155563")])
        df["FORMAT"] = ["GT:GQ:PS"]
        vcf = VCF(
            [df],
            lane_variants={},
            unique_variants=set(),
            sample="s",
            reference_genome="GRCh38",
        )
        for token in ("1:25272548_A_G", "1:25272548_A_T"):
            with self.subTest(token=token):
                self.assertEqual(vcf.variants[token]["GQ"], "30")
                self.assertEqual(vcf.variants[token]["PS"], "159155563")

    def test_the_phase_set_stays_shared(self) -> None:
        """Deliberate - it is one physical locus, so it is one phase set."""
        df = rows_df([("chr1", "25272548", "A", "G,T", "1|2:159155563")])
        df["FORMAT"] = ["GT:PS"]
        vcf = VCF(
            [df],
            lane_variants={},
            unique_variants=set(),
            sample="s",
            reference_genome="GRCh38",
        )
        self.assertEqual(
            vcf.variants["1:25272548_A_G"]["PS"],
            vcf.variants["1:25272548_A_T"]["PS"],
        )

    def test_a_biallelic_row_is_untouched(self) -> None:
        """The regression guard for every existing dataset."""
        vcf = a_vcf([("chr1", "25272548", "A", "G", "0/1")])
        self.assertEqual(vcf.variants["1:25272548_A_G"]["GT"], "0/1")


class TestZygosityOfAMultiAllelicLocus(unittest.TestCase):
    """What the pool ends up saying, which is the point of the whole change."""

    def test_two_different_alternates_are_both_het(self) -> None:
        vcf = a_vcf([("chr1", "25272548", "A", "G,T", "1/2")])
        self.assertEqual(get_ref(vcf.variants["1:25272548_A_G"]), Zygosity.HET)
        self.assertEqual(get_ref(vcf.variants["1:25272548_A_T"]), Zygosity.HET)

    def test_a_homozygous_second_alternate_is_hom(self) -> None:
        vcf = a_vcf([("chr1", "25272548", "A", "G,T", "2/2")])
        self.assertEqual(get_ref(vcf.variants["1:25272548_A_T"]), Zygosity.HOM)

    def test_a_half_called_alternate_is_no_data(self) -> None:
        vcf = a_vcf([("chr1", "25272548", "A", "G,T", "./2")])
        self.assertEqual(get_ref(vcf.variants["1:25272548_A_T"]), Zygosity.NO_DATA)

    def test_a_haploid_multi_allelic_locus_is_hemizygous(self) -> None:
        """The case a gene at one copy produces - state table row D2 meeting a split."""
        vcf = a_vcf([("chr1", "25272548", "A", "G,T", "2")])
        self.assertEqual(
            get_ref(vcf.variants["1:25272548_A_T"], "1:25272548_A_T", locus_copies=1),
            Zygosity.HEM,
        )

    def test_get_ref_never_sees_a_multi_allelic_genotype(self) -> None:
        """Which is why none of its logic had to change."""
        vcf = a_vcf([("chr1", "25272548", "A", "G,T", "1/2")])
        for metrics in vcf.variants.values():
            with self.subTest(GT=metrics["GT"]):
                self.assertLessEqual(set(metrics["GT"]) - {"/", "|", "."}, {"0", "1"})


class TestPhaseOfAMultiAllelicLocus(unittest.TestCase):
    """The silent half. The phased filters compare these strings literally."""

    def test_the_two_alternates_are_on_opposite_sides(self) -> None:
        vcf = a_vcf([("chr1", "25272548", "A", "G,T", "1|2")])
        self.assertEqual(vcf.variants["1:25272548_A_G"]["GT"], "1|0")
        self.assertEqual(vcf.variants["1:25272548_A_T"]["GT"], "0|1")

    def test_they_no_longer_compare_as_being_in_phase(self) -> None:
        """Both carried the literal string '1|2', so `phase == phase2` was True."""
        vcf = a_vcf([("chr1", "25272548", "A", "G,T", "1|2")])
        self.assertNotEqual(
            vcf.variants["1:25272548_A_G"]["GT"],
            vcf.variants["1:25272548_A_T"]["GT"],
        )

    def test_the_sides_are_the_ones_the_filters_test_for(self) -> None:
        """phased.py compares against these exact strings."""
        vcf = a_vcf([("chr1", "25272548", "A", "G,T", "1|2")])
        sides = {m["GT"] for m in vcf.variants.values()}
        self.assertEqual(sides, {"1|0", "0|1"})


class TestTheRaisesStayDistinguishable(unittest.TestCase):
    """A row with one alternate is still refused, and says which check refused it."""

    def test_a_single_alternate_row_with_a_high_index_still_raises(self) -> None:
        """Not recoded, deliberately - '2/2' with one alternate is malformed input."""
        with self.assertRaises(BeyondLogicError) as ctx:
            get_ref({"GT": "1/2"}, "1:25272548_A_G")
        self.assertIn("Multi-allelic", str(ctx.exception))

    def test_the_diploid_raise_names_itself(self) -> None:
        with self.assertRaises(BeyondLogicError) as ctx:
            get_ref({"GT": "1/2"}, "1:25272548_A_G")
        self.assertIn("get_ref/multi_allelic_diploid_GT", str(ctx.exception))

    def test_the_haploid_raise_names_itself_differently(self) -> None:
        """Same wording, different check. This is the pair that was indistinguishable."""
        with self.assertRaises(BeyondLogicError) as ctx:
            get_ref({"GT": "2"}, "X:37686068_G_A", chrom_copies=1)
        self.assertIn("get_ref/multi_allelic_haploid_GT", str(ctx.exception))

    def test_the_two_are_not_the_same_string(self) -> None:
        with self.assertRaises(BeyondLogicError) as diploid:
            get_ref({"GT": "1/2"}, "1:25272548_A_G")
        with self.assertRaises(BeyondLogicError) as haploid:
            get_ref({"GT": "2"}, "X:37686068_G_A", chrom_copies=1)
        self.assertNotEqual(str(diploid.exception), str(haploid.exception))

    def test_the_haploid_rejection_still_names_itself_too(self) -> None:
        """The third raise in get_ref, which shares neither wording nor cause."""
        with self.assertRaises(BeyondLogicError) as ctx:
            get_ref({"GT": "1"}, "1:25272548_G_A", chrom_copies=2)
        self.assertIn(
            "get_ref/haploid_GT_where_neither_count_is_one", str(ctx.exception)
        )

    def test_the_context_still_carries_the_variant(self) -> None:
        with self.assertRaises(BeyondLogicError) as ctx:
            get_ref({"GT": "1/2"}, "1:25272548_A_G")
        self.assertIn("1:25272548_A_G", str(ctx.exception))


class TestBeyondLogicErrorFormatting(unittest.TestCase):
    """The raise site leads the message so it is the first thing in a log line."""

    def test_a_named_error_leads_with_the_name(self) -> None:
        err = BeyondLogicError(message="boom", raised_by="somewhere/some_reason")
        self.assertEqual(str(err), "[somewhere/some_reason] boom")

    def test_name_and_context_together(self) -> None:
        err = BeyondLogicError(message="boom", context="x: 1", raised_by="a/b")
        self.assertEqual(str(err), "[a/b] boom | Context: x: 1")

    def test_an_unnamed_error_is_unchanged(self) -> None:
        """Every pre-existing raise keeps its exact wording."""
        self.assertEqual(str(BeyondLogicError(message="boom")), "boom")
        self.assertEqual(
            str(BeyondLogicError(message="boom", context="x: 1")),
            "boom | Context: x: 1",
        )


if __name__ == "__main__":
    unittest.main()
