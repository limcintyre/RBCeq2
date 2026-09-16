import unittest
from unittest.mock import patch

from rbceq2.core_logic.filter_semantics import (
    FilterMeaning,
    filter_excludes_allele,
    filter_value_meanings,
)

TABLE = (
    "Value\tClassification\tOrigin\tEvidence\n"
    "LowQual\tsample_call\tlong_read\tLow quality variant\n"
    "LowMLSQ\tsite_statistic\tshort_read\tSite with ML Site Quality score lower than 0.1\n"
    "DUP\tnot_correctness\tarray\tAssay for duplicate site defined on the chip\n"
)


def with_table(table: str):
    """Patch the resource loader and clear the cache around one test."""
    filter_value_meanings.cache_clear()
    return patch(
        "rbceq2.core_logic.filter_semantics.load_filter_values", return_value=table
    )


class TestFilterValueMeanings(unittest.TestCase):
    def tearDown(self):
        filter_value_meanings.cache_clear()

    def test_table_is_read_into_a_mapping(self):
        with with_table(TABLE):
            self.assertEqual(
                filter_value_meanings(),
                {
                    "LowQual": FilterMeaning.SAMPLE_CALL,
                    "LowMLSQ": FilterMeaning.SITE_STATISTIC,
                    "DUP": FilterMeaning.NOT_CORRECTNESS,
                },
            )

    def test_unknown_classification_raises_on_load(self):
        """A typo would silently change filtering for every sample, so it is caught here."""
        bad = "Value\tClassification\tOrigin\tEvidence\nLowQual\tsample_cal\tx\ty\n"
        with with_table(bad):
            with self.assertRaises(ValueError) as ctx:
                filter_value_meanings()
        self.assertIn("sample_cal", str(ctx.exception))


class TestFilterExcludesAllele(unittest.TestCase):
    def setUp(self):
        self.patcher = with_table(TABLE)
        self.patcher.start()

    def tearDown(self):
        self.patcher.stop()
        filter_value_meanings.cache_clear()

    def test_pass_never_excludes(self):
        self.assertEqual(filter_excludes_allele("PASS"), (False, ()))

    def test_absent_filter_never_excludes(self):
        """'.' is the VCF specification's 'no filtering applied', not a failure."""
        for field in (".", "", "  "):
            with self.subTest(field=field):
                self.assertEqual(filter_excludes_allele(field), (False, ()))

    def test_a_doubt_about_this_sample_excludes(self):
        self.assertEqual(filter_excludes_allele("LowQual"), (True, ()))

    def test_a_cohort_statistic_does_not_exclude(self):
        """The value behind the jointly called cohort's 256 genotype conflicts.

        A site quality score computed across every sample says nothing about whether this
        sample's call is right, so it is not grounds to drop this sample's allele.
        """
        self.assertEqual(filter_excludes_allele("LowMLSQ"), (False, ()))

    def test_probeset_selection_does_not_exclude(self):
        """On an array, PASS/FAIL can mark which probeset is recommended for a marker."""
        self.assertEqual(filter_excludes_allele("DUP"), (False, ()))

    def test_unknown_value_excludes_and_is_reported(self):
        """Fail closed - a new caller's vocabulary is not assumed benign."""
        self.assertEqual(filter_excludes_allele("NovelValue"), (True, ("NovelValue",)))

    def test_semicolon_joined_values_are_split(self):
        """A row may carry several values; 28 combinations appear in one cohort."""
        self.assertEqual(filter_excludes_allele("LowMLSQ;DUP"), (False, ()))
        self.assertEqual(filter_excludes_allele("DUP;LowQual"), (True, ()))
        self.assertEqual(
            filter_excludes_allele("DUP;NovelValue"), (True, ("NovelValue",))
        )

    def test_any_one_doubt_is_enough(self):
        self.assertEqual(filter_excludes_allele("LowQual;LowMLSQ;DUP"), (True, ()))


class TestShippedTable(unittest.TestCase):
    """Against the real resources/filter_values.tsv, not a double."""

    def tearDown(self):
        filter_value_meanings.cache_clear()

    def test_the_shipped_table_loads_and_classifies(self):
        filter_value_meanings.cache_clear()
        meanings = filter_value_meanings()
        self.assertEqual(meanings["LowQual"], FilterMeaning.SAMPLE_CALL)
        self.assertEqual(meanings["LowMLSQ"], FilterMeaning.SITE_STATISTIC)
        self.assertEqual(meanings["DUP"], FilterMeaning.NOT_CORRECTNESS)

    def test_the_structural_and_copy_number_values_are_recognised(self):
        """Ten values that were reaching the table unclassified.

        Every one is a reason to doubt this sample's call - a CNV quality score, a QUAL
        threshold, a depth or mappability anomaly at a breakend - so all ten classify as
        sample_call and go on excluding exactly as they did while unrecognised. What
        changes is the second half of the tuple: they no longer report as unclassified,
        so the warning stops firing on values somebody has already ruled on and a
        genuinely new vocabulary still stands out.

        Two of them, MaxMQ0Frac and NoPairSupport, describe a statistic over "all
        samples" and so read as site-scoped on a jointly called file and sample-scoped on
        a single-sample one. Both forms occur here and the table is keyed by value alone,
        so it cannot hold both readings; sample_call is the side that keeps the stricter
        behaviour.
        """
        filter_value_meanings.cache_clear()
        for value in (
            "Ploidy",
            "cnvLikelihoodRatio",
            "cnvQual",
            "MinQUAL",
            "MaxDepth",
            "MaxMQ0Frac",
            "NoPairSupport",
            "dinucQual",
            "DRAGENIndelHardQUAL",
            "MosaicLowAF",
        ):
            with self.subTest(value=value):
                self.assertEqual(filter_excludes_allele(value), (True, ()))

    def test_ploidy_agrees_with_the_value_beside_it(self):
        """'Ploidy' and 'PloidyConflict' describe the same kind of inconsistency.

        One is the structural caller saying overlapping genotypes cannot all be right,
        the other the small variant caller saying a genotype disagrees with the
        chromosome ploidy. They were classified years apart and must not diverge.
        """
        filter_value_meanings.cache_clear()
        meanings = filter_value_meanings()
        self.assertEqual(meanings["Ploidy"], meanings["PloidyConflict"])

    def test_long_read_behaviour_is_unchanged(self):
        """LowQual is the only non-PASS value at database loci on both long-read sets.

        It classifies as a doubt about this sample's call, so those datasets keep the
        behaviour they had before the table existed. This is what makes the change move
        zero cells on every dataset that has a gold standard.
        """
        filter_value_meanings.cache_clear()
        self.assertEqual(filter_excludes_allele("LowQual"), (True, ()))


if __name__ == "__main__":
    unittest.main()
