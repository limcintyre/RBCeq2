"""Preserve exact curated SV matches when overlap scores tie."""
import unittest

from rbceq2.core_logic.large_variants import SvEvent, SvMatcher, parse_db_token
from rbceq2.db.db import prepare_db


class TestExactSvMatchTies(unittest.TestCase):
    """An exact match must not lose to a nearby event because of input order."""

    @classmethod
    def setUpClass(cls):
        """Use the curated GRCh38 GCNT2 deletion definition."""
        frame = prepare_db()
        row = frame.loc[frame['Genotype'] == 'GCNT2*01N.06'].iloc[0]
        cls.definition = parse_db_token('6', row['GRCh38'], 'GCNT2*01N.06')

    def assert_exact_selected(self, exact_first):
        """Check that selected coordinates, source token and genotype stay linked."""
        definition = self.definition
        self.assertIsNotNone(definition)

        def event(offset, gt):
            pos = definition.pos + offset
            return SvEvent(
                chrom=definition.chrom, pos=pos, end=pos + definition.length,
                svtype='DEL', svlen=-definition.length, alt='<DEL>', id='.',
                qual='50', variant=f'{definition.chrom}:{pos}_del_41kb', info={},
                sample_fmt='GT', sample_value=gt,
            )

        exact, nearby = event(0, '0/1'), event(11, '1/1')
        events = [exact, nearby] if exact_first else [nearby, exact]
        matches = SvMatcher().match([definition], events)
        self.assertEqual(len(matches), 1)
        chosen = matches[0]
        self.assertIs(chosen.vcf, exact)
        self.assertEqual(chosen.variant, exact.variant)
        self.assertEqual(chosen.vcf.sample_value, '0/1')
        self.assertEqual((chosen.pos_delta, chosen.len_delta), (0, 0))

    def test_exact_event_first(self):
        self.assert_exact_selected(True)

    def test_exact_event_last(self):
        self.assert_exact_selected(False)


    def test_nonexact_tie_preserves_position_and_length_weighting(self):
        """Prefer lower existing geometric error without adding POS-first priority."""
        definition = self.definition

        def event(offset, extra_length, gt):
            pos = definition.pos + offset
            length = definition.length + extra_length
            return SvEvent(
                chrom=definition.chrom, pos=pos, end=pos + length,
                svtype='DEL', svlen=-length, alt='<DEL>', id='.', qual='50',
                variant=f'{definition.chrom}:{pos}_DEL_{length}', info={},
                sample_fmt='GT', sample_value=gt,
            )

        less_error = event(100, 0, '0/1')
        more_error = event(0, 1000, '1/1')
        matcher = SvMatcher()
        self.assertEqual(matcher.score(definition, less_error)[0], 0.0)
        self.assertEqual(matcher.score(definition, more_error)[0], 0.0)
        for events in ([less_error, more_error], [more_error, less_error]):
            with self.subTest(first=events[0].variant):
                matches = matcher.match([definition], events)
                self.assertEqual(len(matches), 1)
                self.assertIs(matches[0].vcf, less_error)
                self.assertEqual(matches[0].score, 0.0)
