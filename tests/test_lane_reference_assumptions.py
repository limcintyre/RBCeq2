"""Keep missing-lane reference assumptions distinct from observed calls."""
import unittest

import pandas as pd

from rbceq2.IO.vcf import VCF
from rbceq2.core_logic.constants import SYNTHESISED_HOM_REF_GT
from rbceq2.db.db import Db, prepare_db
from rbceq2.filters.shared_filter_functionality import carries_phase


class TestLaneReferenceAssumptions(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        """Load a curated lane locus for each build."""
        cls.contexts = {}
        for build in ('GRCh37', 'GRCh38'):
            db = Db(ref=build, df=prepare_db())
            row = db.df.loc[db.df['Genotype'] == 'LU*01'].iloc[0]
            pos, ref, alt = row[build].split('_')
            cls.contexts[build] = db, row['Chrom'], pos, ref, alt

    def read(self, build, calls, alternate=None):
        """Return a VCF and its native rows at the test locus."""
        db, chrom, pos, ref, alt = self.contexts[build]
        if alternate is not None:
            alt = alternate
        rows = [[chrom, pos, '.', ref, alt, '50', 'PASS', '.',
                 'GT:PS', f'{gt}:123'] for gt in calls]
        frame = pd.DataFrame(rows, columns=[
            'CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER',
            'INFO', 'FORMAT', 'SAMPLE',
        ])
        vcf = VCF([frame], db.lane_variants, db.unique_variants,
                  'lane_assumptions', build)
        locus = f"{chrom.removeprefix('chr')}:{pos}"
        return vcf, locus, alt

    def test_wholly_missing_gt_keeps_native_metrics_and_uses_sentinel(self):
        """Missing observations stay missing beside an explicit assumption."""
        for build in self.contexts:
            for gt in ('.', './.', '.|.', './././.', '.|.|.|.'):
                with self.subTest(build=build, gt=gt):
                    vcf, locus, alt = self.read(build, [gt])
                    self.assertIn(f'{locus}_ref', vcf.variants)
                    assumed = vcf.variants[f'{locus}_ref']['GT']
                    self.assertEqual(assumed, SYNTHESISED_HOM_REF_GT)
                    self.assertFalse(carries_phase(assumed))
                    native = vcf.df.loc[vcf.df['ALT'] == alt]
                    self.assertEqual(native['SAMPLE'].tolist(), [f'{gt}:123'])
                    chrom, pos = locus.split(':')
                    self.assertNotIn(int(pos), vcf.haploid_loci.get(chrom, ()))
                    self.assertNotIn(int(pos), vcf.diploid_loci.get(chrom, ()))

    def test_called_alt_prevents_assumption_in_either_row_order(self):
        """A missing row cannot override a measured alternate at the locus."""
        for build in self.contexts:
            for calls in (['./.', '1/1'], ['1/1', './.'],
                          ['./.', '0/.'], ['0/.', './.']):
                with self.subTest(build=build, calls=calls):
                    vcf, locus, _ = self.read(build, calls)
                    self.assertNotIn(f'{locus}_ref', vcf.variants)

    def test_partial_genotypes_do_not_receive_wholly_missing_fallback(self):
        """Partial calls retain their existing handling."""
        for build in self.contexts:
            for gt in ('0/.', './1', '1|.'):
                with self.subTest(build=build, gt=gt):
                    vcf, locus, _ = self.read(build, [gt])
                    self.assertNotIn(f'{locus}_ref', vcf.variants)

    def test_called_unspecified_allele_does_not_become_reference(self):
        """A called non-reference placeholder is not a reference assumption."""
        for build in self.contexts:
            with self.subTest(build=build):
                vcf, locus, _ = self.read(build, ['1/1'], '<NON_REF>')
                self.assertNotIn(f'{locus}_ref', vcf.variants)


if __name__ == '__main__':
    unittest.main()
