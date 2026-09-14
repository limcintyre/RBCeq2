"""Keep the curated FY*01N.09 sequence token usable in both genome builds."""
import unittest

import pandas as pd

from rbceq2.core_logic.constants import AlleleState
from rbceq2.core_logic.large_variants import parse_db_token
from rbceq2.core_logic.utils import sub_alleles_relationships
from rbceq2.db.db import Db, build_antigen_map_for_checks, prepare_db
from rbceq2.main import find_hits, parse_args


class TestFyNullToken(unittest.TestCase):
    """Exercise the existing broad deletion representation without redefining it."""

    @classmethod
    def setUpClass(cls):
        cls.contexts = {}
        for build in ('GRCh37', 'GRCh38'):
            db = Db(ref=build, df=prepare_db())
            kn = {'KN': [a for a in db.make_alleles() if a.blood_group == 'KN']}
            relationships, bg_type = sub_alleles_relationships(kn, 'KN')
            cls.contexts[build] = (
                db, {bg_type: relationships}, build_antigen_map_for_checks(db.df)
            )

    def test_curated_sequence_tokens_parse_in_both_builds(self):
        db = self.contexts['GRCh38'][0]
        row = db.df.loc[db.df['Genotype'] == 'FY*01N.09'].iloc[0]
        for build, position in [('GRCh37', 159175525), ('GRCh38', 159205735)]:
            with self.subTest(build=build):
                token = row[build].split(',')[0]
                definition = parse_db_token('1', token)
                self.assertIsNotNone(definition)
                self.assertEqual(definition.pos, position)
                self.assertEqual(definition.svtype, 'DEL')
                self.assertEqual(definition.length, 192)

    def assert_fy_result(self, build, deletion_gt):
        """Use the existing GRCh37 sequence at each build's curated position."""
        db, relationships, mapping = self.contexts[build]
        row = db.df.loc[db.df['Genotype'] == 'FY*01N.09'].iloc[0]
        _, ref, alt = row['GRCh37'].split(',')[0].split('_')
        position = row[build].split('_')[0]
        lane = db.df.loc[db.df['Genotype'] == 'FY*02', build].iloc[0]
        lane_pos, lane_ref, lane_alt = lane.split('_')
        frame = pd.DataFrame([
            ['chr1', lane_pos, '.', lane_ref, lane_alt, '50', 'PASS', '.', 'GT', '0/0'],
            ['chr1', position, '.', ref, alt, '50', 'PASS', '.', 'GT', deletion_gt],
        ], columns=['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'SAMPLE'])
        result = find_hits(
            db, (frame, f'FY_{build}_{deletion_gt}'),
            args=parse_args(['--reference_genome', build, '--HPAs']),
            allele_relationships=relationships, excluded=['RHD', 'RHCE'],
            ant_mapping=mapping,
        )
        _, geno, numeric, alpha, groups, _ = result
        is_deletion = deletion_gt == '1/1'
        expected = 'FY*01N.09/FY*01N.09' if is_deletion else 'FY*01/FY*01'
        with self.subTest(output='internal pairs'):
            self.assertCountEqual(
                ['/'.join(pair.genotypes) for pair in groups['FY'].alleles[AlleleState.NORMAL] or []],
                [expected],
            )
        with self.subTest(output='genotype'):
            self.assertEqual(geno['FY'], expected)
        with self.subTest(output='numeric phenotype'):
            self.assertEqual(numeric['FY'], 'FY:-1,-2,-3,-5,-6' if is_deletion else 'FY:1,-2')
        with self.subTest(output='alphanumeric phenotype'):
            self.assertEqual(alpha['FY'], 'Fy(a-b-),Fy3-,Fy5-,Fy6-' if is_deletion else 'Fy(a+b-)')

    def test_grch37_deletion(self):
        self.assert_fy_result('GRCh37', '1/1')

    def test_grch38_deletion(self):
        self.assert_fy_result('GRCh38', '1/1')

    def test_grch37_reference_control(self):
        self.assert_fy_result('GRCh37', '0/0')

    def test_grch38_reference_control(self):
        self.assert_fy_result('GRCh38', '0/0')
