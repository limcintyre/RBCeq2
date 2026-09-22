"""Guard strict fallback across pair routes and copy-count output shapes."""
import unittest
from types import SimpleNamespace

import pandas as pd

from rbceq2.core_logic.alleles import Allele, BloodGroup, Pair
from rbceq2.core_logic.constants import AlleleState
from rbceq2.core_logic.utils import Zygosity, sub_alleles_relationships
from rbceq2.db.db import Db, build_antigen_map_for_checks, prepare_db
from rbceq2.main import find_hits, parse_args
from rbceq2.core_logic.data_procesing import (
    explicit_no_data_loci, get_genotypes,
    reference_fallback_blocked_by_no_data,
)

REASON = 'reference_fallback_blocked_by_no_data'


def allele(name, reference=False, variants=frozenset({'1:100_ref'})):
    """Construct synthetic definitions for stage and renderer tests only."""
    return Allele(name, '', '', '', variants, False, reference=reference)


class TestStrictFallbackGuards(unittest.TestCase):
    def setUp(self):
        self.ref = allele('TEST*01', True)
        self.alt = allele('TEST*02', variants=frozenset({'1:200_A_G'}))
        self.references = {'TEST': self.ref}

    def group(self, normal=(), co=None, chrom=2, copies=None):
        """Create a group with explicit pair collections and copy context."""
        return BloodGroup('TEST', {
            AlleleState.NORMAL: list(normal),
            AlleleState.CO: None if co is None else list(co),
        }, 'sample', chrom_copies=chrom, locus_copies=copies)

    def apply(self, bg, strict=True):
        return reference_fallback_blocked_by_no_data(
            {'TEST': bg}, reference_alleles=self.references,
            no_data_loci=frozenset({'1:100'}), strict=strict,
            gene_absent_subtypes={'TEST': 'TEST*ABSENT'},
        )['TEST']

    def test_co_reference_removal_preserves_compound_partner(self):
        compound = allele('TEST*02+TEST*03',
                          variants=frozenset({'1:200_A_G', '1:300_C_T'}))
        pair = Pair(self.ref, compound)
        bg = self.apply(self.group(co=[pair]))
        self.assertIsNone(bg.alleles[AlleleState.CO])
        self.assertEqual(bg.genotypes, ['TEST*02+TEST*03/Undetermined'])
        self.assertIn(pair, bg.filtered_out[REASON])

    def test_callable_pair_remains_without_stale_partial_output(self):
        blocked, allowed = Pair(self.ref, self.alt), Pair(self.alt, self.alt)
        bg = self.group(normal=[blocked, allowed])
        bg.single_slot_genotypes = ['TEST*01/Undetermined']
        self.apply(bg)
        self.assertEqual(bg.alleles[AlleleState.NORMAL], [allowed])
        self.assertEqual(bg.genotypes, ['TEST*02/TEST*02'])
        self.assertEqual(bg.single_slot_genotypes, [])

    def test_unpaired_reference_slot_is_blocked_by_name(self):
        bg = self.group()
        bg.single_slot_genotypes = ['TEST*01/Undetermined']
        self.apply(bg)
        self.assertEqual(bg.genotypes, ['Undetermined/Undetermined'])
        self.assertIn(self.ref, bg.filtered_out[REASON])

    def test_supported_partial_slot_survives(self):
        bg = self.group()
        bg.single_slot_genotypes = ['TEST*02/Undetermined']
        get_genotypes({'TEST': bg})
        self.apply(bg)
        self.assertEqual(bg.genotypes, ['TEST*02/Undetermined'])
        self.assertFalse(bg.filtered_out.get(REASON))

    def test_empty_or_absence_definitions_ignore_unrelated_no_call(self):
        for variants in (frozenset(), frozenset({'1:.'})):
            with self.subTest(variants=variants):
                ref = allele('TEST*01', True, variants)
                bg = self.group(normal=[Pair(ref, ref)])
                reference_fallback_blocked_by_no_data(
                    {'TEST': bg}, reference_alleles={'TEST': ref},
                    no_data_loci=frozenset({'1:100'}), strict=True,
                )
                self.assertTrue(bg.alleles[AlleleState.NORMAL])
                self.assertFalse(bg.filtered_out.get(REASON))

    def test_haploid_refusal_keeps_absent_chromosome_marker(self):
        bg = self.apply(self.group(normal=[Pair(self.ref, self.ref)], chrom=1))
        self.assertEqual(bg.genotypes, ['Undetermined/-'])

    def test_gene_absence_placeholder_is_not_a_reference_claim(self):
        bg = self.group(normal=[Pair(self.ref, self.alt)], copies=1)
        get_genotypes({'TEST': bg}, reference_alleles=self.references,
                      gene_absent_subtypes={'TEST': 'TEST*ABSENT'})
        self.apply(bg)
        self.assertEqual(bg.genotypes, ['TEST*02/TEST*ABSENT'])
        self.assertFalse(bg.filtered_out.get(REASON))

    def test_unresolved_present_gene_keeps_known_absence(self):
        bg = self.apply(self.group(
            normal=[Pair(self.ref, self.ref)], copies=1
        ))
        self.assertEqual(bg.genotypes, ['Undetermined/TEST*ABSENT'])

    def test_no_data_context_does_not_hide_a_called_allele(self):
        vcf = SimpleNamespace(variants={
            '1:100_G_A': {'GT': './.'},
            '1:100_G_T': {'GT': '1/1'},
            '1:200_G_A': {'GT': './.'},
        })
        self.assertEqual(explicit_no_data_loci(vcf), frozenset({'1:200'}))

    def test_late_strict_audit_does_not_change_calls(self):
        """Late strict groups retain missing evidence for audit."""
        db = Db(ref="GRCh38", df=prepare_db())
        row = db.df.loc[db.df["Genotype"] == "CROM*-01"].iloc[0]
        pos, ref, _ = row["GRCh38"].split("_")
        kn = {"KN": [a for a in db.make_alleles() if a.blood_group == "KN"]}
        relations, group = sub_alleles_relationships(kn, "KN")
        results = []
        for debug in (False, True):
            frame = pd.DataFrame([[
                "chr1", pos, ".", ref, "<NON_REF>", ".", ".",
                f"END={pos}", "GT", "./.",
            ]], columns=["CHROM", "POS", "ID", "REF", "ALT", "QUAL",
                         "FILTER", "INFO", "FORMAT", "SAMPLE"])
            args = parse_args(["--reference_genome", "GRCh38",
                               "--strict"]
                              + (["--debug"] if debug else []))
            result = find_hits(
                db, (frame, "audit_sample"), args=args,
                allele_relationships={group: relations},
                excluded=["RHD", "RHCE"],
                ant_mapping=build_antigen_map_for_checks(db.df),
            )
            results.append(result)
        self.assertEqual(results[0][1:4], results[1][1:4])
        bg = results[1][4]["CROM"]
        self.assertEqual(bg.sample, "audit_sample")
        self.assertTrue(any(
            token.startswith(f"1:{pos}_") and state == Zygosity.NO_DATA
            for token, state in bg.unused_pool.items()
        ))
        self.assertTrue(bg.filtered_out[REASON])

    def test_default_mode_does_not_remove_reference(self):
        pair = Pair(self.ref, self.alt)
        bg = self.apply(self.group(normal=[pair]), strict=False)
        self.assertEqual(bg.alleles[AlleleState.NORMAL], [pair])
        self.assertFalse(bg.filtered_out.get(REASON))


if __name__ == '__main__':
    unittest.main()
