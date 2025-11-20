import unittest
import pandas as pd
from rbceq2.core_logic.large_variants import (
    _parse_length_unit,
    parse_db_token,
    SvDef,
    SvEvent,
    SvMatcher,
    MatchResult,
    _parse_info,
    _is_large_indel,
    SnifflesVcfSvReader,
    select_best_per_vcf,
)


class TestHelpers(unittest.TestCase):
    def test_parse_length_unit(self):
        self.assertEqual(_parse_length_unit("100"), 100)
        self.assertEqual(_parse_length_unit("1kb"), 1000)
        self.assertEqual(_parse_length_unit("1.5kb"), 1500)
        self.assertEqual(_parse_length_unit("1mb"), 1_000_000)
        self.assertEqual(_parse_length_unit("  500 bp "), 500)

        with self.assertRaises(ValueError):
            _parse_length_unit("100gb")  # Unhandled unit
        with self.assertRaises(ValueError):
            _parse_length_unit("abc")

    def test_parse_info(self):
        info_str = "SVTYPE=DEL;SVLEN=-500;IMPRECISE"
        res = _parse_info(info_str)
        self.assertEqual(res["SVTYPE"], "DEL")
        self.assertEqual(res["SVLEN"], "-500")
        self.assertEqual(res["IMPRECISE"], "True")
        self.assertEqual(_parse_info("."), {})
        self.assertEqual(_parse_info(""), {})

    def test_is_large_indel(self):
        # Threshold = 10
        # Small insertion
        self.assertFalse(_is_large_indel("A", "ACGT", 10))
        # Large insertion
        self.assertTrue(_is_large_indel("A", "A" * 12, 10))
        # Large deletion (REF is long)
        self.assertTrue(_is_large_indel("A" * 12, "A", 10))
        # Symbolic
        self.assertFalse(_is_large_indel("A", "<DEL>", 10))


class TestSvDef(unittest.TestCase):
    def test_parse_db_token_word_form(self):
        # 1:25272547_del_59kb
        res = parse_db_token("1", "25272547_del_59kb", "id1")
        self.assertIsNotNone(res)
        self.assertEqual(res.svtype, "DEL")
        self.assertEqual(res.length, 59000)
        self.assertEqual(res.pos, 25272547)
        self.assertEqual(res.id, "id1")

    def test_parse_db_token_seq_form_del(self):
        # DEL: REF > ALT
        ref = "A" * 20
        alt = "A"
        token = f"100_{ref}_{alt}"
        res = parse_db_token("1", token)
        self.assertEqual(res.svtype, "DEL")
        self.assertEqual(res.length, 19)

    def test_parse_db_token_seq_form_ins(self):
        # INS: ALT > REF
        ref = "A"
        alt = "A" * 20
        token = f"100_{ref}_{alt}"
        res = parse_db_token("1", token)
        self.assertEqual(res.svtype, "INS")
        self.assertEqual(res.length, 19)

    def test_parse_db_token_invalid(self):
        self.assertIsNone(parse_db_token("1", "invalid"))
        self.assertIsNone(parse_db_token("1", "100_bad_unit"))
        self.assertIsNone(parse_db_token("1", "notadigit_del_100"))

    def test_interval_property(self):
        # DEL: span = length
        d = SvDef("1", 100, "DEL", 50, "raw")
        self.assertEqual(d.interval, (100, 150))

        # INS: span = 1 (point)
        i = SvDef("1", 100, "INS", 50, "raw")
        self.assertEqual(i.interval, (100, 101))


class TestSvMatcher(unittest.TestCase):
    def setUp(self):
        self.matcher = SvMatcher()

    def _make_event(self, pos, end, svtype, length):
        return SvEvent(
            chrom="1",
            pos=pos,
            end=end,
            svtype=svtype,
            svlen=length,
            alt="",
            id=".",
            qual=".",
            variant="var1",
            info={},
        )

    def test_compatibility(self):
        db = SvDef("1", 100, "DEL", 100, "raw")
        ev_del = self._make_event(100, 200, "DEL", -100)
        ev_ins = self._make_event(100, 100, "INS", 100)

        self.assertTrue(self.matcher.compatible(db, ev_del))
        self.assertFalse(self.matcher.compatible(db, ev_ins))

    def test_intervals_overlap(self):
        # DB: 100-200 (Len 100)
        db = SvDef("1", 100, "DEL", 100, "raw")

        # EV: 150-250 (Overlap 50bp -> 50% of DB, 50% of EV) -> Should Pass
        ev_overlap = self._make_event(150, 250, "DEL", -100)
        self.assertTrue(self.matcher._intervals_overlap(db, ev_overlap))

        # EV: 300-400 (No overlap)
        ev_no = self._make_event(300, 400, "DEL", -100)
        self.assertFalse(self.matcher._intervals_overlap(db, ev_no))

    def test_scoring_perfect(self):
        db = SvDef("1", 100, "DEL", 100, "raw")
        ev = self._make_event(100, 200, "DEL", -100)  # Perfect match

        score, pos_d, len_d = self.matcher.score(db, ev)
        self.assertEqual(pos_d, 0)
        self.assertEqual(len_d, 0)
        # Bonus applied for overlap?
        # Score = (0) + (0) - 0.5 = -0.5 -> max(0, -0.5) = 0.0
        self.assertEqual(score, 0.0)

    def test_scoring_mismatch_length(self):
        # DB: 100bp.
        # EV: 300bp (Overlap=True).
        # Max Length = 300.
        # Tolerance (Overlap) = 50% of 300 = 150bp.
        # Delta = 200bp.
        # 200 > 150 -> Reject (inf)
        db = SvDef("1", 100, "DEL", 100, "raw")
        ev = self._make_event(100, 400, "DEL", -300) # Length 300
        
        score, _, _ = self.matcher.score(db, ev)
        self.assertEqual(score, float("inf"))


class TestReader(unittest.TestCase):
    def test_sniffles_reader_basic(self):
        # Mock DataFrame
        data = {
            "CHROM": ["1", "1"],
            "POS": [100, 200],
            "ID": [".", "."],
            "REF": ["A", "A"],
            "ALT": ["<DEL>", "<INS>"],
            "QUAL": [".", "."],
            "INFO": ["SVTYPE=DEL;SVLEN=-50;END=150", "SVTYPE=INS;SVLEN=30"],
            "FORMAT": [".", "."],
            "SAMPLE": [".", "."],
            "variant": ["var1", "var2"],
        }
        df = pd.DataFrame(data)
        reader = SnifflesVcfSvReader(df)
        events = list(reader.events())

        self.assertEqual(len(events), 2)
        self.assertEqual(events[0].svtype, "DEL")
        self.assertEqual(events[0].size, 50)
        self.assertEqual(events[1].svtype, "INS")
        self.assertEqual(events[1].size, 30)

    def test_sniffles_reader_inferred_indel(self):
        # No SVTYPE, but huge ref/alt difference
        data = {
            "CHROM": ["1"],
            "POS": [100],
            "ID": ["."],
            "REF": ["A"],
            "ALT": ["A" * 51],  # +50bp INS
            "QUAL": ["."],
            "INFO": ["."],
            "FORMAT": ["."],
            "SAMPLE": ["."],
            "variant": ["var1"],
        }
        df = pd.DataFrame(data)
        reader = SnifflesVcfSvReader(df, min_size=10)
        events = list(reader.events())

        self.assertEqual(len(events), 1)
        self.assertEqual(events[0].svtype, "INS")
        self.assertEqual(events[0].svlen, 50)

    def test_bnd_cache_pairing(self):
        # BNDs usually come in pairs. The reader tries to yield them together.
        data = {
            "CHROM": ["1", "2"],
            "POS": [100, 500],
            "ID": ["bnd1", "bnd2"],
            "REF": ["N", "N"],
            "ALT": ["N]2:500]", "N]1:100]"],
            "QUAL": [".", "."],
            "INFO": ["SVTYPE=BND;MATEID=bnd2", "SVTYPE=BND;MATEID=bnd1"],
            "FORMAT": [".", "."],
            "SAMPLE": [".", "."],
            "variant": ["v1", "v2"],
        }
        df = pd.DataFrame(data)
        reader = SnifflesVcfSvReader(df)
        events = list(reader.events())

        # It should cache bnd1, then when bnd2 appears, yield both.
        self.assertEqual(len(events), 2)
        ids = {e.id for e in events}
        self.assertEqual(ids, {"bnd1", "bnd2"})

    def test_bnd_orphan(self):
        # TODO If a BND has no mate in the file (or logic fails), currently it is lost
        # (based on the code review bug). This test confirms current behavior
        # so if you fix it, this test will need updating (to expect 1 event).
        data = {
            "CHROM": ["1"],
            "POS": [100],
            "ID": ["bnd1"],
            "REF": ["N"],
            "ALT": ["N[2:500["],
            "QUAL": ["."],
            "INFO": ["SVTYPE=BND;MATEID=bnd2"], # bnd2 missing
            "FORMAT": ["."], # FIXED: length 1
            "SAMPLE": ["."], # FIXED: length 1
            "variant": ["v1"]
        }
        df = pd.DataFrame(data)
        reader = SnifflesVcfSvReader(df)
        events = list(reader.events())
        
        # CURRENT BEHAVIOR: Orphan BNDs sitting in cache are lost on exit
        self.assertEqual(len(events), 0) 


class TestSelection(unittest.TestCase):
    def test_select_best_per_vcf(self):
        db = SvDef("1", 100, "DEL", 100, "raw", id="allele1")
        ev = SvEvent("1", 100, 200, "DEL", -100, "", "", "", "v1", {})

        # Match 1: Perfect
        m1 = MatchResult(db, ev, score=0.0, pos_delta=0, len_delta=0, variant="v1")
        # Match 2: Worse score
        m2 = MatchResult(db, ev, score=0.5, pos_delta=10, len_delta=10, variant="v1")

        matches = [m1, m2]
        best = select_best_per_vcf(matches)

        self.assertEqual(len(best), 1)
        self.assertEqual(best[0], m1)

    def test_select_best_tie_breaking(self):
        db1 = SvDef("1", 100, "DEL", 100, "raw", id="A")
        db2 = SvDef("1", 100, "DEL", 100, "raw", id="B")
        ev = SvEvent("1", 100, 200, "DEL", -100, "", "", "", "v1", {})

        # Tie in score
        # m1 has better POS delta
        m1 = MatchResult(db1, ev, score=0.1, pos_delta=2, len_delta=20, variant="v1")
        # m2 has better LEN delta
        m2 = MatchResult(db2, ev, score=0.1, pos_delta=20, len_delta=2, variant="v1")

        # Default delta_weight=0.5 -> equal weight.
        # m1 combined: 2/20 (0.1) + 20/20 (1.0) = 0.55 avg (roughly)
        # m2 combined: 20/20 (1.0) + 2/20 (0.1) = 0.55 avg

        # If identical, it breaks on DB ID.
        matches = [m2, m1]
        best = select_best_per_vcf(matches, tie_tol=0.0)

        self.assertEqual(len(best), 1)
        # Expect A because sorted by DB ID if math is perfectly tied
        self.assertEqual(best[0].db.id, "A")


if __name__ == "__main__":
    unittest.main()
