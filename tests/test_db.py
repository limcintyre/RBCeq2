import unittest
from unittest.mock import patch
import pandas as pd
from rbceq2.core_logic.alleles import Allele, Line

# Updated imports:
from rbceq2.db.db import (
    Db,
    VariantCountMismatchError,
    compare_antigen_profiles,
    DbDataConsistencyChecker,  # New import
    prepare_db,  # New import
    # Keep existing antigen parsing helpers if they are directly tested or needed for test setup
    # For now, assuming they are not directly tested here, but used by compare_antigen_profiles
    # _NUM_ID_RE, _ALPHA_CANON_RE, _canonical_alpha, Antigen, AntigenParser, NumericParser, AlphaParser
)


class TestVariantCountMismatchError(unittest.TestCase):
    def test_error_message(self):
        err = VariantCountMismatchError("A,B", "X,Y,Z")
        self.assertIn("Number of GRCh37 variants must equal", str(err))
        self.assertIn("A,B", str(err))
        self.assertIn("X,Y,Z", str(err))


class TestDb(unittest.TestCase):
    def setUp(self):
        self.headers = "GRCh37\tGRCh38\tChrom\tGenotype\tPhenotype_change\tGenotype_alt\tPhenotype_alt_change\tLane\tSub_type\tWeight_of_genotype\tReference_genotype\tAntithetical\tPhenotype\tPhenotype_alt\n"

        self.good_csv_content = (
            self.headers
            + """1:100_G_A\t1:100_G_A\tchr1\tBG*01.01\tSYS:1w\tBG*01.01X\tSYS:X+weak\tTrue\tSubA\t100\tNo\tNo\tSYS:1w\tSYS:X+weak
1:200_T_C\t1:200_T_C\tchr1\tREF*01.02\tSYS:2\tREF*01.02X\tSYS:Y+\tFalse\tSubA\t50\tYes\tNo\tSYS:2\tSYS:Y+
1:300_T_C,1:301_A_G\t1:300_T_C,1:301_A_G\tchr1\tANTI*01.03\tSYS:3\tANTI*01.03X\tSYS:Z+\tTrue\tSubB\t1\tNo\tYes\tSYS:3\tSYS:Z+
2:400_A_G\t2:400_A_G\tchr2\tREF*02.01\t.\t.\t.\tFalse\tSubC\t1\tYes\tNo\t.\t.
"""
        )
        self.mismatch_grch_csv_content = (
            self.headers
            + """1:300_T_C,1:301_A_G\t1:300_T_C\tchr1\tBG*01.03\t.\t.\t.\tTrue\tSubB\t1\tNo\tYes\t.\t.
"""
        )
        # KEL:1 (numeric, positive) vs KEL:K- (alpha, negative) for Phenotype_change/Phenotype_alt_change
        # The Phenotype/Phenotype_alt are KEL:1 vs KEL:K+ (consistent) to not trigger mismatch there
        self.antigen_mismatch_csv_content = (
            self.headers
            + """1:100_G_A\t1:100_G_A\tchr1\tBG*01.01\tKEL:1\tBG*01.01X\tKEL:K-\tTrue\tSubA\t100\tNo\tNo\tKEL:1\tKEL:K+
"""
        )
        self.mock_load_db_path = "rbceq2.db.db.load_db"

        global ANTIGEN_MAP  # Use the global ANTIGEN_MAP defined later in the file
        if "SYS" not in ANTIGEN_MAP:
            ANTIGEN_MAP["SYS"] = {"1": "X", "2": "Y", "3": "Z"}
        if "KEL" not in ANTIGEN_MAP:
            ANTIGEN_MAP["KEL"] = {}
        ANTIGEN_MAP["KEL"]["1"] = "K"


    @patch("rbceq2.db.db.load_db")
    def test_db_creation_and_consistency_checks(self, mock_load_db_func):
        # Scenario 1: GRCh37/38 mismatch
        mock_load_db_func.return_value = self.mismatch_grch_csv_content
        df_mismatch = prepare_db()
        with self.assertRaises(VariantCountMismatchError):
            DbDataConsistencyChecker.run_all_checks(
                df=df_mismatch, ref_genome_name="GRCh37"
            )

        # Scenario 2: Antigen profile mismatch for Phenotype_change/Phenotype_alt_change
        mock_load_db_func.return_value = self.antigen_mismatch_csv_content
        df_antigen_mismatch = prepare_db()
        with self.assertRaises(
            AssertionError
        ):  # Expecting mismatch leading to AssertionError
            DbDataConsistencyChecker.run_all_checks(
                df=df_antigen_mismatch, ref_genome_name="GRCh37"
            )

        # Scenario 3: Good data
        mock_load_db_func.return_value = self.good_csv_content
        df_good = prepare_db()
        DbDataConsistencyChecker.run_all_checks(
            df=df_good, ref_genome_name="GRCh37"
        )  # Should pass
        db_obj = Db(ref="Genotype", df=df_good)

        self.assertIsInstance(db_obj.df, pd.DataFrame)
        self.assertEqual(db_obj.df.shape[0], 4)
        self.assertIn("chr1", db_obj.lane_variants)
        self.assertEqual(len(db_obj.lane_variants["chr1"]), 0)
        self.assertIn("REF", db_obj.reference_alleles)
        self.assertEqual(len(db_obj.reference_alleles), 1)
        self.assertEqual(db_obj.reference_alleles["REF"].genotype, "REF*02.01")
        self.assertIn("ANTI", db_obj.antitheticals)
        self.assertIn("ANTI*01.03", db_obj.antitheticals["ANTI"])

    def test_check_grch37_38_variant_counts_direct(self):
        df_mismatch = pd.DataFrame(
            {
                "GRCh37": ["1:300_T_C,1:301_A_G", "1:1_A_T"],
                "GRCh38": ["1:300_T_C", "1:1_A_T,1:2_C_G"],
            }
        )
        with self.assertRaises(VariantCountMismatchError) as cm:
            DbDataConsistencyChecker.check_grch37_38_variant_counts(df_mismatch)
        self.assertIn("1:300_T_C,1:301_A_G", str(cm.exception))
        self.assertIn("1:300_T_C", str(cm.exception))

        df_good = pd.DataFrame(
            {
                "GRCh37": ["1:300_T_C,1:301_A_G", "1:1_A_T"],
                "GRCh38": ["1:300_T_C,1:301_A_G", "1:1_A_T"],
            }
        )
        DbDataConsistencyChecker.check_grch37_38_variant_counts(df_good)

        df_empty_vs_nonempty = pd.DataFrame({"GRCh37": ["."], "GRCh38": ["1:1_A_T"]})
        with self.assertRaises(VariantCountMismatchError):
            DbDataConsistencyChecker.check_grch37_38_variant_counts(
                df_empty_vs_nonempty
            )

        df_both_empty = pd.DataFrame({"GRCh37": ["."], "GRCh38": ["."]})
        DbDataConsistencyChecker.check_grch37_38_variant_counts(df_both_empty)

    @patch("rbceq2.db.db.load_db")
    def test_line_generator_and_make_alleles(self, mock_load_db_func):
        mock_load_db_func.return_value = self.good_csv_content
        df = prepare_db()
        DbDataConsistencyChecker.run_all_checks(df=df, ref_genome_name="GRCh37")
        db_obj = Db(ref="Genotype", df=df)

        ref_df_subset = db_obj.df.query('Reference_genotype == "Yes"')
        lines = list(db_obj.line_generator(ref_df_subset))
        self.assertEqual(len(lines), 2)
        self.assertIsInstance(lines[0], Line)
        self.assertEqual(lines[0].geno, "REF*01.02")

        all_alleles = list(db_obj.make_alleles())
        self.assertEqual(len(all_alleles), 4)
        for a in all_alleles:
            self.assertIsInstance(a, Allele)

        allele_bg0101 = next(a for a in all_alleles if a.genotype == "BG*01.01")
        # The defining_variants construction adds chrom: before the variant string
        # and assumes GRCh37 column from self.ref being Genotype, and Line.allele_defining_variants
        # getting self.df[self.ref] as input, which is 'Genotype'.
        # This seems like a slight mix-up. `make_alleles` actually uses `line.allele_defining_variants`
        # which is populated from `df[self.ref]`. If `self.ref` is "Genotype", then `line.allele_defining_variants`
        # will be "BG*01.01", not "1:100_G_A".
        # The Line dataclass has `allele_defining_variants: str` which gets `df[self.ref]`.
        # This should be df['GRCh37'] or df['GRCh38'] depending on the desired reference for defining variants.
        # Let's assume `self.ref` is intended to be the column containing the defining variants like "1:100_G_A"
        # For `make_alleles` to work as intended (creating defining_variants like "chr1:1:100_G_A"),
        # `db_obj` should be initialized with `ref="GRCh37"` (or "GRCh38").
        # Let's re-initialize db_obj for this specific test part.

        db_obj_for_alleles = Db(
            ref="GRCh37", df=df
        )  # Use GRCh37 as the source for defining variants
        all_alleles_corrected = list(db_obj_for_alleles.make_alleles())
        allele_bg0101_corrected = next(
            a for a in all_alleles_corrected if a.genotype == "BG*01.01"
        )
        # line.allele_defining_variants will be "1:100_G_A" from df["GRCh37"]
        # and line.chrom will be "chr1"
        self.assertEqual(
            allele_bg0101_corrected.defining_variants, frozenset({"1:1:100_G_A"})
        )


# ANTIGEN_MAP: dict[str, dict[str, str]] = {
#     "ABCC1": {"1": "WLF"}, "ATP11C": {"1": "LIL"},
#     "AUG": {"1": "AUG1", "2": "At(a", "3": "ATML", "4": "ATAM"},
#     "CD59": {"1": "CD59.1"},
#     "CO": {"1": "Co(a", "2": "Co(b", "3": "Co3", "4": "Co4"},
#     "CROM": {"1": "Cr(a", "10": "UMC", "11": "GUTI", "12": "SERF", "13": "ZENA", "14": "CROV", "15": "CRAM", "16": "CROZ", "17": "CRUE", "18": "CRAG", "19": "CROK", "2": "Tc(a", "20": "CORS", "3": "Tc(b", "4": "Tc(c", "5": "Dr(a", "6": "Es(a", "7": "IFC", "8": "WES(a", "9": "WES(b"},
#     "CTL2": {"1": "VER", "2": "RIF", "3": "Cs(a", "4": "Cs(b", "5": "BROS"},
#     "DI": {"1": "Di(a", "10": "Bp(a", "11": "Mo(a", "12": "Hg(a", "13": "Vg(a", "14": "Sw(a", "15": "BOW", "16": "NFLD", "17": "Jn(a", "18": "KREP", "19": "Tr(a", "2": "Di(b", "20": "Fr(a", "21": "SW1", "22": "DISK", "23": "DIST", "3": "Wr(a", "4": "Wr(b", "5": "Wd(a", "6": "Rb(a", "7": "WARR", "8": "ELO", "9": "Wu"},
#     "DO": {"1": "Do(a", "10": "DODE", "2": "Do(b", "3": "Gy(a", "4": "Hy", "5": "Jo(a", "6": "DOYA", "7": "DOMR", "8": "DOLG", "9": "DOLC"},
#     "EMM": {"1": "Emm"},
#     "ER": {"1": "Er(a", "2": "Er(b", "3": "Er3", "4": "ERSA", "5": "ERAMA"},
#     "FORS": {"1": "FORS"},
#     "FY": {"1": "Fy(a", "2": "Fy(b", "3": "Fy3", "5": "Fy5", "6": "Fy6"},
#     "GE": {"10": "GEPL", "11": "GEAT", "12": "GETI", "13": "GECT", "14": "GEAR", "2": "Ge2", "3": "Ge3", "4": "Ge4", "5": "Wb", "6": "Ls(a", "7": "An(a", "8": "Dh(a", "9": "GEIS"},
#     "GIL": {"1": "GIL"}, "GLOB": {"1": "P"}, "I": {"1": "I"},
#     "IN": {"1": "In(a", "2": "In(b", "3": "INFI", "4": "INJA", "5": "INRA", "6": "INSL"},
#     "JK": {"1": "Jk(a", "2": "Jk(b", "3": "Jk3"},
#     "JMH": {"1": "JMH", "2": "JMHK", "3": "JMHL", "4": "JMHG", "5": "JMHM", "6": "JMHQ", "7": "JMHN", "8": "JMHA"},
#     "KANNO": {"1": "KANNO1"},
#     "KEL": {"1": "K", "10": "Ul(a", "11": "K11", "12": "K12", "13": "K13", "14": "K14", "16": "K16", "17": "K17", "18": "K18", "19": "K19", "2": "k", "20": "Km", "21": "Kp(c", "22": "K22", "23": "K23", "24": "K24", "25": "VLAN", "26": "TOU", "27": "RAZ", "28": "VONG", "29": "KALT", "3": "Kp(a", "30": "KTIM", "31": "KYO", "32": "KUCI", "33": "KANT", "34": "KASH", "35": "KELP", "36": "KETI", "37": "KHUL", "38": "KYOR", "39": "KEAL", "4": "Kp(b", "40": "KHIZ", "41": "KHOZ", "5": "Ku", "6": "Js(a", "7": "Js(b"},
#     "KN": {"1": "Kn(a", "10": "KDAS", "11": "DACY", "12": "YCAD", "13": "KNMB", "2": "Kn(b", "3": "McC(a", "4": "Sl1", "5": "Yk(a", "6": "McC(b", "7": "Vil", "8": "Sl3", "9": "KCAM"},
#     "LU": {"1": "Lu(a", "11": "Lu11", "12": "Lu12", "13": "Lu13", "14": "Lu14", "16": "Lu16", "17": "Lu17", "18": "Au(a", "19": "Au(b", "2": "Lu(b", "20": "Lu20", "21": "Lu21", "22": "LURC", "23": "LUIT", "24": "LUGA", "25": "LUAC", "26": "LUBI", "27": "LUYA", "28": "LUNU", "29": "LURA", "3": "Lu3", "30": "LUOM", "4": "Lu4", "5": "Lu5", "6": "Lu6", "7": "Lu7", "8": "Lu8", "9": "Lu9"},
#     "LW": {"5": "LWa", "6": "LWab", "7": "LWb", "8": "LWEM"},
#     "MAM": {"1": "MAM"},
#     "MNS": {"1": "M", "10": "Mur", "11": "Mg", "12": "Vr", "13": "Me", "14": "Mt(a", "15": "St(a", "16": "Ri(a", "17": "Cl(a", "18": "Ny(a", "19": "Hut", "2": "N", "20": "Hil", "21": "Mv", "22": "Far", "23": "sD", "24": "Mit", "25": "Dantu", "26": "Hop", "27": "Nob", "28": "En(a", "29": "ENKT", "3": "S", "30": "`N'", "31": "Or", "32": "DANE", "33": "TSEN", "34": "MINY", "35": "MUT", "36": "SAT", "37": "ERIK", "38": "Os(a", "39": "ENEP", "4": "s", "40": "ENEH", "41": "HAG", "42": "ENAV", "43": "MARS", "44": "ENDA", "45": "ENEV", "46": "MNTD", "47": "SARA", "48": "KIPP", "49": "JENU", "5": "U", "50": "SUMI", "6": "He", "7": "Mia", "8": "Mc", "9": "Vw"},
#     "OK": {"1": "Ok(a", "2": "OKGV", "3": "OKVM"},
#     "PEL": {"1": "PEL"}, "RAPH": {"1": "MER2"},
#     "RH": {"10": "V", "11": "Ew", "12": "G", "17": "Hro", "18": "Hr", "19": "hrS", "2": "C", "20": "VS", "21": "CG", "22": "CE", "23": "Dw", "26": "c_like", "27": "cE", "28": "hrH", "29": "Rh29", "3": "E", "30": "Goa", "31": "hrB", "32": "Rh32", "33": "Rh33", "34": "HrB", "35": "Rh35", "36": "Bea", "37": "Evans", "39": "Rh39", "4": "c", "40": "Tar", "41": "Rh41", "42": "Rh42", "43": "Crawford", "44": "Nou", "45": "Riv", "46": "Sec", "47": "Dav", "48": "JAL", "49": "STEM", "5": "e", "50": "FPTT", "51": "MAR", "52": "BARC", "53": "JAHK", "54": "DAK", "55": "LOCR", "56": "CENR", "57": "CEST", "58": "CELO", "59": "CEAG", "6": "f", "60": "PARG", "61": "CEVF", "62": "CEWA", "63": "CETW", "7": "Ce", "8": "Cw", "9": "Cx"},
#     "RHAG": {"1": "Duclos", "2": "Ol(a", "3": "DSLK", "5": "Kg", "6": "SHER", "7": "THIN"},
#     "RHD": {"1": "D"},
#     "SC": {"1": "Sc1", "2": "Sc2", "3": "Sc3", "4": "Rd", "5": "STAR", "6": "SCER", "7": "SCAN", "8": "SCAR", "9": "SCAC"},
#     "SID": {"1": "Sd(a"}, "VEL": {"1": "Vel"}, "XK": {"1": "kx"},
#     "YT": {"1": "Yt(a", "2": "Yt(b", "3": "YTEG", "4": "YTLI", "5": "YTOT"},
# }


class TestAntigenProfileComparison(unittest.TestCase):
    def _cmp(self, num: str, alpha: str, system: str, strict: bool = True) -> bool:
        return compare_antigen_profiles(num, alpha, ANTIGEN_MAP, system, strict=strict)

    def test_basic_positive_negative(self):
        self.assertTrue(self._cmp("RH:-2,3", "C-,E+", "RH"))
        self.assertFalse(self._cmp("RH:-2,3", "C-,E-", "RH"))

    def test_case_sensitivity_in_alpha_name(self):
        self.assertTrue(self._cmp("RH:-2,5", "C-,e+", "RH"))
        self.assertFalse(self._cmp("RH:-2,3", "C-,e+", "RH"))

    def test_alpha_name_validity(self):
        self.assertTrue(self._cmp("SC:-4,5", "Rd-,STAR+", "SC"))
        self.assertFalse(
            self._cmp("SC:-4,5", "rd-,STAR+", "SC")
        )  # 'rd' is not the canonical 'Rd'

    def test_correctly_signed_alpha_expression(self):
        self.assertTrue(self._cmp("RHAG:3,5", "DSLK+,Kg+", "RHAG"))
        self.assertTrue(self._cmp("OK:2,3", "OKGV+,OKVM+", "OK"))

    def test_order_sensitivity(self):
        self.assertTrue(self._cmp("OK:2,3", "OKGV+,OKVM+", "OK"))
        self.assertTrue(
            self._cmp("OK:2,3", "OKVM+,OKGV+", "OK")
        )  # Set comparison, order insensitive

    def test_modifier_weak(self):
        self.assertTrue(self._cmp("RH:4w", "c+weak", "RH"))
        self.assertFalse(self._cmp("RH:4w", "c+", "RH"))

    def test_modifier_partial(self):
        self.assertTrue(self._cmp("RH:5p", "e+partial", "RH"))

    def test_modifier_negative_transition(self):
        self.assertTrue(self._cmp("RH:5n", "e+positive_to_neg", "RH"))
        self.assertTrue(self._cmp("RH:5n", "e+negative", "RH"))
        self.assertFalse(self._cmp("RH:5n", "e+positive", "RH"))

    def test_modifier_monoclonal(self):
        self.assertTrue(self._cmp("RH:58m", "CELO+monoclonal", "RH"))

    def test_modifier_inferred(self):
        self.assertTrue(self._cmp("RH:31i", "hrB+inferred", "RH"))

    def test_modifier_robust(self):
        self.assertTrue(self._cmp("VEL:1r", "Vel+robust", system="VEL"))

    def test_modifier_strong(self):
        self.assertTrue(self._cmp("VEL:1s", "Vel+strong", system="VEL"))

    def test_modifier_very_weak(self):
        ANTIGEN_MAP["KEL"]["2"] = "k"
        self.assertTrue(self._cmp("KEL:2v", "k+very_weak", system="KEL"))
        self.assertTrue(self._cmp("KEL:2v", "k+v", system="KEL"))

    def test_modifier_weak_partial(self):
        self.assertTrue(self._cmp("RH:4wp", "c+weak_partial", "RH"))
        self.assertTrue(self._cmp("RH:4wp", "c+partial weak", "RH"))
        self.assertTrue(self._cmp("RH:4wp", "c+wp", "RH"))

    def test_modifier_weak_partial_negative(self):
        self.assertTrue(self._cmp("RH:5pwn", "e+partial_weak_to_neg", "RH"))
        self.assertTrue(self._cmp("RH:5pwn", "e+pwn", "RH"))

    def test_real_world_examples(self):
        self.assertTrue(
            self._cmp("RH:-2,-3,4wp,5wp", "C-,E-,c+weak_partial,e+weak_partial", "RH")
        )
        self.assertTrue(
            self._cmp(
                "RH:-2,-3,4,5n,-18,-19,49w",
                "C-,E-,c+,e+positive_to_neg,Hr-,hrS-,STEM+weak",
                "RH",
            )
        )
        self.assertFalse(
            self._cmp(
                "RH:-2,-3,4,5n,-18,-19,49",
                "C-,E-,c+,e+positive_to_neg,Hr-,hrS-,STEM+weak",
                "RH",
            )
        )
        self.assertFalse(
            self._cmp(
                "RH:-2,-3,4,5n,-18,-19,49w",
                "C-,E-,c+,e+positive,Hr-,hrS-,STEM+weak",
                "RH",
            )
        )

    def test_extra_alpha_antigen_strict(self):
        self.assertFalse(self._cmp("RH:-2,3", "C-,E+,hrB+", "RH", strict=True))

    def test_extra_alpha_antigen_lenient(self):
        self.assertTrue(self._cmp("RH:-2,3", "C-,E+,hrB+", "RH", strict=False))

    def test_missing_required_modifier(self):
        self.assertFalse(self._cmp("RH:5wp", "e+weak", "RH"))


ANTIGEN_MAP: dict[str, dict[str, str]] = {
    "ABCC1": {"1": "WLF"},
    "ATP11C": {"1": "LIL"},
    "AUG": {"1": "AUG1", "2": "At(a", "3": "ATML", "4": "ATAM"},
    "CD59": {"1": "CD59.1"},
    "CO": {"1": "Co(a", "2": "Co(b", "3": "Co3", "4": "Co4"},
    "CROM": {
        "1": "Cr(a",
        "10": "UMC",
        "11": "GUTI",
        "12": "SERF",
        "13": "ZENA",
        "14": "CROV",
        "15": "CRAM",
        "16": "CROZ",
        "17": "CRUE",
        "18": "CRAG",
        "19": "CROK",
        "2": "Tc(a",
        "20": "CORS",
        "3": "Tc(b",
        "4": "Tc(c",
        "5": "Dr(a",
        "6": "Es(a",
        "7": "IFC",
        "8": "WES(a",
        "9": "WES(b",
    },
    "CTL2": {"1": "VER", "2": "RIF", "3": "Cs(a", "4": "Cs(b", "5": "BROS"},
    "DI": {
        "1": "Di(a",
        "10": "Bp(a",
        "11": "Mo(a",
        "12": "Hg(a",
        "13": "Vg(a",
        "14": "Sw(a",
        "15": "BOW",
        "16": "NFLD",
        "17": "Jn(a",
        "18": "KREP",
        "19": "Tr(a",
        "2": "Di(b",
        "20": "Fr(a",
        "21": "SW1",
        "22": "DISK",
        "23": "DIST",
        "3": "Wr(a",
        "4": "Wr(b",
        "5": "Wd(a",
        "6": "Rb(a",
        "7": "WARR",
        "8": "ELO",
        "9": "Wu",
    },
    "DO": {
        "1": "Do(a",
        "10": "DODE",
        "2": "Do(b",
        "3": "Gy(a",
        "4": "Hy",
        "5": "Jo(a",
        "6": "DOYA",
        "7": "DOMR",
        "8": "DOLG",
        "9": "DOLC",
    },
    "EMM": {"1": "Emm"},
    "ER": {"1": "Er(a", "2": "Er(b", "3": "Er3", "4": "ERSA", "5": "ERAMA"},
    "FORS": {"1": "FORS"},
    "FY": {"1": "Fy(a", "2": "Fy(b", "3": "Fy3", "5": "Fy5", "6": "Fy6"},
    "GE": {
        "10": "GEPL",
        "11": "GEAT",
        "12": "GETI",
        "13": "GECT",
        "14": "GEAR",
        "2": "Ge2",
        "3": "Ge3",
        "4": "Ge4",
        "5": "Wb",
        "6": "Ls(a",
        "7": "An(a",
        "8": "Dh(a",
        "9": "GEIS",
    },
    "GIL": {"1": "GIL"},
    "GLOB": {"1": "P"},
    "I": {"1": "I"},
    "IN": {
        "1": "In(a",
        "2": "In(b",
        "3": "INFI",
        "4": "INJA",
        "5": "INRA",
        "6": "INSL",
    },
    "JK": {"1": "Jk(a", "2": "Jk(b", "3": "Jk3"},
    "JMH": {
        "1": "JMH",
        "2": "JMHK",
        "3": "JMHL",
        "4": "JMHG",
        "5": "JMHM",
        "6": "JMHQ",
        "7": "JMHN",
        "8": "JMHA",
    },
    "KANNO": {"1": "KANNO1"},
    "KEL": {
        "1": "K",
        "10": "Ul(a",
        "11": "K11",
        "12": "K12",
        "13": "K13",
        "14": "K14",
        "16": "K16",
        "17": "K17",
        "18": "K18",
        "19": "K19",
        "2": "k",
        "20": "Km",
        "21": "Kp(c",
        "22": "K22",
        "23": "K23",
        "24": "K24",
        "25": "VLAN",
        "26": "TOU",
        "27": "RAZ",
        "28": "VONG",
        "29": "KALT",
        "3": "Kp(a",
        "30": "KTIM",
        "31": "KYO",
        "32": "KUCI",
        "33": "KANT",
        "34": "KASH",
        "35": "KELP",
        "36": "KETI",
        "37": "KHUL",
        "38": "KYOR",
        "39": "KEAL",
        "4": "Kp(b",
        "40": "KHIZ",
        "41": "KHOZ",
        "5": "Ku",
        "6": "Js(a",
        "7": "Js(b",
    },
    "KN": {
        "1": "Kn(a",
        "10": "KDAS",
        "11": "DACY",
        "12": "YCAD",
        "13": "KNMB",
        "2": "Kn(b",
        "3": "McC(a",
        "4": "Sl1",
        "5": "Yk(a",
        "6": "McC(b",
        "7": "Vil",
        "8": "Sl3",
        "9": "KCAM",
    },
    "LU": {
        "1": "Lu(a",
        "11": "Lu11",
        "12": "Lu12",
        "13": "Lu13",
        "14": "Lu14",
        "16": "Lu16",
        "17": "Lu17",
        "18": "Au(a",
        "19": "Au(b",
        "2": "Lu(b",
        "20": "Lu20",
        "21": "Lu21",
        "22": "LURC",
        "23": "LUIT",
        "24": "LUGA",
        "25": "LUAC",
        "26": "LUBI",
        "27": "LUYA",
        "28": "LUNU",
        "29": "LURA",
        "3": "Lu3",
        "30": "LUOM",
        "4": "Lu4",
        "5": "Lu5",
        "6": "Lu6",
        "7": "Lu7",
        "8": "Lu8",
        "9": "Lu9",
    },
    "LW": {"5": "LWa", "6": "LWab", "7": "LWb", "8": "LWEM"},
    "MAM": {"1": "MAM"},
    "MNS": {
        "1": "M",
        "10": "Mur",
        "11": "Mg",
        "12": "Vr",
        "13": "Me",
        "14": "Mt(a",
        "15": "St(a",
        "16": "Ri(a",
        "17": "Cl(a",
        "18": "Ny(a",
        "19": "Hut",
        "2": "N",
        "20": "Hil",
        "21": "Mv",
        "22": "Far",
        "23": "sD",
        "24": "Mit",
        "25": "Dantu",
        "26": "Hop",
        "27": "Nob",
        "28": "En(a",
        "29": "ENKT",
        "3": "S",
        "30": "`N'",
        "31": "Or",
        "32": "DANE",
        "33": "TSEN",
        "34": "MINY",
        "35": "MUT",
        "36": "SAT",
        "37": "ERIK",
        "38": "Os(a",
        "39": "ENEP",
        "4": "s",
        "40": "ENEH",
        "41": "HAG",
        "42": "ENAV",
        "43": "MARS",
        "44": "ENDA",
        "45": "ENEV",
        "46": "MNTD",
        "47": "SARA",
        "48": "KIPP",
        "49": "JENU",
        "5": "U",
        "50": "SUMI",
        "6": "He",
        "7": "Mia",
        "8": "Mc",
        "9": "Vw",
    },
    "OK": {"1": "Ok(a", "2": "OKGV", "3": "OKVM"},
    "PEL": {"1": "PEL"},
    "RAPH": {"1": "MER2"},
    "RH": {
        "10": "V",
        "11": "Ew",
        "12": "G",
        "17": "Hro",
        "18": "Hr",
        "19": "hrS",
        "2": "C",
        "20": "VS",
        "21": "CG",
        "22": "CE",
        "23": "Dw",
        "26": "c_like",
        "27": "cE",
        "28": "hrH",
        "29": "Rh29",
        "3": "E",
        "30": "Goa",
        "31": "hrB",
        "32": "Rh32",
        "33": "Rh33",
        "34": "HrB",
        "35": "Rh35",
        "36": "Bea",
        "37": "Evans",
        "39": "Rh39",
        "4": "c",
        "40": "Tar",
        "41": "Rh41",
        "42": "Rh42",
        "43": "Crawford",
        "44": "Nou",
        "45": "Riv",
        "46": "Sec",
        "47": "Dav",
        "48": "JAL",
        "49": "STEM",
        "5": "e",
        "50": "FPTT",
        "51": "MAR",
        "52": "BARC",
        "53": "JAHK",
        "54": "DAK",
        "55": "LOCR",
        "56": "CENR",
        "57": "CEST",
        "58": "CELO",
        "59": "CEAG",
        "6": "f",
        "60": "PARG",
        "61": "CEVF",
        "62": "CEWA",
        "63": "CETW",
        "7": "Ce",
        "8": "Cw",
        "9": "Cx",
    },
    "RHAG": {
        "1": "Duclos",
        "2": "Ol(a",
        "3": "DSLK",
        "5": "Kg",
        "6": "SHER",
        "7": "THIN",
    },
    "RHD": {"1": "D"},
    "SC": {
        "1": "Sc1",
        "2": "Sc2",
        "3": "Sc3",
        "4": "Rd",
        "5": "STAR",
        "6": "SCER",
        "7": "SCAN",
        "8": "SCAR",
        "9": "SCAC",
    },
    "SID": {"1": "Sd(a"},
    "VEL": {"1": "Vel"},
    "XK": {"1": "kx"},
    "YT": {"1": "Yt(a", "2": "Yt(b", "3": "YTEG", "4": "YTLI", "5": "YTOT"},
}

if __name__ == "__main__":
    unittest.main()
