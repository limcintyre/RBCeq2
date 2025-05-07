import unittest
from io import StringIO
from unittest.mock import MagicMock, patch

import pandas as pd
from rbceq2.core_logic.alleles import Allele, Line
from rbceq2.db.db import Db, VariantCountMismatchError, compare_antigen_profiles


class TestVariantCountMismatchError(unittest.TestCase):
    """Tests the custom VariantCountMismatchError exception."""

    def test_error_message(self):
        """Ensure the error message is formatted correctly."""
        err = VariantCountMismatchError("A,B", "X,Y,Z")
        self.assertIn("Number of GRCh37 variants must equal", str(err))
        self.assertIn("A,B", str(err))
        self.assertIn("X,Y,Z", str(err))


class TestDb(unittest.TestCase):
    """Tests for the Db dataclass and its methods."""

    def setUp(self):
        """Create minimal CSV data in-memory and patch pd.read_csv to load it."""
        # Note: We include columns used by `Db.prepare_db` and other methods:
        #  - GRCh37, GRCh38, Chrom, Genotype, Phenotype_change,
        #    Genotype_alt, Phenotype_alt_change, Lane, Sub_type,
        #    Weight_of_genotype, Weight_of_phenotype, Reference_genotype, ...
        #  - db code also calls df.query('Lane == True'), so we have a "Lane" col.
        #  - antithetical is read from "Antithetical" col if `Antithetical == 'Yes'`.
        #  - We produce 2 valid rows + 1 reference row + 1 mismatch row for tests.

        self.csv_data = StringIO(
            """GRCh37\tGRCh38\tChrom\tGenotype\tPhenotype_change\tGenotype_alt\tPhenotype_alt_change\tLane\tSub_type\tWeight_of_genotype\tWeight_of_phenotype\tReference_genotype\tAntithetical
1:100_G_A\t1:100_G_A\tchr1\tBG*01.01\tphenoX\tBG*01.01X\tphenoX_alt\tTrue\tSubA\t100\t200\tNo\tNo
1:200_T_C\t1:200_T_C\tchr1\tBG*01.02\tphenoY\tBG*01.02X\tphenoY_alt\tFalse\tSubA\t50\t100\tYes\tNo
1:300_T_C,1:301_A_G\t1:300_T_C\tchr1\tBG*01.03\tphenoZ\t.\t.\tTrue\tSubB\t1\t1\tNo\tYes
2:400_A_G\t2:400_A_G\tchr2\tBG*02.01\t.\t.\t.\tFalse\tSubC\t1\t1\tYes\tNo
"""
        )
        # Explanation of each row:
        # 1) Row0 => Lane=True, Sub_type=SubA, no antithetical, not reference
        # 2) Row1 => Lane=False, same SubA, reference genotype => "Yes"
        # 3) Row2 => MISMATCH in GRCh37 vs GRCh38 => "1:300_T_C,1:301_A_G" vs "1:300_T_C"
        #           also Antithetical="Yes" => tests get_antitheticals
        # 4) Row3 => Lane=False, Sub_type=SubC, also reference => "Yes", no mismatch

        # Build the DataFrame to store separately for manual checks if needed
        self.expected_df = pd.read_csv(self.csv_data, sep="\t")
        self.csv_data.seek(0)  # reset pointer

    @patch("pandas.read_csv")
    def test_prepare_db_and_init(self, mock_read):
        """Test the entire __post_init__ logic, including prepare_db,
        get_antitheticals, get_lane_variants, get_reference_allele,
        and the mismatch check.
        """
        # Mock read_csv to return our expected DataFrame
        mock_read.return_value = self.expected_df.copy()

        # Attempting to init -> should raise VariantCountMismatchError because row2
        # has mismatch
        with self.assertRaises(VariantCountMismatchError):
            Db(ref="Genotype")

        # Now let's remove the mismatch from row2 and try again.
        # We'll unify them so row2 is "1:300_T_C,1:301_A_G" in both GRCh37 & GRCh38:
        good_df = self.expected_df.copy()
        good_df.loc[2, "GRCh38"] = "1:300_T_C,1:301_A_G"
        mock_read.return_value = good_df

        # Now creation should succeed (no mismatch error)
        db_obj = Db(ref="Genotype")

        # Check the final processed DataFrame
        self.assertIsInstance(db_obj.df, pd.DataFrame)
        # We expect 4 rows
        self.assertEqual(db_obj.df.shape[0], 4)

        # Check that the lane_variants is built
        # Row0 => Lane=True => "Chrom=1" => check if the code extracts "BG*01.01"
        # Row2 => Lane=True => "Chrom=1" => check if "BG*01.03" is recognized
        # But note the code only extracts those with "_ref" in them for lane.
        # In our example, none have that substring. So it might produce an empty set.
        # We'll just confirm the dict is created for Chrom=1 though.
        self.assertIn(
            "chr1",
            db_obj.lane_variants,
            "Lane dict should have chromosome '1' due to row0/row2.",
        )
        # The code specifically checks for if variant.endswith('_ref'),
        # so likely the set is empty unless the genotype had `_ref`.
        self.assertEqual(
            len(db_obj.lane_variants["chr1"]),
            0,
            "No '_ref' present, so this set is empty.",
        )

        # Check that reference_alleles is built
        # Row1 => Reference_genotype=Yes => genotype=BG*01.02
        # Row3 => Reference_genotype=Yes => genotype=BG*02.01
        self.assertIn("BG", db_obj.reference_alleles)
        # Actually, key is line.geno.split("*")[0], so e.g. "BG" from "BG*01.02"
        self.assertTrue(len(db_obj.reference_alleles) >= 1)

        # Check that antitheticals is built
        # Row2 => Antithetical=Yes => type = ???. Wait, type is derived from 'Genotype.apply(...)' in prepare_db?
        # Actually the code sets 'df["type"] = df.Genotype.apply(lambda x: x.split("*")[0])'
        # So row2 => genotype="BG*01.03" => type="BG"
        # => antitheticals["BG"] => [list of Genotype from row2 if "Antithetical==Yes"]
        self.assertIn("BG", db_obj.antitheticals)
        self.assertIn("BG*01.03", db_obj.antitheticals["BG"])

        # Confirm no mismatch error was raised
        # => row2 now has matching comma-split values in GRCh37 & GRCh38

    def test_grch37_38_def_var_count_equal_direct(self):
        """Direct test of mismatch logic without full init."""
        # If row 2 has mismatch => raises
        with self.assertRaises(VariantCountMismatchError):
            db_obj = MagicMock()
            db_obj.df = self.expected_df
            # We call the method directly:
            Db.grch37_38_def_var_count_equal(db_obj)

        # Fix mismatch, now no raise
        good_df = self.expected_df.copy()
        good_df.loc[2, "GRCh38"] = "1:300_T_C,1:301_A_G"

        db_obj = MagicMock()
        db_obj.df = good_df
        # Should not raise:
        Db.grch37_38_def_var_count_equal(db_obj)

    @patch("pandas.read_csv")
    def test_line_generator_and_make_alleles(self, mock_read):
        """Test line_generator and make_alleles, verifying the data
        is turned into Allele objects properly.
        """
        # We'll fix row2 mismatch from the start:
        df_local = self.expected_df.copy()
        df_local.loc[2, "GRCh38"] = "1:300_T_C,1:301_A_G"
        mock_read.return_value = df_local

        db_obj = Db(ref="Genotype")

        # line_generator is used inside get_reference_allele and make_alleles,
        # but let's call it ourselves for coverage:
        small_df = db_obj.df.query('Reference_genotype == "Yes"')
        lines = list(db_obj.line_generator(small_df))
        self.assertTrue(len(lines) >= 1)
        self.assertIsInstance(lines[0], Line)

        # Now test make_alleles => yields Allele objects
        all_alleles = list(db_obj.make_alleles())
        # row0 => Genotype=BG*01.01 => Lane=True =>
        # row1 => Genotype=BG*01.02 => reference => included if
        # allele_defining_variants != '.'
        # row2 => ...
        # row3 => ...
        # Some rows might have '.' in "Phenotype_alt_change" => no effect
        # Some might have '.' in "Genotype_alt" => code ignores that for
        #  allele_defining_variants if it's '.'?

        # By default, if line.allele_defining_variants == '.', skip.
        # Let's see how many are skipped.
        # Row3 => has '.' in Genotype_alt => => line.allele_defining_variants => ?

        # Just ensure we got a non-empty list:
        self.assertTrue(
            len(all_alleles) >= 1, "At least row0 or row2 has valid variants."
        )
        for a in all_alleles:
            self.assertIsInstance(a, Allele)


# ── static mapping fixture (small subset, extensible) ─────────────────────
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


class TestAntigenProfileComparison(unittest.TestCase):
    """Unit‑tests exercising every recognised modifier & edge case."""

    # ── helper ------------------------------------------------------------
    def _cmp(self, num: str, alpha: str, system: str, strict: bool = True) -> bool:
        """Wrapper around SUT for brevity inside test bodies."""
        return compare_antigen_profiles(num, alpha, ANTIGEN_MAP, system, strict=strict)

    # ── baseline expression matches --------------------------------------
    def test_basic_positive_negative(self):
        self.assertTrue(self._cmp("RH:-2,3", "C-,E+", "RH"))
        self.assertFalse(self._cmp("RH:-2,3", "C-,E-", "RH"))

    def test_case(self):  # e is real but dif to E
        self.assertTrue(self._cmp("RH:-2,5", "C-,e+", "RH"))
        self.assertFalse(self._cmp("RH:-2,3", "C-,e-", "RH"))

    def test_real(self):  # rd is not real
        self.assertTrue(self._cmp("SC:-4,5", "Rd-,STAR+", "SC"))
        self.assertFalse(self._cmp("SC:-4,5", "rd-,STAR+", "SC"))

    def test_expression(self):  # DSLK needs a +
        self.assertTrue(self._cmp("RHAG:3,5", "DSLK+,Kg+", "RHAG"))
        self.assertFalse(self._cmp("RHAG:3,5", "DSLK,Kg+", "RHAG"))

    def test_order(self):  # "2": "OKGV", "3": "OKVM"},
        self.assertTrue(self._cmp("OK:2,3", "OKGV,OKVM", "OK"))
        self.assertFalse(self._cmp("OK:2,3", "OKVM,OKGV", "OK"))

    # ── individual modifiers ---------------------------------------------
    def test_modifier_weak(self):
        self.assertTrue(self._cmp("RH:4w", "c+weak", "RH"))
        self.assertFalse(  # missing weak in alpha
            self._cmp("RH:4w", "c+", "RH")
        )

    def test_modifier_partial(self):
        self.assertTrue(self._cmp("RH:5p", "e+partial", "RH"))

    def test_modifier_negative(self):
        self.assertTrue(self._cmp("RH:5n", "e+positive_to_neg", "RH"))
        # "positive" alone is *not* accepted → should fail
        self.assertFalse(self._cmp("RH:5n", "e+positive", "RH"))

    def test_modifier_monoclonal(self):
        self.assertTrue(self._cmp("RH:58m", "CELO+monoclonal", "RH"))

    def test_modifier_inferred(self):
        self.assertTrue(self._cmp("RH:31i", "hrB+inferred", "RH"))

    def test_modifier_robust(self):
        self.assertTrue(self._cmp("VEL:1r", "Vel+robust", system="VEL"))

    def test_modifier_strong(self):
        self.assertTrue(self._cmp("VEL:1s", "Vel+strong", system="VEL"))

    # ── multiple modifiers -----------------------------------------------
    def test_modifier_weak_partial(self):
        self.assertTrue(self._cmp("RH:4wp", "c+weak_partial", "RH"))
        # Accept space‑separated order variants
        self.assertTrue(self._cmp("RH:4wp", "c+partial weak", "RH"))

    def test_modifier_weak_partial_negative(self):
        self.assertTrue(self._cmp("RH:5pwn", "e+partial_weak_to_neg", "RH"))

    # ── real‑world composite examples ------------------------------------
    def test_real_world_examples(self):
        self.assertTrue(
            self._cmp(
                "RH:-2,-3,4wp,5wp",  # numeric
                "C-,E-,c+weak_partial,e+weak_partial",
                "RH",
            )
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
                "RH:-2,-3,4,5n,-18,-19,49",
                "C-,E-,c+,e+positive,Hr-,hrS-,STEM+weak",
                "RH",
            )
        )

    # ── strict/extras handling -------------------------------------------
    def test_extra_alpha_antigen_strict(self):
        # α string has an extra antigen (hrB) not in numeric → strict=True should fail
        self.assertFalse(self._cmp("RH:-2,3", "C-,E+,hrB+", "RH", strict=True))

    def test_extra_alpha_antigen_lenient(self):
        self.assertTrue(self._cmp("RH:-2,3", "C-,E+,hrB+", "RH", strict=False))

    # ── failure when required modifier missing ---------------------------
    def test_missing_required_modifier(self):
        # numeric says weak+partial; alpha only weak → should fail
        self.assertFalse(self._cmp("RH:5wp", "e+ weak", "RH"))


if __name__ == "__main__":
    unittest.main()
