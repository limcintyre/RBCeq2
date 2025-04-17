import unittest
from io import StringIO
from unittest.mock import MagicMock, patch

import pandas as pd
from rbceq2.core_logic.alleles import Allele, Line
from rbceq2.db.db import Db, VariantCountMismatchError


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

if __name__ == "__main__":
    unittest.main()
