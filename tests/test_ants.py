import unittest

from rbceq2.phenotype.antigens import (
    AlphaNumericAntigen,
    AlphaNumericAntigenABO,
    AlphaNumericAntigenMNS,
    AlphaNumericAntigenVel,
    AlphaNumericAntigenXG,
    NumericAntigen,
)

# For the purpose of this test, I'll define the classes inline
# (Include your provided classes here or import them from your module)
# used chatGPT o1 preview in october 2024 (10 mins work, which seems amazing at this
# point in time!)



class TestAntigenClasses(unittest.TestCase):
    def test_numeric_antigen_initialization(self):
        antigen = NumericAntigen(
            given_name="3",
            expressed=True,
            homozygous=True,
            antithetical_relationships={"3": ["4"]},
        )
        self.assertEqual(antigen.base_name, "3")
        self.assertFalse(antigen.weak)
        self.assertEqual(antigen.weight, 1)
        self.assertTrue(antigen.expressed)
        self.assertEqual(antigen.name, "3")
        self.assertEqual(antigen.antithetical_antigen[0].given_name, "4")

    def test_numeric_antigen_weight_with_w(self):
        """Test that _set_weight returns 2 when given_name contains 'w'."""
        antigen = NumericAntigen(
            given_name="Aw",
            expressed=True,
            homozygous=True,
            antithetical_relationships={},
        )
        self.assertEqual(antigen.weight, 2)

    def test_numeric_antigen_weight_with_dash(self):
        """Test that _set_weight returns 3 when given_name contains '-'."""
        antigen = NumericAntigen(
            given_name="A-",
            expressed=True,
            homozygous=True,
            antithetical_relationships={},
        )
        self.assertEqual(antigen.weight, 3)


    def test_alphanumeric_antigen_initialization(self):
        antigen = AlphaNumericAntigen(
            given_name="K+",
            expressed=True,
            homozygous=False,
            antithetical_relationships={"K": ["k"]},
        )
        self.assertEqual(antigen.base_name, "K")
        self.assertFalse(antigen.weak)
        self.assertEqual(antigen.weight, 1)
        self.assertEqual(antigen.name, "K+")
        self.assertEqual(antigen.antithetical_antigen[0].given_name, "k")

    def test_alphanumeric_antigen_weak(self):
        antigen = AlphaNumericAntigen(
            given_name="K+w",
            expressed=True,
            homozygous=False,
            antithetical_relationships={},
        )
        self.assertTrue(antigen.weak)
        self.assertEqual(antigen.weight, 2)
        self.assertEqual(antigen.name, "K+w")

    def test_alphanumeric_antigen_negative(self):
        antigen = AlphaNumericAntigen(
            given_name="K-",
            expressed=False,
            homozygous=False,
            antithetical_relationships={},
        )
        self.assertFalse(antigen.weak)
        self.assertEqual(antigen.weight, 3)
        self.assertEqual(antigen.name, "K-")

    def test_alphanumeric_antigen_ABO(self):
        antigen = AlphaNumericAntigenABO(
            given_name="A1",
            expressed=True,
            homozygous=False,
            antithetical_relationships={},
        )
        self.assertEqual(antigen.base_name, "A")
        self.assertFalse(antigen.weak)
        self.assertEqual(antigen.name, "A1")

    def test_alphanumeric_antigen_ABO_invalid(self):
        with self.assertRaises(ValueError):
            AlphaNumericAntigenABO(
                given_name="Invalid",
                expressed=True,
                homozygous=False,
                antithetical_relationships={},
            )

    def test_alphanumeric_antigen_XG(self):
        antigen = AlphaNumericAntigenXG(
            given_name="Xg(a+)",
            expressed=True,
            homozygous=False,
            antithetical_relationships={},
        )
        expected_name = "XgaCD99"
        self.assertEqual(antigen.name, expected_name)

    def test_alphanumeric_antigen_MNS(self):
        antigen = AlphaNumericAntigenMNS(
            given_name="M+",
            expressed=True,
            homozygous=False,
            antithetical_relationships={"M": ["N"]},
        )
        self.assertEqual(antigen.base_name, "M")
        self.assertEqual(antigen.name, "M+")
        self.assertEqual(antigen.antithetical_antigen[0].given_name, "N")

    def test_alphanumeric_antigen_Vel(self):
        antigen = AlphaNumericAntigenVel(
            given_name="VelSTRONG",
            expressed=True,
            homozygous=False,
            antithetical_relationships={},
        )
        self.assertEqual(antigen.base_name.strip(), "Vel")
        self.assertEqual(antigen.weight, 1)
        self.assertEqual(antigen.name, "VelSTRONG")

    def test_comparisons(self):
        antigen1 = NumericAntigen(
            given_name="3",
            expressed=True,
            homozygous=True,
            antithetical_relationships={},
        )
        antigen2 = NumericAntigen(
            given_name="3w",
            expressed=True,
            homozygous=True,
            antithetical_relationships={},
        )
        self.assertTrue(antigen1 < antigen2)  # 1 < 2 (weight comparison)
        self.assertTrue(antigen2 > antigen1)
        self.assertFalse(antigen1 == antigen2)

    def test_antigen_equality(self):
        antigen1 = AlphaNumericAntigen(
            given_name="K+",
            expressed=True,
            homozygous=True,
            antithetical_relationships={},
        )
        antigen2 = AlphaNumericAntigen(
            given_name="K+",
            expressed=True,
            homozygous=True,
            antithetical_relationships={},
        )
        self.assertTrue(antigen1 == antigen2)

    def test_antigen_repr(self):
        antigen = AlphaNumericAntigen(
            given_name="K+k-",
            expressed=True,
            homozygous=False,
            antithetical_relationships={"K": ["k"]},
        )
        repr_str = repr(antigen)
        self.assertIn("Antigen(name=", repr_str)
        self.assertIn("Antigen(name='K+k-'", repr_str)
        self.assertIn("base_name='Kk'", repr_str)


class TestAlphaNumericAntigenVelWeight(unittest.TestCase):
    def test_weight_strong(self):
        """Test that if 'STRONG' is present (case-insensitive), weight is 1."""
        antigen = AlphaNumericAntigenVel(
            given_name="VelSTRONG",
            expressed=True,
            homozygous=False,
            antithetical_relationships={},
        )
        self.assertEqual(antigen._set_weight(), 1)

    def test_weight_default(self):
        """Test that if no special tokens are present, weight is 2."""
        antigen = AlphaNumericAntigenVel(
            given_name="Vel",
            expressed=True,
            homozygous=False,
            antithetical_relationships={},
        )
        self.assertEqual(antigen._set_weight(), 2)

    def test_weight_weak(self):
        """Test that if '+w' or 'weak' is in the name, weight is 3."""
        antigen_plusw = AlphaNumericAntigenVel(
            given_name="Vel+w",
            expressed=True,
            homozygous=False,
            antithetical_relationships={},
        )
        self.assertEqual(antigen_plusw._set_weight(), 3)

        antigen_weak = AlphaNumericAntigenVel(
            given_name="Velweak",
            expressed=True,
            homozygous=False,
            antithetical_relationships={},
        )
        self.assertEqual(antigen_weak._set_weight(), 3)

    def test_weight_negative(self):
        """Test that if '-' is in the name, weight is 4."""
        antigen = AlphaNumericAntigenVel(
            given_name="Vel-",
            expressed=True,
            homozygous=False,
            antithetical_relationships={},
        )
        self.assertEqual(antigen._set_weight(), 4)



class TestAlphaNumericAntigenABOBaseName(unittest.TestCase):
    def test_get_base_name_A(self):
        """Test base_name returns 'A' when given_name starts with 'A'."""
        antigen = AlphaNumericAntigenABO(
            given_name="A1",
            expressed=True,
            homozygous=False,
            antithetical_relationships={},
        )
        self.assertEqual(antigen.base_name, "A")

    def test_get_base_name_B(self):
        """Test base_name returns 'B' when given_name starts with 'B'."""
        antigen = AlphaNumericAntigenABO(
            given_name="B12",
            expressed=True,
            homozygous=False,
            antithetical_relationships={},
        )
        self.assertEqual(antigen.base_name, "B")

    def test_get_base_name_CIS(self):
        """Test base_name returns 'cis' when given_name starts with 'CIS' (case-insensitive)."""
        antigen = AlphaNumericAntigenABO(
            given_name="CIS-XYZ",
            expressed=True,
            homozygous=False,
            antithetical_relationships={},
        )
        self.assertEqual(antigen.base_name, "cis")

    def test_get_base_name_O(self):
        """Test base_name returns 'O' when given_name starts with 'O'."""
        antigen = AlphaNumericAntigenABO(
            given_name="O3",
            expressed=True,
            homozygous=False,
            antithetical_relationships={},
        )
        self.assertEqual(antigen.base_name, "O")

    def test_get_base_name_invalid(self):
        """Test that an invalid given_name raises a ValueError."""
        with self.assertRaises(ValueError) as context:
            AlphaNumericAntigenABO(
                given_name="X1",
                expressed=True,
                homozygous=False,
                antithetical_relationships={},
            )
        self.assertIn("ABO given name wrong: X1", str(context.exception))



class TestAlphaNumericAntigenName(unittest.TestCase):
    def test_special_case_returns_given_name(self):
        """If given_name contains a special word, property 'name' returns given_name."""
        antigen = AlphaNumericAntigen(
            given_name="Kcommon",
            expressed=True,
            homozygous=True,
            antithetical_relationships={}
        )
        self.assertEqual(antigen.name, "Kcommon")

    def test_not_special_not_weak_expressed_true(self):
        """For non-special, non-weak antigen expressed True, returns base_name + '+'."""
        antigen = AlphaNumericAntigen(
            given_name="K+",
            expressed=True,
            homozygous=True,
            antithetical_relationships={}
        )
        # base_name is "K" so name becomes "K+"
        self.assertEqual(antigen.name, "K+")

    def test_not_special_not_weak_expressed_false(self):
        """For non-special, non-weak antigen expressed False, returns base_name + '-'."""
        antigen = AlphaNumericAntigen(
            given_name="K+",
            expressed=False,
            homozygous=True,
            antithetical_relationships={}
        )
        # base_name is "K" so name becomes "K-"
        self.assertEqual(antigen.name, "K-")

    def test_weak_expressed_true(self):
        """For weak antigen expressed True, returns base_name + '+w'."""
        antigen = AlphaNumericAntigen(
            given_name="K+w",
            expressed=True,
            homozygous=True,
            antithetical_relationships={}
        )
        # base_name is "K" so name becomes "K+w"
        self.assertEqual(antigen.name, "K+w")

    def test_weak_expressed_false(self):
        """For weak antigen expressed False, returns base_name + '-'."""
        antigen = AlphaNumericAntigen(
            given_name="K+w",
            expressed=False,
            homozygous=True,
            antithetical_relationships={}
        )
        # Despite weak state, not expressed returns "K-"
        self.assertEqual(antigen.name, "K-")

    def test_parentheses_expressed_true(self):
        """If given_name contains '(', use alternate generation; expressed True returns modified name."""
        antigen = AlphaNumericAntigen(
            given_name="K(+)",
            expressed=True,
            homozygous=True,
            antithetical_relationships={}
        )
        # base_name is computed via translation: "K(+)" -> "K()"
        # Since '(' in given_name, generate_name returns base_name.replace(")", "+)")
        # "K()" becomes "K(+)"
        self.assertEqual(antigen.name, "K(+)")

    def test_parentheses_expressed_false(self):
        """If given_name contains '(', use alternate generation; expressed False returns modified name."""
        antigen = AlphaNumericAntigen(
            given_name="K(+)",
            expressed=False,
            homozygous=True,
            antithetical_relationships={}
        )
        # Expected: "K()".replace(")", "-)") -> "K(-)"
        self.assertEqual(antigen.name, "K(-)")



if __name__ == "__main__":
    unittest.main()
