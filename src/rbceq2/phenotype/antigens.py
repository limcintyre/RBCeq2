from __future__ import annotations

import operator
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from typing import Callable


@dataclass(slots=True, frozen=False)
class Antigen(ABC):
    """Abstract base class representing an antigen.

    Attributes:
        given_name (str): The original name provided for the antigen.
        expressed (bool): Indicates whether the antigen is expressed.
        homozygous (bool): Indicates whether the antigen is homozygous.
        antithetical_relationships (dict[str, str]): Relationships mapping the base name
            to a string of antithetical antigen names.
        antithetical_antigen (list[Antigen] | None): A list of antithetical antigen objects,
            initialized in __post_init__.
        base_name (str): The processed base name of the antigen, initialized in __post_init__.
        weak (bool): Flag indicating whether the antigen is weak, determined in __post_init__.
        weight (int): A numerical weight used for ranking antigens, set in __post_init__.
    """

    given_name: str
    expressed: bool  # can freeze if not for this...
    homozygous: bool
    antithetical_antigen: list[Antigen] | None = field(init=False)
    base_name: str = field(init=False)
    weak: bool = field(init=False)
    weight: int = field(init=False)
    antithetical_relationships: dict[str, str]

    def __post_init__(self) -> None:
        """Perform post-initialization tasks.

        This method is automatically called after the object is initialized.
        """
        object.__setattr__(self, "base_name", self._get_base_name())
        object.__setattr__(self, "weak", self._is_weak())
        object.__setattr__(self, "weight", self._set_weight())
        object.__setattr__(self, "antithetical_antigen", self._get_antithetical())

    @abstractmethod
    def _is_weak(self):
        """Determine whether the antigen is considered weak.

        Returns:
            bool: True if the antigen is weak, False otherwise.
        """
        pass

    @abstractmethod
    def _get_base_name(self):
        """Extract the base name of the antigen from its given name.

        Returns:
            str: The base name of the antigen.
        """
        pass

    @abstractmethod
    def _set_weight(self):
        """Determine the weight of the antigen based on its properties.

        Returns:
            int: The weight assigned to the antigen.
        """
        pass

    @abstractmethod
    def name(self):
        """Generate or retrieve the formatted antigen name.

        Returns:
            str: The formatted name of the antigen.
        """
        pass

    def _get_antithetical(self) -> list[Antigen] | None:
        """Retrieve antithetical antigens based on the antigen's relationships.

        Uses the antithetical_relationships dictionary to determine if any
        antithetical antigens exist for the current antigen's base name.

        Returns:
            list[Antigen] | None: A list of antithetical Antigen objects if available,
            otherwise None.
        """
        if self.antithetical_relationships is not None:
            names = self.antithetical_relationships.get(self.base_name)
            if names is not None:
                return [
                    type(self)(
                        given_name=name,
                        expressed="-" not in name,
                        homozygous=False,
                        antithetical_relationships={},
                    )
                    for name in names
                ]
        return None

    def _rank(self, other: Antigen, operator_func: Callable[[int, int], bool]) -> bool:
        """Compare the weight of this antigen with another antigen using a given operator.

        Args:
            other (Antigen): Another antigen to compare with.
            operator_func (Callable[[int, int], bool]): A function that compares two
                integers (e.g., operator.gt or operator.lt).

        Returns:
            bool: The result of applying operator_func to the weights of self and other.
        """
        return operator_func(self.weight, other.weight)

    def __gt__(self, other: Antigen) -> bool:
        """Check if this antigen is greater than another antigen based on weight.

        Args:
            other (Antigen): Another antigen to compare.

        Returns:
            bool: True if this antigen's weight is greater than the other antigen's weight,
            False otherwise.
        """
        return self._rank(other, operator.gt)

    def __lt__(self, other: Antigen) -> bool:
        """Check if this antigen is less than another antigen based on weight.

        Args:
            other (Antigen): Another antigen to compare.

        Returns:
            bool: True if this antigen's weight is less than the other antigen's weight,
            False otherwise.
        """
        return self._rank(other, operator.lt)

    def __eq__(self, other: Antigen) -> bool:
        """Check if two antigens are equal based on their weight.

        Args:
            other (Antigen): Another antigen to compare.

        Returns:
            bool: True if the weights of both antigens are equal, False otherwise.
        """
        return self.weight == other.weight

    def __repr__(self) -> str:
        """Generate a string representation of the antigen.

        The representation excludes the antithetical_relationships attribute and includes
        the antigen's type, given name, antithetical antigen (base names), base name, weak flag,
        homozygous flag, expressed flag, and weight.

        Returns:
            str: A formatted string representation of the antigen.
        """
        anti_ant = (
            [ant.base_name for ant in self.antithetical_antigen]
            if self.antithetical_antigen is not None
            else None
        )
        return (
            "\n"
            f"Antigen(type={type(self)!r}, \n"
            f"Antigen(name={self.given_name!r}, \n"
            f"antithetical_antigen={anti_ant!r}, \n"
            f"base_name={self.base_name!r}, \n"
            f"weak={self.weak!r}, \n"
            f"homozygous={self.homozygous!r}, \n"
            f"expressed={self.expressed!r}, \n"
            f"weight={self.weight!r})"
        )


class NumericAntigen(Antigen):
    """A concrete Antigen subclass representing numeric antigens."""

    def _is_weak(self) -> bool:
        """Determine if the NumericAntigen is weak based on its given name.

        Returns:
            bool: True if 'w' or 'weak' is found in the given name, False otherwise.
        """
        return "w" in self.given_name or "weak" in self.given_name

    def _get_base_name(self) -> str:
        """Extract the base name from the given name by removing certain characters.

        Characters removed include '-', '+', and 'w'.

        Returns:
            str: The base name of the NumericAntigen.
        """
        translation_table = str.maketrans("", "", "-+w")
        return self.given_name.translate(translation_table)  # .replace("var", "")

    def _set_weight(self) -> int:
        """Set the weight of the NumericAntigen based on its given name.

        Returns:
            int: The weight of the antigen. A lower number indicates a stronger antigen.
                 Returns 1 if the antigen is strong, 2 if it is weak, and 3 if it has a '-'
                 modifier.
        """
        if (
            "w" not in self.given_name
            and "weak" not in self.given_name
            and "-" not in self.given_name
        ):
            return 1
        elif "w" in self.given_name or "weak" in self.given_name:
            return 2
        elif "-" in self.given_name:
            return 3
        else:
            raise ValueError("weight")

    @property
    def name(self) -> str:
        """Generate the name for the NumericAntigen.

        If the antigen is weak, append a 'w' to the base name. If the antigen is not
        expressed, prepend a '-' to the base name.

        Returns:
            str: The formatted name of the NumericAntigen.
        """
        name = self.base_name if not self.weak else f"{self.base_name}w"

        return name if self.expressed else f"-{name}"


class AlphaNumericAntigen(Antigen):
    """A concrete Antigen subclass representing alphanumeric antigens."""

    def _is_weak(self) -> bool:
        """Determine if the AlphaNumericAntigen is weak based on its given name.

        Returns:
            bool: True if '+w' or 'weak' is present in the given name, False otherwise.
        """
        return "+w" in self.given_name or "weak" in self.given_name

    def _get_base_name(self) -> str:
        """Extract the base name from the given name by removing specific characters.

        Characters removed include '-', '+', and 'w'. Additionally, the substring 'var'
        is removed.

        Returns:
            str: The base name of the AlphaNumericAntigen.
        """
        translation_table = str.maketrans("", "", "-+w")
        return self.given_name.translate(translation_table).replace("var", "")

    def _set_weight(self) -> int:
        """Set the weight of the AlphaNumericAntigen based on its given name.

        Returns:
            int: The weight assigned to the antigen. Returns 1 if strong, 2 if weak, and 3
            if a '-' modifier is present.
        """
        if (
            "+w" not in self.given_name
            and "weak" not in self.given_name
            and "-" not in self.given_name
        ):
            return 1
        elif "+w" in self.given_name or "weak" in self.given_name:
            return 2
        elif "-" in self.given_name:
            return 3
        else:
            raise ValueError("weight")

    @property
    def name(self) -> str:
        """Generate the name for the AlphaNumericAntigen.

        For special cases containing substrings like 'common', 'partial', 'altered', or
        'var', the given name is returned directly. Otherwise, a suffix is generated
        based on whether the antigen is weak and/or expressed.

        Returns:
            str: The formatted name of the AlphaNumericAntigen.
        """
        if any(
            special_case in self.given_name.lower()
            for special_case in ["common", "partial", "altered", "var"]
        ):
            return self.given_name

        def generate_name(suffix: str) -> str:
            """Helper function to generate the name with a given suffix.

            Args:
                suffix (str): The suffix to append or replace in the base name.

            Returns:
                str: The formatted allele name.
            """
            return (
                f"{self.base_name}{suffix}"
                if "(" not in self.given_name
                else self.base_name.replace(")", f"{suffix})")
            )

        # Determine the suffix based on the allele's state
        if self.weak:
            name_expressed = generate_name("+w")
            name_not_expressed = generate_name("-")
        else:
            name_expressed = generate_name("+")
            name_not_expressed = generate_name("-")

        # Return the appropriate name based on whether the allele is expressed
        return name_expressed if self.expressed else name_not_expressed


class AlphaNumericAntigenABO(AlphaNumericAntigen):
    """An AlphaNumericAntigen subclass for ABO blood group antigens."""

    @property
    def name(self) -> str:
        """Return the name for the AlphaNumericAntigenABO.

        For ABO antigens, the given name is returned directly without modification.

        Returns:
            str: The original given name.
        """
        return self.given_name

    def _get_base_name(self) -> str:
        """Determine the base name for an ABO antigen based on the given name.

        The base name is determined by checking the starting character(s) of the
        given name. It must start with 'A', 'B', 'CIS', or 'O'. If none match, a
        ValueError is raised.

        Returns:
            str: The base name for the ABO antigen.

        Raises:
            ValueError: If the given name does not start with a recognized ABO type.
        """
        if self.given_name.upper().startswith("A"):
            return "A"
        elif self.given_name.upper().startswith("B"):
            return "B"
        elif self.given_name.upper().startswith("CIS"):
            return "cis"
        elif self.given_name.upper().startswith("O"):
            return "O"
        else:
            raise ValueError(f"ABO given name wrong: {self.given_name}")


class AlphaNumericAntigenXG(AlphaNumericAntigen):
    """An AlphaNumericAntigen subclass for XG blood group antigens."""

    @property
    def name(self) -> str:
        """Generate the name for the AlphaNumericAntigenXG.

        Returns:
            str: The given name with parentheses and plus signs removed, with 'CD99'
            appended.
        """
        return (
            self.given_name.replace("(", "").replace(")", "").replace("+", "") + "CD99"
        )


class AlphaNumericAntigenMNS(AlphaNumericAntigen):
    """An AlphaNumericAntigen subclass for MNS blood group antigens."""

    def _get_base_name(self) -> str:
        """Extract the base name for the AlphaNumericAntigenMNS by removing specific
        characters and substrings.

        Characters removed include '-', '+', and 'w'. Additionally, substrings like
        'partial', 'alteredGPB', 'alteredU/GPB', and 'var' are removed. The result is
        then stripped of leading and trailing whitespace.

        Returns:
            str: The base name for the MNS antigen.
        """
        translation_table = str.maketrans("", "", "-+w")
        return (
            self.given_name.translate(translation_table)
            .replace("partial", "")
            .replace("alteredGPB", "")
            .replace("alteredU/GPB", "")
            .replace("var", "")
            .strip()
        )


class NumericAntigenVel(NumericAntigen):
    def _set_weight(self) -> int:
        """Set the weight for the AlphaNumericAntigenVel based on its given name.

        For Vel antigens, if 'STRONG' is present in the given name (case-insensitive),
        the weight is set to 1. Otherwise, the weight is determined by the presence of
        '+w', 'weak', or '-' modifiers.

        Returns:
            int: The weight assigned to the Vel antigen.
        """
        if "s" in self.given_name:
            return 1
        elif "+w" not in self.given_name and "-" not in self.given_name:
            return 2
        elif "+w" in self.given_name:
            return 3
        elif "-" in self.given_name:
            return 4
        else:
            raise ValueError("Vel weight")

    def _get_base_name(self) -> str:
        """Extract the base name from the given name by removing certain characters.

        Characters removed include '-', '+', and 'w' etc.
        weak = w, strong = s

        Returns:
            str: The base name of the NumericAntigen.
        """
        translation_table = str.maketrans("", "", "-+ws")
        return self.given_name.translate(translation_table)


class AlphaNumericAntigenVel(AlphaNumericAntigen):
    """An AlphaNumericAntigen subclass for Vel blood group antigens."""

    def _get_base_name(self) -> str:
        """Extract the base name for the AlphaNumericAntigenVel by removing specific
        characters and substrings.

        Characters removed include '-', '+', and 'w'. Additionally, substrings like 'var',
        'strong', and 'STRONG' are removed.

        Returns:
            str: The base name for the Vel antigen.
        """
        translation_table = str.maketrans("", "", "-+w")
        return (
            self.given_name.translate(translation_table)
            .replace("var", "")
            .replace("strong", "")
            .replace("STRONG", "")
        )

    def _set_weight(self) -> int:
        """Set the weight for the AlphaNumericAntigenVel based on its given name.

        For Vel antigens, if 'STRONG' is present in the given name (case-insensitive),
        the weight is set to 1. Otherwise, the weight is determined by the presence of
        '+w', 'weak', or '-' modifiers.

        Returns:
            int: The weight assigned to the Vel antigen.
        """
        if "STRONG" in self.given_name.upper():
            return 1
        elif (
            "+w" not in self.given_name
            and "weak" not in self.given_name
            and "-" not in self.given_name
        ):
            return 2
        elif "+w" in self.given_name or "weak" in self.given_name:
            return 3
        elif "-" in self.given_name:
            return 4
        else:
            raise ValueError("vel weight")

    @property
    def name(self) -> str:
        """Return the name for the AlphaNumericAntigenVel.

        For Vel antigens, the given name is returned directly without modification.

        Returns:
            str: The original given name.
        """
        return self.given_name


class AlphaNumericAntigenRHCE(AlphaNumericAntigen):
    """An AlphaNumericAntigen subclass for RHCE blood group antigens.
    neg = n, partial = p, weak = w, monoclonal = m, infered = i
    RH: -2,-3,4,5,6,-7 ,-8 ,-9, -10,-11,-12,17 ,18,19 ,-20,-21,-22,-23,26    ,
        C-,E-,c,e,f,Ce-,Cw-,Cx-, V-,Ew-, G-,Hro,Hr,hrS,VS-,CG-,CE-,Dw-,c-like,
        -27,-28 ,29  ,-30 ,31 ,-32  ,-33 , 34 ,-35 , -36, -37   ,-39  ,-40 ,-41 ,
        cE-,hrH-,Rh29,Goa-,hrB,Rh32-,Rh33-,HrB,Rh35-,Bea-,Evans-,Rh39-,Tar-,Rh41-,
        -42 , -43      ,44 ,-45, 46 ,47 ,-48,-49   ,-50  ,51 ,-52 , -53 , -54,
        Rh42-,Crawford-,Nou,Riv-,Sec,Dav,JAL-,STEM-,FPPT-,MAR,BARC-,JAHK-,DAK-,
        -55 , -56 , 57  ,58  ,59  ,-60 , 61  ,62  ,-63
        LOCR-,CENR-,CEST,CELO,CEAG,PARG-,CEVF,CEWA,CETW-
    """

    def _is_weak(self) -> bool:
        """Determine if the AlphaNumericAntigen is weak based on its given name.

        Returns:
            bool: True if 'weak' is present in the given name, False otherwise.
        """
        return (
            "weak" in self.given_name.lower() or "very_weak" in self.given_name.lower()
        )

    # def _is_partial(self):
    #     """Determine if the AlphaNumericAntigen is weak based on its given name.

    #     Returns:
    #         bool: True if 'partial' is present in the given name, False otherwise.
    #     """
    #     return "partial" in self.given_name.lower()

    # def _is_negative(self):
    #     """Determine if the AlphaNumericAntigen is negative based on its given name.

    #     Returns:
    #         bool: True if 'negative' or 'neg' is present in the given name, False otherwise.
    #     """
    #     return "negative" in self.given_name.lower() or "neg" in self.given_name.lower()

    # def _is_monoclonal(self):
    #     """Determine if the AlphaNumericAntigen is monoclonal based on its given name.

    #     Returns:
    #         bool: True if 'monoclonal' is present in the given name, False otherwise.
    #     """
    #     return "monoclonal" in self.given_name.lower()

    # def _is_infered(self):
    #     """Determine if the AlphaNumericAntigen is infered based on its given name.

    #     Returns:
    #         bool: True if 'infered' is present in the given name, False otherwise.
    #     """
    #     return "infered" in self.given_name.lower()

    # is robust
    # these are just for repr??

    def _get_base_name(self) -> str:
        """Extract the base name for the AlphaNumericAntigenRHCE by removing specific
        characters and substrings:
            partial
            robust
            expression
            weak
            Some
            monoclonal
            anti-D
            cross-
            react
            cross-react
            very
            neg


        Returns:
            str: The base name for the RHCE antigen.
        """
        translation_table = str.maketrans("", "", "-+_")
        return (
            self.given_name.translate(translation_table)
            .replace("partial", "")
            .replace("robust", "")
            .replace("expression", "")
            .replace("weak", "")
            .replace("some", "")
            .replace("anti_d", "")
            .replace("cross_", "")
            .replace("react", "")
            .replace("cross_react", "")
            .replace("very", "")
            .replace("neg", "")
            .replace("as", "")
            .replace("probable", "")
            .replace("trans", "")
            .replace("positive", "")
            .replace("in", "")
        )

    def _set_weight(self) -> int:
        """
        Set the weight for the AlphaNumericAntigenRHCE based on its given name.
        Lower weight = stronger/more expressed.
        1: Robust ("ROBUST")
        2: Normal/Strong (no other specific flags)
        3: Partial ("PARTIAL")
        4: Weak ("WEAK")
        5: Weak Partial ("PARTIAL" and "WEAK")
        6: Very Weak ("VERY_WEAK")
        7: Negative ("NEG", including "WEAK TO NEG", "VERY_WEAK TO NEG")
        8: Not Expressed/Null ("-")
        """
        name_upper = self.given_name.upper()

        # Highest priority: "NEG"
        if "NEG" in name_upper:  # This catches "WEAK TO NEG", "PARTIAL NEG", etc.
            return 7
        # Second highest priority: "-" for not expressed/null
        if "-" in self.given_name:  # Check the original string for the '-' character
            return 8

        # Check for specific characteristics
        is_robust = "ROBUST" in name_upper
        is_partial = "PARTIAL" in name_upper
        is_weak = "WEAK" in name_upper  # Note: "VERY_WEAK" also contains "WEAK"
        is_very_weak = "VERY_WEAK" in name_upper

        if is_robust:
            return 1

        # Case 1: It's "PARTIAL" and also some form of "WEAK"
        if is_partial:
            if is_weak:  # And not is_very_weak (because that implies "VERY_WEAK")
                return 5  # Weak Partial
            else:  # Just Partial, not weak
                return 3  # Partial

        # Case 2: Not "PARTIAL" or "PARTIAL" was handled
        if (
            is_very_weak
        ):  # This will catch "VERY_WEAK" and "PARTIAL VERY_WEAK" if not handled above
            return 6  # Very Weak

        if is_weak:  # And not is_very_weak (because that would have been caught)
            return 4  # Weak

        if (
            is_partial
        ):  # This means it's "PARTIAL" alone (no weak, no very_weak, no robust)
            return 3  # Partial

        return 2  # Normal / Strong

    @property
    def name(self) -> str:
        """Return the name for the AlphaNumericAntigenVel.

        For RHCE antigens, the given name is returned directly without modification.

        Returns:
            str: The original given name.
        """
        return self.given_name


class NumericAntigenRHCE(NumericAntigen):
    def _is_weak(self) -> bool:
        """Determine if the NumericAntigen (RHCE) is weak.
        Considers 'w' (weak) or 'v' (very weak).
        """
        name_lower = self.given_name.lower()
        return "w" in name_lower or "v" in name_lower

    def _set_weight(self) -> int:
        """
        Set the weight of the NumericAntigenRHCE based on its given name.
        Lower weight = stronger/more expressed.
        Suffixes: r (robust), p (partial), w (weak), v (very_weak), n (negative).
                  '-' for not expressed/null.

        1: Robust ('r')
        2: Normal/Strong (no other specific flags)
        3: Partial ('p')
        4: Weak ('w')
        5: Weak Partial ('p' and 'w')
        6: Very Weak ('v')
        7: Negative ('n', including if combined with w, v, p)
        8: Not Expressed/Null ('-')
        """
        name_lower = self.given_name.lower()  # Use lower for char checks

        # Highest priority: 'n' for negative
        if "n" in name_lower:
            return 7
        # Second highest priority: '-' for not expressed/null
        if "-" in self.given_name:  # check original string for '-'
            return 8

        # Check for specific characteristic flags
        is_robust = "r" in name_lower
        is_partial = "p" in name_lower
        is_weak = "w" in name_lower
        is_very_weak = "v" in name_lower

        if is_robust:
            return 1

        if is_partial:
            if is_weak:  # 'p' and 'w' are present
                return 5  # Weak Partial
            else:  # Just 'p'
                return 3  # Partial

        if is_very_weak:  # 'v' is present (could be 'v' alone or 'vp')
            return 6  # Very Weak

        if is_weak:  # 'w' is present (could be 'w' alone, 'pw' was handled)
            return 4  # Weak

        if is_partial:  # 'p' is present alone (other 'p' combos handled)
            return 3  # Partial

        return 2  # Normal / Strong

    def _get_base_name(self) -> str:
        """Extract the base name from the given name by removing certain characters.

        Characters removed include '-', '+', and 'w' etc.
        neg = n, partial = p, weak = w, very_weak = v, robust = r

        Returns:
            str: The base name of the NumericAntigen.
        """
        translation_table = str.maketrans("", "", "-+wpimnrv")
        return self.given_name.translate(translation_table)


class AlphaNumericAntigenRHD(AlphaNumericAntigen):
    """An AlphaNumericAntigen subclass for RHD blood group antigens."""

    def _is_weak(self) -> bool:
        """Determine if the AlphaNumericAntigen is weak based on its given name.

        Returns:
            bool: True if 'weak' is present in the given name, False otherwise.
        """
        return "weak" in self.given_name.lower()

    def _get_base_name(self) -> str:
        """Extract the base name for the AlphaNumericAntigenRHD

        Returns:
            str: D
        """
        return "D"

    def _set_weight(self) -> int:
        """Set the weight for the AlphaNumericAntigenRHD based on its given name."""
        if "-" in self.given_name:
            return 4  # null
        elif "Del" in self.given_name or "WEAK" in self.given_name.upper():
            return 3  # very weak
        elif "Type" in self.given_name or "WEAK" in self.given_name.upper():
            return 2  # weak
        else:
            return 1

    @property
    def name(self) -> str:
        """Return the name for the AlphaNumericAntigenRHD.

        For RHD antigens, the given name is returned directly without modification.

        Returns:
            str: The original given name.
        """
        return self.given_name
