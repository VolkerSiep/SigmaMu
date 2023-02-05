"""This module defines a parser for formulae that can be used to establish
generic element balances, and to calculate molecular weights."""

from re import compile as re_compile
from yaml import safe_load

from .structures import MCounter
from .quantity import Quantity
from ..data import DATA_DIR


class FormulaParser:
    """This class implements the functionality to analyse chemical sum
    formulae. The atomic composition and molecular weight can be obtained.
    """
    VALID_REG = re_compile(r"^[A-Za-z0-9.()=\-+]+(:\d+[+-])?$")
    EL_REG = re_compile(r"([A-Z][a-z]*|[(])")
    NUM_REG = re_compile(r"\d+(\.\d*)?")
    STRUC_REG = re_compile(r"[=\-+]")
    CHARGE_REG = re_compile(r":\d+[+-]$")

    def __init__(self):
        filename = DATA_DIR / "atomic_weights.yml"
        with open(filename, encoding="utf-8") as file:
            data = safe_load(file)
        self._atomic_weights = {
            sym: Quantity(mw, "g/mol")
            for sym, [_, mw] in data.items()
        }
        self._element_counters = {
            s: MCounter({s: 1})
            for s in self._atomic_weights
        }

    @property
    def atomic_weights(self) -> dict[str, Quantity]:
        """Return a dictionary with all elements mapped to their molecular
        weight"""
        return self._atomic_weights

    def parse(self, formula: str) -> dict[str, int]:
        """Parse a formula and return a Counter object with the atomic symbols
        as keys and the number of occurances as values:

            >>> parser = FormulaParser()
            >>> parser.parse("H3PO4")
            MCounter({'O': 4, 'H': 3, 'P': 1})
            >>> parser.parse("FISH")
            MCounter({'F': 1, 'I': 1, 'S': 1, 'H': 1})
            >>> parser.parse("(NH4)2HPO4")
            MCounter({'H': 9, 'O': 4, 'N': 2, 'P': 1})
            >>> parser.parse("CH3-(CH2)3-CH=O")
            MCounter({'H': 10, 'C': 5, 'O': 1})
            >>> parser.parse("SO4:2-")
            MCounter({'O': 4, 'S': 1})

        .. warning::

            This code uses the evil ``eval`` function that yields potential and
            very real security risks in certain contexts. Here it is considered
            safe for the following reasons.

            - the expression must come from a valid formula and is modified
                into the expression itself.
            - global variable dictionary is set to be empty, and local
                variables are only the element counter objects
            - SiMu is running on the PC of the user. Even if the high level
                models are published via network, the formulae are not part of
                the remote user interface.

        However, If this code was used to offer a web-based formula parser,
        some evil genius might in theory find a way to harm the host system
        this way.
        """

        def repl_el(arg):
            """Replace element ``El`` with ``+El``"""
            return f"+{arg.group(0)}"

        def repl_fac(arg):
            """Replace coefficient ``3`` with ``*3``"""
            return f"*{arg.group(0)}"

        if self.VALID_REG.match(formula) is None:
            raise ValueError(f"Invalid syntax of formula: {formula}")

        # Let's turn the formula into an equation with elements as symbols:
        # E.g. (NH4)2SO4 to (+N+H*4)*2+S+O*4
        # This is done by prepending '+' to each element and '*' to each
        # coefficient.

        f_mod = self.CHARGE_REG.sub("", formula)
        f_mod = self.STRUC_REG.sub("", f_mod)
        f_mod = self.EL_REG.sub(repl_el, f_mod)
        f_mod = self.NUM_REG.sub(repl_fac, f_mod)

        # now it can just be evaluated.
        scope = self._element_counters
        try:
            return eval(f_mod, {}, scope)  # pylint: disable=eval-used
        except NameError as err:
            raise ValueError(
                f"Formula contains invalid element: {err}") from None

    def molecular_weight(self, formula: str) -> Quantity:
        """Return the molecular weight for a given formula. The atomic weights
        are originally obtained from [CoplenShrestha_2016]_.

            >>> parser = FormulaParser()
            >>> mw = parser.molecular_weight("CH3-(CH2)24-CH3")
            >>> print(f"{mw:~.2f}")
            366.72 g / mol

        """
        elements = self.parse(formula)
        weights = self.atomic_weights
        w_0 = Quantity(0, "g/mol")
        return sum((weights[sym] * nu for sym, nu in elements.items()), w_0)

    def charge(self, formula: str) -> int:
        """Return the charge associated to the given formula.

            >>> parser = FormulaParser()
            >>> parser.charge("H2SO4")
            0
            >>> parser.charge("SO4:2-")
            -2
            >>> parser.charge("Al:3+")
            3
            >>> parser.charge("S6:12-")
            -12
        """
        if self.VALID_REG.match(formula) is None:
            raise ValueError(f"Invalid syntax of formula: {formula}")

        match = self.CHARGE_REG.search(formula)
        if match is None:
            return 0
        raw = match.group(0)
        return int(f"{raw[-1]}{raw[1:-1]}")
