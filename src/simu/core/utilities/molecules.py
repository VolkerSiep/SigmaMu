"""This module defines a parser for formulae that can be used to establish
generic element balances, and to calculate molecular weights."""

# stdlib
from re import compile as re_compile
from yaml import safe_load
from ast import (AST, expr, parse, BinOp, UnaryOp, Constant, Name, dump,
                 Add, Mult, UAdd)
from operator import add, mul
from typing import Type, Callable

# internal
from simu.core.data import DATA_DIR
from .structures import MCounter, Map
from .quantity import Quantity


class FormulaParser:
    """This class implements the functionality to analyse chemical sum
    formulae. The atomic composition and molecular weight can be obtained.
    """
    VALID_REG = re_compile(r"^[A-Za-z0-9.·\[{(<>)}\]=≡|\-+]+(:\d+[+-])?$")
    EL_REG = re_compile(r"([A-Z][a-z]*|[(])")
    NUM_REG = re_compile(r"\d+(\.\d*)?")
    CRYSTAL_REG = re_compile(r"·(\d+)?([^·]+)")
    STRUC_REG = re_compile(r"[=\-<>≡|+]")
    CHARGE_REG = re_compile(r":\d+[+-]$")
    OPERATORS: dict[Type, Callable] = {Add: add, Mult: mul}

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
    def atomic_weights(self) -> Map[Quantity]:
        """Return a dictionary with all elements mapped to their molecular
        weight"""
        return self._atomic_weights

    def parse(self, formula: str) -> MCounter:
        """Parse a formula and return a Counter object with the atomic symbols
        as keys and the number of occurances as values:

        Plain formulae:
            >>> parser = FormulaParser()
            >>> parser.parse("H3PO4")
            MCounter({'O': 4, 'H': 3, 'P': 1})
            >>> parser.parse("KMnO4")
            MCounter({'O': 4, 'K': 1, 'Mn': 1})
            >>> parser.parse("FISH")
            MCounter({'F': 1, 'I': 1, 'S': 1, 'H': 1})

        With parantheses:
            >>> parser.parse("(NH4)2HPO4")
            MCounter({'H': 9, 'O': 4, 'N': 2, 'P': 1})

        With structure:
            >>> parser.parse("CH3-(CH2)3-CH=O>")
            MCounter({'H': 10, 'C': 5, 'O': 1})
            >>> parser.parse("|N≡N|")
            MCounter({'N': 2})
            >>> parser.parse("<O=O>")
            MCounter({'O': 2})

        With charge:
            >>> parser.parse("SO4:2-")
            MCounter({'O': 4, 'S': 1})

        With a complex:
            >>> parser.parse("Na(UO2)3[Zn(H2O)6](CH3CO2)9")
            MCounter({'H': 39, 'O': 30, 'C': 18, 'U': 3, 'Na': 1, 'Zn': 1})

        With crystal water:
            >>> parser.parse("CuSO4·5H2O")
            MCounter({'H': 10, 'O': 9, 'Cu': 1, 'S': 1})

        With crystal water and another solvent:
            >>> parser.parse("CuSO4·3H2O·2(CH3)-COOH")
            MCounter({'H': 14, 'O': 11, 'C': 4, 'Cu': 1, 'S': 1})
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
        for c in "{( [( }) ])".split():
            f_mod =  f_mod.replace(c[0], c[1])
        f_mod = self.CRYSTAL_REG.sub(r"(\2)\1", f_mod)
        f_mod = self.STRUC_REG.sub("", f_mod)
        f_mod = self.EL_REG.sub(repl_el, f_mod)
        f_mod = self.NUM_REG.sub(repl_fac, f_mod)

        # now it can just be evaluated.
        scope = self._element_counters
        try:
            return self._evaluate(f_mod, scope)
        except NameError as err:
            raise ValueError(
                f"Formula contains invalid element: {err}") from None

    def molecular_weight(self, formula: str) -> Quantity:
        """Return the molecular weight for a given formula. The atomic weights
        are originally obtained from :cite:p:`CoplenShrestha_2016`.

            >>> parser = FormulaParser()
            >>> mw = parser.molecular_weight("CH3-(CH2)24-CH3")
            >>> print(f"{mw:~.2f}")
            366.72 g / mol

        """
        elements = self.parse(formula)
        weights = self.atomic_weights
        w_0 = Quantity(0, "g/mol")
        return sum((weights[sym] * nu for sym, nu in elements.items()), w_0)

    def charge(self, formula: str) -> Quantity:
        """Return the charge associated to the given formula.

            >>> parser = FormulaParser()
            >>> for n in "H2SO4 SO4:2- Al:3+ S6:12-".split():
            ...     print(f"{n}: {parser.charge(n):~}")
            H2SO4: 0 e / mol
            SO4:2-: -2 e / mol
            Al:3+: 3 e / mol
            S6:12-: -12 e / mol
        """
        if self.VALID_REG.match(formula) is None:
            raise ValueError(f"Invalid syntax of formula: {formula}")

        match = self.CHARGE_REG.search(formula)
        if match is None:
            return Quantity(0, "e / mol")
        raw = match.group(0)
        return Quantity(int(f"{raw[-1]}{raw[1:-1]}"), "e / mol")

    @staticmethod
    def _evaluate(expression: str, variables: Map[MCounter]) -> MCounter:
        def _eval(node: AST | expr):
            match node:
                case BinOp(left=left, op=op, right=right):
                    ops = FormulaParser.OPERATORS
                    return ops[type(op)](_eval(left), _eval(right))
                case UnaryOp(op=UAdd(), operand=child):
                    return _eval(child)
                case Constant(value=value):
                    return value
                case Name(id=name):
                    return variables[name]
            raise ValueError(f"Unsupported node: {dump(node)}")

        tree = parse(expression, mode="eval")
        return _eval(tree.body)