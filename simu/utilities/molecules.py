"""This module defines a parser for formulae that can be used to establish
generic element balances, and to calculate moecular weights."""

from pathlib import Path
from re import compile as re_compile
from collections import Counter
from yaml import safe_load


def _load_file():
    filename = Path(__file__).resolve().parent / "atomic_weights.yml"
    with open(filename, encoding="utf-8") as file:
        data = safe_load(file)
        return {sym: mw for sym, [_, mw] in data.items()}


ATOMIC_WEIGHTS = _load_file()

__VALID_REG = re_compile(r"^[A-Za-z0-9.()]+$")
__EL_REG = re_compile(r"([A-Z][a-z]*|[(])")
__NUM_REG = re_compile(r"\d+(\.\d*)?")


class MCounter(Counter):
    """This is a slight extention of the ``Collections.Counter`` class
    to also allow multiplication with integers."""

    def __mul__(self, other):
        if not isinstance(other, int):
            raise TypeError("Non-int factor")
        return MCounter({k: other * v for k, v in self.items()})

    def __rmul__(self, other):
        return self * other  # call __mul__

    def __add__(self, other):
        return MCounter(super().__add__(other))

    def __pos__(self):
        return self

    @classmethod
    def fromkeys(cls, iterable, v=None):
        raise NotImplementedError()


__ELEM_COUNTERS = {sym: MCounter({sym: 1}) for sym in ATOMIC_WEIGHTS}


def parse_formula(formula: str) -> dict[str, int]:
    """Parse a formula and return a Counter object with the atomic symbols
    as keys and the number of occurances as values

    .. warning::

        This code uses the evil ``eval`` function that yields potential and
        very real security risks in certain contexts. Here it is considered
        safe for the following reasons.

          - the expression must come from a valid formula and is modified
            into the expression itself.
          - global variable dictionary is set to be empty, and local variables
            are only the element counter objects
          - SiMu is running on the PC of the user. Even if the high level
            models are published via network, the formulae are not part of the
            remote user interface.

    However, If this code was used to offer a web-based formula parser, some
    evil genius might in theory find a way to harm the host system this way.
    """

    def repl_el(arg):
        return f"+{arg.group(0)}"

    def repl_fac(arg):
        return f"*{arg.group(0)}"

    if __VALID_REG.match(formula) is None:
        raise ValueError(f"Invalid syntax of formula: {formula}")
    f_mod = __EL_REG.sub(repl_el, formula)
    f_mod = __NUM_REG.sub(repl_fac, f_mod)

    try:
        return eval(f_mod, {}, __ELEM_COUNTERS)  # pylint: disable=eval-used
    except NameError as err:
        # don't show eval line in trace of error. People could try to be evil.
        raise ValueError(f"Formula contains invalid element: {err}") from None


def molecular_weight(formula: str) -> float:
    """Return the molecular weight for a given formula"""
    elements = parse_formula(formula)
    return sum(ATOMIC_WEIGHTS[sym] * nu for sym, nu in elements.items())
