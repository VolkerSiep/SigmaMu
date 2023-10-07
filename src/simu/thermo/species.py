"""This module defines data structures to host the global species list"""

from dataclasses import dataclass, field
from collections.abc import Mapping, Iterator

from ..utilities import Quantity
from ..utilities.molecules import FormulaParser

_PARSER = FormulaParser()

@dataclass
class SpeciesDefinition:
    """
    This class holds a definition of a species and provides, based on the
    formula, the molecular weight, charge, and a dictionary of element
    composition.

    >>> a = SpeciesDefinition("H3PO4")
    >>> print(f"{a.molecular_weight:~.3f}")
    97.993 g / mol
    >>> a.elements
    {'H': 3, 'P': 1, 'O': 4}
    >>> a = SpeciesDefinition("PO4:3-")
    >>> a.charge
    -3
    """
    formula: str
    molecular_weight: Quantity = field(init=False)
    charge: int = field(init=False)
    elements: dict[str, int] = field(init=False)

    def __post_init__(self):
        self.elements = dict(_PARSER.parse(self.formula))
        self.molecular_weight = _PARSER.molecular_weight(self.formula)
        self.charge = _PARSER.charge(self.formula)


class SpeciesDB(Mapping[str, SpeciesDefinition]):
    """Based on a dictionary of species names to formulae, this class
    represents a dictionary of the species names to species definitions.

    .. note::
        That meeting could have been an email, and this class could have
        been a function!?
    """
    def __init__(self, formulae: Mapping[str, str]):
        self.__species = {n: SpeciesDefinition(f) for n, f in formulae.items()}

    def __getitem__(self, key: str) -> SpeciesDefinition:
        return self.__species[key]

    def __len__(self) -> int:
        return len(self.__species)

    def __iter__(self) -> Iterator[str]:
        return iter(self.__species)
