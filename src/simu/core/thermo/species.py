"""This module defines data structures to host the global species list"""

from typing import Self
from dataclasses import dataclass, field
from collections.abc import Mapping, Iterator, Sequence

from ..utilities import Quantity
from ..utilities.structures import Map
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
    >>> print(f"{a.charge:~}")
    -3 e / mol
    """
    formula: str
    """The formula as it was given in the constructor. The admitted formula 
    syntax is described with examples for the
    :class:`FormulaParser <simu.core.utilities.molecules.FormulaParser>` class.
    """
    molecular_weight: Quantity = field(init=False)
    """The molecular weight determined by summing up the atomic weights of 
    the contained atoms; quantum effects and electron masses are neglected."""
    charge: Quantity = field(init=False)
    """The electronic charge of the species' molecule"""
    elements: Map[int] = field(init=False)
    """A dictionary, mapping the occurring atoms to their amount in the 
    species' molecule"""

    def __post_init__(self):
        self.elements = dict(_PARSER.parse(self.formula))
        self.molecular_weight = _PARSER.molecular_weight(self.formula)
        self.charge = _PARSER.charge(self.formula)


class SpeciesDB(Mapping[str, SpeciesDefinition]):
    """Based on a dictionary of species names to formulae, this class
    represents a dictionary of the species names to species definitions.

    .. note::
        For now, this class is quite primitive, but might be extended to handle
        more meta-data, such as CAS registry numbers and species aliases.
    """
    def __init__(self, formulae: Mapping[str, str]):
        """Create a species collection based on a mapping of species names
        to their formulae."""
        self.__species = {n: SpeciesDefinition(f) for n, f in formulae.items()}

    def __getitem__(self, key: str) -> SpeciesDefinition:
        return self.__species[key]

    def __len__(self) -> int:
        return len(self.__species)

    def __iter__(self) -> Iterator[str]:
        return iter(self.__species)


    def get_sub_db(self, keys: Sequence[str]) -> Self:
        """Get a subset of species definitions as a new SpeciesDB object."""
        return {n: self[n] for n in keys}
