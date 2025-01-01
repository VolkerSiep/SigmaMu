"""Module defining classes related to thermodynamic state representation"""
# stdlib modules
from typing import Self
from abc import ABC, abstractmethod
from collections.abc import Sequence
from dataclasses import dataclass

from pint import DimensionalityError

# internal modules
from ..utilities import Quantity, qvertcat
from ..utilities.types import Map, NestedMap
from .species import SpeciesDefinition

@dataclass
class InitialState:
    """Dataclass describing an initial state, which is always defined in terms
    of temperature, pressure, and molar quantities.

    Temperature and pressure are scalar quantities of respective physical
    dimensions. the ``mol_vector`` quantity is vectorial and can arbitrarily
    be defined in units compatible with ``mol`` or ``mol/s``.
    The stored value will however be as a state, not a flow, i.e. compatible to
    ``mol``.
    """

    temperature: Quantity
    pressure: Quantity
    mol_vector: Quantity

    def __post_init__(self):
        """Check that attributes have compatible units, and convert mole flows
        to mole quantities (per 1 second).
        Raises ``DimensionalityError`` if dimensions do not match.
        """
        self.temperature.to("K")  # compatible unit?
        self.pressure.to("Pa")  # compatible unit?
        try:
            self.mol_vector.to("mol")
        except DimensionalityError:
            n = self.mol_vector.to("mol/s")
            self.mol_vector =  Quantity(n.magnitude, "mol")
        self.__setattr__ = self.__setattr_prototype  # freeze this object

    def __setattr_prototype(self, name, value):
        raise AttributeError(f"Cannot modify frozen instance: {name}")

    @classmethod
    def from_si(cls, temperature: float, pressure: float,
                mol_vector: Sequence[float]) -> Self:
        """Construct an initial state based on SI units, i.e. K, Pa and mol."""
        return cls(temperature=Quantity(temperature, "K"),
                   pressure=Quantity(pressure, "Pa"),
                   mol_vector=Quantity(mol_vector, "mol"))

    @classmethod
    def from_cbar(cls, temperature: float, pressure: float,
                  mol_vector: Sequence[float]) -> Self:
        """Construct an initial state based on degC, bar and mol as units."""
        return cls(temperature=Quantity(temperature, "degC"),
                   pressure=Quantity(pressure, "bar"),
                   mol_vector=Quantity(mol_vector, "mol"))

    @classmethod
    def from_std(cls, num_species: int):
        """Construct an initial state at 25 degC, 1 bar and one mol for each
        species."""
        return cls(temperature=Quantity(25, "degC"),
                   pressure=Quantity(1, "atm"),
                   mol_vector=Quantity([1.0] * num_species, "mol"))

    @classmethod
    def from_dict(cls, struct: NestedMap[Quantity], species: Sequence[str]):
        """Convert a nested Quantity map as obtained by :meth:`to_dict` to an
        initial state. The map must contain the top-level keys ``T``, ``p`` and
        ``n``, and a mapping from species to molar quantities as value of the
        ``n`` key. The correct sequence of ``species`` is provided separately,
        to be consistent with the targeted thermodynamic model.

        The method raises a ``KeyError``, if the elements of ``struct`` are not
        as expected, in particular considering the given ``species`` argument.
         """
        for s in struct["n"].keys():
            if s not in species:
                raise KeyError(f"Species {s} not expected")
        mole_list = [struct["n"][s] for s in species]
        units = mole_list[0].units
        mole_vector = [m.magnitude for m in mole_list]
        return cls(temperature=struct["T"],
                   pressure=struct["p"],
                   mol_vector=Quantity(mole_vector, units))


    def to_dict(self, species: Sequence[str]) -> NestedMap[Quantity]:
        """Convert initial state into dictionary of quantities.
        The mole vector is described as a sub-dictionary, mapping species
        to partial molar quantities"""
        mag = self.mol_vector.magnitude
        unit = self.mol_vector.units
        if len(species) != mag.size:
            raise ValueError("Incompatible length of species list")
        moles = {s: Quantity(mag[i], unit) for i, s in enumerate(species)}
        return {"T": self.temperature, "p": self.pressure, "n": moles}


class StateDefinition(ABC):
    """This class defines the interpretation of the state vector in terms of
    physical properties. This interpretation is then consumed by the
    contributions as input for their calculations towards the complete
    thermodynamic model."""

    @abstractmethod
    def prepare(self, result: dict, flow: bool = False):
        """This method can assume to find the state vector ``x`` in the
        ``result`` dictionary, and is expected to add the physical
        interpretation of its elements to the same dictionary. It is entirely
        up to the contributions that rely on this state.

        For the Gibbs state, the new elements would be ``T``, ``p``, and ``n``,
        denoting temperature, pressure and quantities respectively.

        The parameter ``flow`` impacts the definition of units of measurement.
        When ``True``, all extensive variables are divided by time.
        """
        ...

    def declare_vector_keys(
            self, species: Map[SpeciesDefinition]) -> Map[Sequence[str]]:
        """Declare the keys of vectorial state properties. In most cases,
        this will be the species names for the mole vector."""
        return {}

    @abstractmethod
    def reverse(self, state: InitialState) -> Sequence[float]:
        """Return the state vector as complete as possible with given
        temperature, pressure and quantities. The task of the contributions'
        :meth:`ThermoContribution.initial_state` method is it then to
        complete it. Missing elements shall be filled with None.
        """
        ...
