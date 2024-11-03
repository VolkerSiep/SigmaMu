"""Module defining classes related to thermodynamic state representation"""
# stdlib modules
from typing import Self
from abc import ABC, abstractmethod
from collections.abc import Sequence
from dataclasses import dataclass

# external modules
from casadi import vertsplit, vertcat

# internal modules
from ..utilities import Quantity, base_unit, base_magnitude
from ..utilities.types import MutMap


@dataclass
class InitialState:
    """Dataclass describing an initial state, which is always defined in terms
    of temperature, pressure, and molar quantities.

    Temperature and pressure are scalar quantities of respective physical
    dimensions. the ``mol_vector`` quantity is vectorial and can arbitrarily
    be defined in units compatible with ``mol`` or ``mol/s``."""

    temperature: Quantity
    pressure: Quantity
    mol_vector: Quantity

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

    @abstractmethod
    def reverse(self, state: InitialState) -> Sequence[float]:
        """Return the state vector as complete as possible with given
        temperature, pressure and quantities. The task of the contributions'
        :meth:`ThermoContribution.initial_state` method is it then to
        complete it. Missing elements shall be filled with None.
        """
        ...
