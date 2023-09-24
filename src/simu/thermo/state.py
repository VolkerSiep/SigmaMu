"""Module defining classes related to thermodynamic state representation"""
# stdlib modules
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


class HelmholtzState(StateDefinition):
    """This definition interprets the state as being temperature, volume,
    and mole numbers. Accordingly, it defines:

    ======== ============================
    Property Description
    ======== ============================
    ``T``    Temperature
    ``V``    Volume
    ``n``    Mole vector
    ======== ============================
    """

    def prepare(self, result: MutMap[Quantity], flow: bool = False):
        state = result["_state"].magnitude
        result["T"], result["V"], *n_vec = vertsplit(state, 1)
        result["n"] = vertcat(*n_vec)
        s = "/s" if flow else ""
        for name, unit in [("T", "K"), ("V", f"m**3{s}"), ("n", f"mol{s}")]:
            result[name] = Quantity(result[name], base_unit(unit))

    def reverse(self, state):
        return [base_magnitude(state.temperature), None] + \
            list(base_magnitude(state.mol_vector))


class GibbsState(StateDefinition):
    """This definition interprets the state as being temperature, pressure,
    and mole numbers. Accordingly, it defines:

    ======== ============================
    Property Description
    ======== ============================
    ``T``    Temperature
    ``p``    Pressure
    ``n``    Mole vector
    ======== ============================
    """

    def prepare(self, result: MutMap[Quantity], flow: bool = False):
        state = result["_state"].magnitude
        result["T"], result["p"], *n_vec = vertsplit(state, 1)
        result["n"] = vertcat(*n_vec)
        q_unit = "mol/s" if flow else "mol"
        for name, unit in [("T", "K"), ("p", "Pa"), ("n", q_unit)]:
            result[name] = Quantity(result[name], base_unit(unit))

    def reverse(self, state):
        return [base_magnitude(state.temperature),
                base_magnitude(state.pressure)] + \
                list(base_magnitude(state.mol_vector))
