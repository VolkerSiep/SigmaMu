# -*- coding: utf-8 -*-
"""This module defines the class :class:`ThermoContribution`, which defines
the building blocks of a :class:`ThermoFrame` function object."""

# stdlib modules
from abc import ABC, abstractmethod
from typing import Collection

# internal modules
from ..utilities import Quantity, ParameterDictionary

QuantityDict = dict[str, Quantity]


class ThermoContribution(ABC):
    """This abstract class defines the interface of a contribution to a
    thermodynamic state function, as collected in :class:`ThermoFrame` objects.

    A contribution is a reusable part of a thermodynamic model that can be
    recombined meaningfully with other contributions. Examples are standard
    state, ideal mix, ideal gas, Gibbs excess contributions, and
    Helmholtz residual functions (equations of state).

    The definition is based on the ``casadi`` library, and its definition
    is required to build a ``casadi`` evaluation structure as the
    implementation of the belonging equations.

    The usage of this class is mainly indirect by instantiation via the
    :cls`ThermoFactory` objects and parametrisation via the provided parameter
    structures.

    A contribution can overwrite the class attribute ``provides``, helping
    the user to identify feasible contributions if a downstream contribution
    does not find a symbol.
    """

    provides: list[str] = []

    def __init__(self, species, options):
        self.species = species
        self.options = options

    @abstractmethod
    def define(self, res: QuantityDict, par: ParameterDictionary):
        """Abstract method to implement the ``casadi`` expressions
        that make up this contribution.

        :param res: A dictionary with already calculated properties that is to
          be supplemented by the properties calculated in this contribution.
          The values of the dictionaries are of type ``casadi.SX``.
        :param par: A dictionary with parameters for this contribution. This
          dictionary can be nested. All values are scalar symbols of type
          ``casadi.SX``.

        .. todo::

            - describe in dedicated section the standard property names

        """

    def relax(self, current_result: dict,
              delta_state: Collection[float]) -> float:
        """Virtual function to report the maximal allowable step size in
        the state variables.

        :param current_result: The numeric results based on the current state
        :param delta_state: The given direction

        :seealso: :meth:`ThermoFrame.relax`
        """
        # default implementation
        del current_result, delta_state  # unused
        return 999  # a number greater than 1 / gamma for practical gamma

    def initial_state(self, temperature: Quantity, pressure: Quantity,
                      quantities: Quantity,
                      properties: dict[str, Quantity]) -> Quantity | None:
        """When the :class:`ThermoFrame` object is queried for an initial state
        representation and deviates from Gibbs coordinates, The upper-most
        contribution that implements this method and does not return ``None``
        takes the responsibility of calculating that state.

        Hence, normally only Helmholtz models need to implement this method.
        The true model coordinates can however be entirely unconventionally,
        such that it is solely up to the contributions on how to obtain the
        initial state.

        :param temperature: Temperature [K]
        :param pressure: Pressure [Pa]
        :param quantities: Quantities [mol]
        :param properties: The property structure, mapping strings to floats
          or list of floats.
        :return: The initial state or ``None``

        .. seealso:: :meth:`ThermoFrame.initial_state`
        """
        del temperature, pressure, quantities, properties  # unused
        return None


class StateDefinition(ABC):
    """This class defines the interpretation of the state vector in terms of
    physical properties. This interpretation is then consumed by the
    contributions as input for their calculations towards the complete
    thermodynamic model."""

    @abstractmethod
    def prepare(self, result: dict):
        """This method can assume to find the state vector ``x`` in the
        ``result`` dictionary, and is expected to add the physical
        interpretation of its elements to the same dictionary. It is entirely
        up to the contributions that rely on this state.

        For the Gibbs state, the new elements would be ``T``, ``p``, and ``n``,
        denoting temperature, pressure and quantities respectively.
        """

    @abstractmethod
    def reverse(self, temperature: Quantity, pressure: Quantity,
                quantities: Quantity) -> list[float]:
        """Return the state vector as complete as possible with given
        temperature, pressure and quantities. The task of the contributions'
        :meth:`ThermoContribution.initial_state` method is it then to
        complete it. Missing elements shall be filled with None.
        """
