# -*- coding: utf-8 -*-
"""This module defines the class :class:`ThermoContribution`, which defines
the building blocks of a :class:`ThermoFrame` function object."""

# stdlib modules
from abc import ABC, abstractmethod
from collections.abc import Sequence, MutableSequence
from typing import Any

# internal modules
from .species import SpeciesDefinition
from .state import InitialState
from ..utilities import Quantity, ParameterDictionary
from ..utilities.residual import ResidualHandler, ResidualProxy
from ..utilities.types import Map, MutMap


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
    :class:`ThermoFactory` objects and parametrisation via the provided parameter
    structures.

    A contribution can overwrite the class attribute ``provides``, helping
    the user to identify feasible contributions if a downstream contribution
    does not find a symbol.
    """

    provides: list[str] = []

    species_definitions: Map[SpeciesDefinition]
    """The map of species definition objects"""

    options: Map[Any]
    """The map of species definition objects"""

    def __init__(self, species: Map[SpeciesDefinition], options):
        self.species_definitions: Map[SpeciesDefinition] = species
        self.options = options
        self.__residuals = ResidualHandler()

    @property
    def species(self) -> Sequence[str]:
        """Returns a list of species names"""
        return list(self.species_definitions.keys())

    @abstractmethod
    def define(self, res: MutMap[Quantity],
               bounds: MutMap[Quantity],
               par: ParameterDictionary):
        """Abstract method to implement the ``casadi`` expressions
        that make up this contribution.

        See :ref:`standard property names` for a guideline on how to name
        standard properties generated in the contribution implementations.

        :param res: A dictionary with already calculated properties that is to
          be supplemented by the properties calculated in this contribution.
          The values of the dictionaries are of type ``casadi.SX``.
        :param bounds: A dictionary including properties of which the base_unit
          magnitude must stay positive. Solvers can use this information to
          stay within the mathematical domain of the model.
          By convention, if the property is also a result, the same name shall
          be used, and it is not a problem if prior entries are over-written.
          For instance multiple contributions will not allow negative ``T``,
          and all of them shall declare this, as they cannot rely on the others
          being used in the same model.
        :param par: A dictionary with parameters for this contribution. This
          dictionary can be nested. All values are scalar symbols of type
          ``casadi.SX``.

        .. todo::

            - describe in dedicated section the standard property names

        """

    def initial_state(self, state: InitialState, properties: Map[Quantity]) \
            -> MutableSequence[float] | None:
        """When the :class:`ThermoFrame` object is queried for an initial state
        representation and deviates from Gibbs coordinates, The uppermost
        contribution that implements this method and does not return ``None``
        takes the responsibility of calculating that state.

        Hence, normally only Helmholtz models need to implement this method.
        The true model coordinates can however be entirely unconventionally,
        such that it is solely up to the contributions, how to obtain the
        initial state.

        :param state: The default state
        :param properties: The property structure, mapping strings to floats
          or list of floats.
        :return: The initial state or ``None``

        .. seealso:: :meth:`ThermoFrame.initial_state`
        """
        ...

    def add_residual(self, name: str, residual: Quantity,
                     tol_unit: str, tol: float = 1e-7):
        """Define a residual that represents an implicit constraint in the
        thermodynamic model itself. Typical examples are equilibrium
        constraints on apparent species systems and any implicit thermodynamic
        models."""
        self.__residuals.add(name, residual, tol_unit,  tol)

    @property
    def residuals(self) -> ResidualProxy:
        """Return the defined residuals of this contribution"""
        return self.__residuals

    def declare_vector_keys(self) -> Map[Sequence[str]]:
        """Declare the keys of newly introduced vectorial properties. In most
        cases, this will be the species names for the mole vector, and the
        implementation will look as follows::

            def declare_vector_keys(self):
                return {"my_vector_property": self.species}

        """
        return {}