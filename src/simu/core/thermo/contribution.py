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
from ..utilities import Quantity, ParameterDictionary, QuantityDict
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
        self.reset()

    def reset(self):
        """Reset the object's state by clearing defined residuals."""
        self.__residuals = ResidualHandler()
        self.__parameters = ParameterDictionary()
        self.__bounds = {}
        self.__vector_props = {}

    @property
    def species(self) -> Sequence[str]:
        """Returns a list of species names"""
        return list(self.species_definitions.keys())

    @abstractmethod
    def define(self, res: MutMap[Quantity]):
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

            - refer to dedicated section the standard property names

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

    def add_bound(self, name: str, bound: Quantity,
                  keys: Sequence[str] = None):
        """Add a domain bound to the contribution. This is a property that
        is required to be truly positive.

        If ``bound`` is a vectorial property, its keys must first be registered
        under given ``name`` via :meth:`declare_vector_keys`.

        :param name: Name of the bounded variable
        :param bound: A Scalar or vectorial quantity to remain truly positive
        :param keys: If the quantity is vectorial or ``keys`` are given, they
           will be used to identify the individual element(s). If ``keys`` are
           not given, the vector quantity will use any previously defined
           :meth:`declare_vector_keys` definition instead. If that entry also
           does not exist, species names are assumed to be keys.
        """
        if bound.magnitude.size()[0] == 1 and keys is None:
            self.__bounds[name] = bound
        else:
            if keys is None:
                keys = self.__vector_props.get(name, self.species)
            self.__bounds[name] = QuantityDict.from_vector_quantity(bound, keys)

    @property
    def par_scalar(self):
        """Shortcut method for ``self.parameters.register_scalar``

        .. seealso::
            :class:`~simu.core.utilities.qstructures.ParameterDictionary`
        """
        return self.__parameters.register_scalar

    @property
    def par_vector(self):
        """Shortcut method for ``self.parameters.register_vector``

        .. seealso::
            :class:`~simu.core.utilities.qstructures.ParameterDictionary`
        """
        return self.__parameters.register_vector

    @property
    def par_sparse_matrix(self):
        """Shortcut method for ``self.parameters.register_sparse_matrix``

        .. seealso::
            :class:`~simu.core.utilities.qstructures.ParameterDictionary`
        """
        return self.__parameters.register_sparse_matrix

    @property
    def bounds(self) -> Map[Quantity]:
        return self.__bounds

    @property
    def parameters(self) -> ParameterDictionary:
        return self.__parameters

    @property
    def residuals(self) -> ResidualProxy:
        """Return the defined residuals of this contribution"""
        return self.__residuals

    @property
    def vectors(self) -> Map[Sequence[str]]:
        return self.__vector_props

    def declare_vector_keys(self, name: str, keys: Sequence[str] = None):
        """
        Register a property as a vector property with given keys.
        This way it can be handled correctly by the numeric handler.

        :param name: The name of the property
        :param keys: The names of the keys. If ``None`` (default), the species
          names are used.
        """
        self.__vector_props[name] = self.species if keys is None else keys
