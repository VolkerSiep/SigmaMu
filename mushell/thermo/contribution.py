# -*- coding: utf-8 -*-

"""This module defines the class :class:`ThermoContribution`, which defines
the building blocks of a :class:`ThermoFrame` function object."""

from abc import ABC, abstractmethod
from typing import Collection, List

# external modules
from casadi import SX, vertcat


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

    The static attributes ``name``, ``category``, and ``requires`` are
    introduced to prevent invalid configurations, such as incompatible
    contributions in one model.

    .. todo::
        - refer to general section about categories
        - write that section
    """

    name = ""
    """The name of the contribution (str), used to define direct
    dependencies. The name only needs to be unique within all contributions
    within the same category. If the category is unique in itself, the name
    can remain empty."""

    category = ""
    """The category of the contribution (str), used to define general
    dependencies"""

    requires = []
    """The dependencies as (name, category) pairs or single category
    entries. If elements of ``requires`` are of type ``str``, any contribution
    of that category is accepted. If an element is a list, it defines a
    direct dependency of ``[category, name]``."""

    def __init__(self, species, options):
        self.species = species
        self.options = options

    @abstractmethod
    def define(self, res: dict, par:dict):
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

    @property
    def parameter_structure(self) -> dict:
        """Virtual parameter to propose a parameter structure, given the
        configuration of the contribution, in particular ``species``.
        The default implementation is empty, such that contributions with
        no parameters need not to overwrite this method.

        By convention, all (leaf) values of this data structure are to be
        ``None``.
        """
        return {}

    def relax(self,  # pylint: disable=R0201
              current_result: dict,
              delta_state: Collection[float]) -> float:
        """Virtual function to report the maximal allowable step size in
        the state variables.

        :current_result: The results based on the current state
        :delta_state: The given direction

        :seealso: :meth:`ThermoFrame.relax`
        """
        # default implementation
        del current_result, delta_state  # unused
        return 100  # a number greater than 1 / gamma for practical gamma

    def initial_state(self,  # pylint: disable=R0201
                      T: float,
                      p: float,
                      n: Collection[float],
                      parameters: dict) -> List[float]:
        """When the :class:`ThermoFrame object is queried for an initial state
        representation, The upper-most contribution that implements this
        method and does not return ``None`` takes the responsibility of
        calculating that state. If no contribution implements this method, the
        state is assumed to be in Gibbs coordinates, trivially ``T``, ``p``,
        and ``n``.

        Hence, normally only Helmholtz models need to implement this method.
        The true model coordinates can however be entirely unconventionally,
        such that it is solely up to the contributions on how to obtain the
        initial state.

        :param T: Temperature [K]
        :param p: Pressure [Pa]
        :param n: Molar quantity [mol]
        :param parameters: The parameter structure, equivalent as provided to
          :meth:`relax`.
        :return: The initial state

        .. seealso:: :meth:`ThermoFrame.initial_state`
        """
        del T, p, n, parameters  # unused

    @staticmethod
    def _tensor_structure(*keys) -> dict:
        """A helper method to create a nested rectabgular structure as often
        required for implementing :meth:`parameter_structure`.

        For instance:

        .. code-block::

            >>> struct = ThermoContribution._tensor_structure
            >>> print(struct(["A", "B"], ["C", "D", "E"])

        yields ::

            {
              "A": {"C": None, "D": None, "E": None},
              "B": {"C": None, "D": None, "E": None}
            }

        :param keys: Each element of ``keys`` is a list of strings, describing
          the keys to define at each level.

        :return: The nested dictionary with ``None`` as leaf values. Note that
          the equal inner dictionaries are individual copies. This to prevent
          a potential bug, if the client code just replaces the ``None`` values
          with actual data.

        """
        if not keys:
            return None
        current, *rest = keys
        func = ThermoContribution._tensor_structure
        return {k: func(*rest) for k in current}

    def _vector(self, dictionary: dict, keys: Collection[str] = None) -> SX:
        """During :meth:`define`, sub-directories often require to be converted
        into a single ``casadi.SX`` vector. This method provides such
        functionality.

        .. code-block::

            >>> vector = ThermoContribution._vector
            >>> struct = {"A": SX.sym("x"),
                          "C": SX.sym("y"),
                          "B": SX.sym("z")}
            >>> print(vector(struct, ["A", "B", "C"])

        yields ::

            SX([x, z, y])

        :param dictionary: The parameter structure, pointing to the node where
          to extract the vector from
        :param keys: The keys of the sub-dictionary to collect the ``casadi``
          symbols from. If left as ``None``, the species names of the
          contribution are assumed.
        :return: The ``casadi.SX`` object containing the collected symbols
        """
        if keys is None:
            keys = self.species
        return vertcat(*[dictionary[key] for key in keys])
