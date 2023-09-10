# -*- coding: utf-8 -*-
"""
This module defines the :cls:`ThermoFrame` and the :cls:`ThermoFactory`
classes, of which the latter one is a factory for the former.

:cls:`ThermoFrame` objects represent thermodynamic models as state function
objects, calculating thermochemical properties as function of their
state and the model parameters.
"""
# stdlib modules
from typing import Type
from collections.abc import Mapping, Collection
from logging import getLogger

# internal modules
from .contribution import ThermoContribution, StateDefinition
from ..utilities import (Quantity, ParameterDictionary, QFunction,
                         SymbolQuantity, extract_units_dictionary)
from ..utilities.types import NestedMap, Map

ThermoContributionDict = Map[tuple[Type["Contribution"], Map]]

logger = getLogger(__name__)


class ThermoFrame:
    """This class represents the thermodynamic model, which defines a casadi
    Function object :attr:`function` of the state vector and the parameter
    vector, and calculates a set of thermodynamic properties.
    """

    def __init__(self, species: Collection[str],
                 state_definition: StateDefinition,
                 contributions: ThermoContributionDict):
        """This constructor establishes a thermo frame function object
        (casadi) with given species and contributions.

        .. note::

            This constructor shall not be called directly, but is invoked by
            the :meth:`ThermoFactory.create_frame` method that handles the
            house-keeping of contribution classes. For this reason, it is not
            further documented here.
        """
        # need to instantiate the contributions
        contributions = {
            name: cls_(species, options)
            for name, (cls_, options) in contributions.items()
        }

        parameters = {}

        # define thermodynamic state (th, mc, [ch])
        state = SymbolQuantity("x", "dimless", len(species) + 2)
        self.__species = list(species)

        # call the contributions
        result = {"_state": state}
        state_definition.prepare(result)
        for name, contribution in contributions.items():
            params = ParameterDictionary()
            contribution.define(result, params)
            logger.debug(f"Defining contribution '{name}'")
            if params:
                parameters[name] = params

        args = {"state": state, "parameters": parameters}
        self.__function = QFunction(args, result, "thermo_frame")

        self.__contributions = contributions
        self.__state_definition = state_definition
        self.__default = None
        self.__param_struct = extract_units_dictionary(parameters)

    def __call__(self, state: Quantity, parameters: NestedMap[Quantity],
                 squeeze_results: bool = True):
        """Shortcut: Call to the function object :attr:`function`.

        :return: A list of property collections, representing the thermodynamic
          properties, in the sequence as defined by :attr:`property_names`.
        """
        return self.__function({
                "state": state,
                "parameters": parameters
            }, squeeze_results)

    @property
    def species(self) -> Collection[str]:
        """Returns a list of species names"""
        return list(self.__species)

    @property
    def property_structure(self) -> NestedMap[str]:
        """Returns a recursive structure properties, defining the calculated
        properties from :meth:`__call__` """
        return self.__function.result_structure

    @property
    def parameter_structure(self) -> NestedMap[str]:
        """This property is to aid the process of parametrizing a model.
        It returns the structure of all required model parameters. Initially,
        the returned object contains units of measurements that must be
        replaced with actual quantities (symbolic or not) before the
        function can be called or :meth:`initialise` invoked. For the latter,
        float quantities have to be provided to the parameter object.
        """
        return self.__param_struct

    def relax(self, current_result: NestedMap[Quantity],
              delta_state: Collection[float]) -> float:
        """As a thermodynamic function, this object's contributions hold
        information on the domain limits of the state variables. This is mostly
        as trivial as demanding positive temperature, pressure, or quantities,
        as these variables occur for instance as logarithm arguments. In
        some cases, elaborate constraints apply, such as the minimum allowed
        volume for an equation of state.

        This method, after querying the contributions, returns a floating
        point number representing the maximal accepted step size.

        :param current_result: The result objects from the calculation with
          the current state
        :param delta_state: The direction vector (typically suggested by a
          solver)
        :return: The maximal allowed step size as a pre-factor to
          ``delta_state``. Hence, a value of one describes full step length.
        """
        return float(
            min(
                cont.relax(current_result, delta_state)
                for cont in self.__contributions.values()))

    def initial_state(self, temperature: Quantity, pressure: Quantity,
                      quantities: Quantity, parameters) -> list[float]:
        """Return a state estimate for given temperature, pressure and
        molar quantities - at given parameter set.

        This method queries all contributions top-down for an implementation
        of the initialise method. If not overwritten by any (and thus returning
        a state), the model is expected to be in Gibbs coordinates already,
        and (T, p, n) is the initial state.

        """
        # call own function (assert that it is defined yet)
        #   defining the state as [T, NaN] + n.

        # then start from top and call contributions to try to initialise,
        #  based on results that might contain NaNs (if they depend on volume).
        st_def = self.__state_definition
        state = st_def.reverse(temperature, pressure, quantities)
        if None not in state:
            return state

        state = Quantity([float("NaN") if x is None else x for x in state])
        # calculate all properties ... accept NaNs
        properties = self(state, parameters)

        for cont in reversed(self.__contributions.values()):
            result = cont.initial_state(temperature, pressure, quantities,
                                        properties)
            if result:
                result[0] = float(result[0])
                result[1] = float(result[1])
                return result
        msg = "No initialisation found despite of non-Gibbs surface"
        raise NotImplementedError(msg)

    # TODO: Should default state be better in terms of quantities? This
    #  because its in the interpretation of T, p, n.

    @property
    def default(self):
        """The definition of the object can optionally contain a default state.
        If this is applied, the given default state is stored in this
        property. Its interpretation is always in ``T, p, n`` coordinates."""
        return self.__default

    @default.setter
    def default(self, state):
        if state is not None:
            num_species = len(self.species)
            if len(state) != 3:
                raise ValueError(
                    "Default state must contain three elements, " +
                    f"found {len(state)} instead.")
            if len(state[2]) != num_species:
                raise ValueError(f"Default state must cover {num_species} " +
                                 f"species, found {len(state[2])} instead.")
        self.__default = state


class ThermoFactory:
    """The ``ThermoFactory`` class hosts the definitions for the *model
    contributions*, enabling it to create instances of thermodynamic models of
    class :class:`ThermoFrame`.

    The class is largely meant to be a singleton, but to keep doors open,
    static attributes are avoided."""

    def __init__(self):
        """Parameter-less constructor, initialising the data structure
        to host contribution definitions"""
        self.__contributions = {}
        self.__state_definitions = {}

    def register_state_definition(self, definition_cls: Type[StateDefinition]):
        """Register a new state definition with the name of its class.

        :param definition_cls: The state definition to register
        """
        name = definition_cls.__name__
        if name in self.__state_definitions:
            raise ValueError(f"State definition '{name}' already defined.")
        self.__state_definitions[name] = definition_cls

    def register(self, *contributions: Type[ThermoContribution]):
        """Registers contributions under the name of the class.
        The contributions must be a concrete subclass of
        :class:`ThermoContribution`.

        :param contributions: The contributions to register, being classes
          (not instances)
        """
        for class_ in contributions:
            name = class_.__name__
            if name in self.__contributions:
                raise ValueError(f"Contribution '{name}' already defined.")
            self.__contributions[name] = class_

    @property
    def contribution_names(self) -> Collection[str]:
        """This property contains the full names of all so long registered
        properties"""
        return sorted(self.__contributions.keys())

    def create_frame(self, configuration: Mapping) -> ThermoFrame:
        """This factory method creates a :class:`ThermoFrame` object from the
        given ``configuration``.

        :param configuration: A nested dictionary with the
          following root entries:

            - ``species``: A list of strings, representing the names of the
              chemical species. These names are used as identifiers in later
              use of the model, and they define the names of thermodynamic
              parameters required by the model.
            - ``state``: An identifier that represents the type of state as
              defined by :meth:`register_state_definition`.
            - ``contributions``: A list of strings, representing the names
              of the contributions to stack. These identifiers must have been
              defined upfront by calls to :meth:`register_contribution`.
              A contribution entry can also be a dictionary with the following
              keys:

                - ``cls``: The contribution class (required). If this is the
                  only key defined, its effect is as if the sole string of the
                  class name was provided.
                - ``name``: The name of the contribution in the created
                  :class:`ThermoFrame`` object. This name must be unique within
                  the frame definition. It will also be used to define the
                  contribution parameter structure.
                  If skipped, the name will be the same as ``cls``. To be
                  unique, this doesn't work if one contribution class is used
                  multiple times.
                - ``options``: Any data structure that is accepted (and
                  hopefully documented) for the particular contribution.
                  If skipped, an empty dictionary is used.

        :return: The thermodynamic model object
        """
        contributions = {}
        for item in configuration["contributions"]:
            if isinstance(item, dict):
                class_ = item["cls"]
                name = item.get("name", class_)
                options = item.get("options", {})
            else:
                name, class_, options = item, item, {}
            class_ = self.__contributions[class_]
            if name in contributions:
                raise ValueError(f"Duplicate contribution name '{name}'")
            contributions[name] = class_, options

        species = configuration["species"]
        state_def_cls = self.__state_definitions[configuration["state"]]
        result = ThermoFrame(species, state_def_cls(), contributions)

        # set default state
        default = configuration.get("default_state", None)
        if default is not None:
            # make sure the values are float
            default = [
                float(default[0]),
                float(default[1]),
                list(map(float, default[2]))
            ]
        result.default = default
        return result
