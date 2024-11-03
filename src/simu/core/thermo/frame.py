# -*- coding: utf-8 -*-
"""
This module defines the :cls:`ThermoFrame` and the :cls:`ThermoFactory`
classes, of which the latter one is a factory for the former.

:cls:`ThermoFrame` objects represent thermodynamic models as state function
objects, calculating thermochemical properties as function of their
state and the model parameters.
"""
# stdlib modules
from typing import Type, Optional
from collections.abc import Mapping, Sequence
from logging import getLogger

# external modules
from casadi import SX

# internal modules
from .contribution import ThermoContribution
from .state import StateDefinition, InitialState
from .species import SpeciesDefinition
from ..utilities import (Quantity, ParameterDictionary, QFunction,
                         SymbolQuantity, extract_units_dictionary)
from ..utilities.types import NestedMap, Map, MutMap

ThermoContributionDict = Map[tuple[Type["ThermoContribution"], Map]]
"""
A dictionary whose keys are the names of the contributions, and the values are
tuples of the belonging classes and the options belonging to the definition.
An example could be:

.. code-block::

    contribs = {
        "LinearHeatCapacity": (LinearHeatCapacity, {}),
        ...
        "MixingRule_A": (NonSymmetricMixingRule, {"target": "_ceos_a"}
    } 
    
This kind of structure is used to define the sequence of contributions in a 
:cls:`ThermoFrame` object. Here, the class ``LinearHeatCapacity`` is defined 
straight forward with no options. For the mixing rule of the ``A`` contribution,
called ``MixingRule_A``, it uses the ``NonSymmetricMixingRule`` class, 
configured with ``ceos_a`` being the ``target``.

It is up to each individual :cls:`ThermoContribution` implementation to support
and document their set of options.
"""

logger = getLogger(__name__)


class ThermoFrame:
    """This class represents the thermodynamic model, which defines a casadi
    Function object :attr:`function` of the state vector and the parameter
    vector, and calculates a set of thermodynamic properties.

    The object should not be constructed via the class constructor by
    client code, but created by the :meth:`ThermoFactory.create_frame`
    method that handles the house-keeping of contribution classes.
    """

    def __init__(self, species: Mapping[str, SpeciesDefinition],
                 state_definition: StateDefinition,
                 contributions: ThermoContributionDict):
        """This constructor establishes a thermo frame function object
        (casadi) with given species and contributions.
        """
        # need to instantiate the contributions
        species_list = list(species.keys())
        contribs: Map[ThermoContribution] = {
            name: cls_(species_list, options)
            for name, (cls_, options) in contributions.items()
        }
        self.__species: Mapping[str, SpeciesDefinition] = species

        def create_function(flow: bool = False):
            """Build up the result dictionary and define function"""
            params: MutMap[ParameterDictionary] = {}
            # define thermodynamic state (th, mc, [ch])
            state = SymbolQuantity("x", "dimless", len(species) + 2)
            # call the contributions; build up result dictionary
            result = {"_state": state}
            state_definition.prepare(result, flow)
            for name, contribution in contribs.items():
                new_params = ParameterDictionary()
                contribution.define(result, new_params)
                logger.debug(f"Defining contribution '{name}'")
                if new_params:
                    params[name] = new_params
            # create function
            args = {"state": state, "parameters": params}
            return QFunction(args, result, "thermo_frame"), params

        self.__function, parameters = create_function(flow=False)
        function, _ = create_function(flow=True)

        self.__res_units = {
            False: self.__function.res_units,
            True: function.res_units
        }

        self.__contributions: Map[ThermoContribution] = contribs
        self.__state_definition: StateDefinition = state_definition
        self.__default: Optional[InitialState] = None
        self.__param_struct: NestedMap[str] = \
            extract_units_dictionary(parameters)

    def __call__(self, state: SX | Sequence[float],
                 parameters: NestedMap[Quantity],
                 squeeze_results: bool = True, flow: bool = False):
        """Shortcut operator to call to the underlying function object.

        The function call can be a stand-alone evaluation of the thermodynamic
        model, given floating point quantities for the state and the
        parameters. Alternatively, called with ``CasADi`` ``SX`` based
        quantities to become part of a larger functional.

        :param state: A ``CasADi`` ``SX`` object or a sequence of floats,
          representing the thermodynamic state of the model. This is to be seen
          as a purely numerical object, as the physical interpretation, for
          instance as temperature, volume and mole flows, is first happening
          within the model.

        :param parameters: A nested dictionary with string keys and
          :class:`simu.Quantity` leaves. Depending on the application, these
          quantities hold float or ``SX`` type magnitudes.

        :return: A list of property collections, representing the thermodynamic
          properties, in the sequence as defined by :attr:`property_names`.

        Depending on the ``flow`` parameter, extensive properties will be
        returned as flow quantities, such as ``W`` or ``mol/s`` or as stagnant
        state properties, such as ``J`` or ``mol``.
        """
        self.__function.res_units = self.__res_units[flow]
        return self.__function({
                "state": Quantity(state),
                "parameters": parameters
            }, squeeze_results)

    @property
    def species(self) -> Sequence[str]:
        """Returns a list of species names"""
        return list(self.__species.keys())

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
        function (:meth:`__call__`) or :meth:`initial_state` can be called .
        For the latter, float quantities have to be provided to the parameter
        object.
        """
        return self.__param_struct

    def create_symbol_state(self) -> SX:
        """Create a symbol state that can be used to call the object functional
        symbolically."""
        return SX.sym("x", len(self.__species) + 2)

    def relax(self, current_result: Map[Quantity],
              delta_state: Sequence[float]) -> float:
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

    def initial_state(self, state: InitialState,
                      parameters: NestedMap[Quantity]) -> Sequence[float]:
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
        state_flat: Sequence[float] = st_def.reverse(state)
        if None not in state_flat:
            return state_flat

        state_flat = [float("NaN") if x is None else x for x in state_flat]
        # calculate all properties ... accept NaNs
        properties = self(state_flat, parameters)

        for cont in reversed(self.__contributions.values()):
            result = cont.initial_state(state, properties)
            if result:
                result[0] = float(result[0])
                result[1] = float(result[1])
                return result
        msg = "No initialisation found for non-Gibbs surface"
        raise NotImplementedError(msg)
