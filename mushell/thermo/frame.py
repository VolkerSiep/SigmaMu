# -*- coding: utf-8 -*-
"""
This module defines the :cls:`ThermoFrame` and the :cls:`ThermoFactory`
classes, of which the latter one is a factory for the former.

:cls:`ThermoFrame` objects represent thermodynamic models as state function
objects, calculating thermophysical properties as function of their
state and the model parameters.
"""
# stdlib modules
from copy import deepcopy
from typing import Type, List, Collection

# external modules
from casadi import Function, SX

# internal modules
from ..utilities import flatten_dictionary, unflatten_dictionary
from .contribution import ThermoContribution, StateDefinition


class ThermoFrame:
    """This class represents the thermodynamic model, which defines a casadi
    Function object :attr:`function` of the state vector and the parameter
    vector, and calculates a set of thermodynamic properties.
    """
    def __init__(self,
                 species: List[str],
                 state_definition: StateDefinition,
                 contributions: dict):
        """This constructor establishes a thermo frame function object
        (casadi) with given species and contributions.

        .. note::

            This constructor shall not be called directly, but is ivoked by
            the :meth:`ThermoFactory.create_frame` method that handles the
            house-keeping of contribution classes. For this reason, it is not
            further documented here.
        """
        ThermoFrame.__check_dependencies(contributions.values())

        # need to instantiate the contributions
        contributions = {name: cls_(species, options)
                         for name, (cls_, options) in contributions.items()}

        self.__parameter_structure = {name: con.parameter_structure
                                      for name, con in contributions.items()}
        flat_params = flatten_dictionary(self.__parameter_structure)
        self.__parameter_names = flat_params.keys()

        # flat parameter vector for casadi function
        parameters = SX.sym('parameters', len(flat_params))

        # recreate dictionary for the contributions to read
        param_symbols = {name: parameters[k]
                         for k, name in enumerate(self.__parameter_names)}
        param_symbols = unflatten_dictionary(param_symbols)

        # define thermodynamic state (th, mc, [ch])
        state = SX.sym("x", len(species) + 2)
        self.__species = species

        # call the contributions
        result = {"state": state}
        state_definition.prepare(result)
        for name, contribution in contributions.items():
            contribution.define(result, param_symbols.get(name, {}))

        # extract properties of interest
        property_names = list(result.keys())
        properties = list(result.values())

        self.__function = Function("thermo_frame",
                                   [state, parameters], properties,
                                   ["state", "parameters"], property_names)

        self.__contributions = contributions
        self.__state_definition = state_definition

    @property
    def function(self) -> Function:
        """The ``casadi`` function object representing the thermodynamic model,
        having two arguments: ``state`` and ``parameters``.

        The ``state`` represents the independent variables of the state
        function. These are initially obtained from the :meth:`initial_state`
        method and then normally updated by some numerical scheme to solve the
        problem at hand.

        The ``parameters`` are the set of thermodynamic parameters as required
        by the model contributions. The :meth:`__call__` calls this function
        with ``parameters`` defined as ``flatten_dictionary(self.parameters)``.

        Instead of calling the call-operator, directly accessing the
        ``casadi.Function`` object allows to query further derivatives of the
        function with respect to the state and parameters.
        """
        return self.__function

    def __call__(self, state: Collection[float]) -> List[Collection[float]]:
        """Call to the function object :attr:`function`, in particular with
        current parameter set :attr:`parameters` and using evaluation with
        float typed variables.

        :param state: A collection of floats representing the state of the
          model, initially obtained by :meth:`initial_state`.

        :return: A list of property collections, representing the thermodynamic
          properties, in the sequence as defined by :attr:`property_names`.
        """
        parameters = flatten_dictionary(self.parameters).values()
        return self.function(state, parameters)

    @property
    def species(self) -> List[str]:
        """Returns a list of species names"""
        return deepcopy(self.__species)

    @property
    def property_names(self) -> List[str]:
        """Returns a list of property names, defining the content and sequence
        of returned properties from :meth:`__call__` and the function object
        :attr:`function`."""
        return self.__function.name_out()

    @property
    def parameters(self) -> dict:
        """This property is to aid the process of parametrising a model.
        It returns the structure of all required model parameters as a nested
        dictionary. Initially, all parameters are set to ``None``.
        This property must be set with actual values before the object can be
        called or :meth:`initialise` invoked.
        """
        return self.__parameter_structure

    @parameters.setter
    def parameters(self, value: dict):
        self.__parameter_structure = value

    @property
    def parameter_names(self) -> List[float]:
        """Returns the names of parameters to be set. The parameters are
        dot-separated paths, containing contribution name, parameter name(s),
        and species name(s), e.g. ``Showmate heat capacity.Cp.A.H2O``.
        On the last level, all parameters are scalar and floating point.
        Discrete parameters (e.g. integers and strings) must be defined as
        configuration parameters.
        """
        return list(self.__parameter_names)

    def relax(self, current_result: List[Collection[float]],
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
          ``delta_state``. Hence a value of one describes full step length.
        """
        result = dict(zip(self.property_names, current_result))
        return min([cont.relax(result, delta_state)
                    for name, cont in self.__contributions.items()])

    def initial_state(self,
                      temperature: float,
                      pressure: float,
                      quantities: Collection[float]) -> List[float]:
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
        state, non_state = st_def.reverse(temperature, pressure, quantities)
        if None not in state:
            return state

        state = [float("NaN") if x is None else x for x in state]
        # calculate all propeties ... accept NaNs
        properties = zip(self.property_names, self(state))
        properties = {name: value for name, value in properties}

        for name, cont in reversed(self.__contributions.items()):
            result = cont.initial_state(temperature, pressure, quantities,
                                        properties)
            if result:
                return result
        msg = "No initialisation found despite of non-Gibbs surface"
        raise NotImplementedError(msg)

    @staticmethod
    def __check_dependencies(contributions):
        """check compatibility and dependency of contributions"""
        cat_names, full_names = [], []
        for cls_, _ in contributions:
            for item in cls_.requires:
                msg = f"Dependency '{item}' not fulfilled in " \
                      f"contribution '{cls_.name}'"
                assert item in cat_names \
                    if type(item) is str else full_names, msg
            cat_names.append(cls_.category)
            full_names.append([cls_.category, cls_.name])


class ThermoFactory:
    """The ``ThermoFactory`` class hosts the definitions for the *model
    contributions*, enabling it to create instances of thermodynamic models of
    class :class:`ThermoFrame`.

    The class is largely meant to be a singelton, but to keep doors open,
    multiple instances can be created."""

    CATEGORY_SEPERATOR = "#"
    """The separator character to join category and name. It can be changed
    before any contributions are registered - otherwise yields inconisistent
    behaviour. The separator cannot be a dot ``.``, as this would conflict
    with the structuring of thermodynamic parameters."""

    def __init__(self):
        """Parameter-less constructor, initialising the data structure
        to host contribution definitions"""
        self.__contributions = {}
        self.__state_definitions = {}

    def register_state_definition(self, definition_cls: Type[StateDefinition]):
        name = definition_cls.name
        if name in self.__state_definitions:
            raise ValueError(f"State definition '{name}' already defined.")
        self.__state_definitions[name] = definition_cls

    def register_contribution(self,
                              contribution_cls: Type[ThermoContribution]):
        """Registers a contribution under the name that concatenates
        ``contribution_cls.category`` and ``contribution_cls.name``
        attribute. The contribution must be a concrete subclass of
        :class:`ThermoContribution`.
        :param contribution_cls: The contribution class to register
        """
        name = contribution_cls.category
        if contribution_cls.name:
            name += ThermoFactory.CATEGORY_SEPERATOR + contribution_cls.name
        if name in self.__contributions:
            raise ValueError(f"Contribution '{name}' already defined.")
        self.__contributions[name] = contribution_cls

    @property
    def contribution_names(self) -> List[str]:
        """This property contains the full names of all so long registered
        properties"""
        return sorted(self.__contributions.keys())


    def create_frame(self, configuration: dict) -> ThermoFrame:
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
            - ``options``: This entry can be omitted if not applicable, but
              some contributions might support options. For addressing these,
              the sub-dictionary contains a key for such contribution, and the
              value is whatever data structure the contribution accepts.

        :return: The thermodynamic model object
        """
        species = configuration["species"]
        options = configuration.get("options", {})
        state_def_cls = self.__state_definitions[configuration["state"]]
        state_definition = state_def_cls(species)
        contributions = {name: [self.__contributions[name],
                                options.get(name, {})]
                         for name in configuration["contributions"]}
        return ThermoFrame(species, state_definition, contributions)
