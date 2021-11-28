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
from .contribution import ThermoContribution


class ThermoFrame:
    """This class represents the thermodynamic model, which defines a casadi
    Function object :attr:`function` of the state vector and the parameter
    vector, and calculates a set of thermodynamic properties.
    """
    def __init__(self,
                 species: List[str],
                 contributions: dict,
                 property_names: List[str] = None):
        """This constructor establishes a thermo frame function object
        (casadi) with given species and contributions.

        .. note::

            This constructor shall not be called directly, but is ivoked by
            the :meth:`ThermoFactory.create_frame` method that handles the
            house-keeping of contribution classes. For this reason, it is not
            further documented here.
        """

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
        for name, contribution in contributions.items():
            contribution.define(result, param_symbols.get(name, {}))

        # extract properties of interest
        if property_names is None:
            property_names = list(result.keys())
        properties = [result[name] for name in property_names]

        self.__function = Function("thermo_frame",
                                   [state, parameters], properties,
                                   ["state", "parameters"], property_names)

        self.__contributions = contributions

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
        return self.function(state, flatten_dictionary(self.parameters))


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
        called, or the methods :meth:`relax` and :meth:`initialise` be invoked.
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

    def relax(self,
              state: Collection[float],
              delta_state: Collection[float]) -> float:
        """As a thermodynamic function, this object's contributions hold
        information on the domain limits of the state variables. This is mostly
        as trivial as demanding positive temperature, pressure, or quantities,
        as these variables occur for instance as logarithm arguments. In
        some cases, elaborate constraints apply, such as the minimum allowed
        volume for an equation of state.

        This method, after querying the contributions, returns a floating
        point number representing the maximal accepted step size.

        :param state: The current state (starting point)
        :param delta_state: The direction vector (typically suggested by a
          solver)
        :return: The maximal allowed step size as a pre-factor to
          ``delta_state``. Hence a value of one describes full step length.
        """
        parameters = self.__parameter_structure
        return min([cont.relax(state, delta_state, parameters.get(name, {}))
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
        parameters = self.__parameter_structure
        for name, cont in reversed(self.__contributions.items()):
            result = cont.initial_state(temperature, pressure, quantities,
                                        parameters.get(name, {}))
            if result:
                return result
        return [temperature, pressure] + list(quantities)


class ThermoFactory:
    """The ``ThermoFactory`` class hosts the definitions for the *model
    contributions*, enabling it to create instances of thermodynamic models of
    class :class:`ThermoFrame`.

    The class is largely meant to be a singelton, but to keep doors open,
    multiple instances can be created."""
    def __init__(self):
        """Parameter-less constructor, initialising the data structure
        to host contribution definitions"""
        self.__contributions = {}

    def register_contribution(self, name: str,
                              contribution_cls: Type[ThermoContribution]):
        """Registers a contribution under a given name. The contribution must
        be a concrete subclass of :class:`ThermoContribution`.

        :param name: The name identifier of the contribution. A ``ValueError``
          is thrown if the identifier was already defined before.
        :param contribution_cls: The contribution class to register
        """
        if name in self.__contributions:
            raise ValueError(f"Contribution of name '{name}' already defined.")
        self.__contributions[name] = contribution_cls

    def create_frame(self, configuration: dict) -> ThermoFrame:
        """This factory method creates a :class:`ThermoFrame` object from the
        given ``configuration``.

        :param configuration: A nested dictionary with the
          following root entries:

            - ``species``: A list of strings, representing the names of the
              chemical species. These names are used as identifiers in later
              use of the model, and they define the names of thermodynamic
              parameters required by the model.
            - ``contributions``: A list of strings, representing the names
              of the contributions to stack. These identifiers must have been
              defined upfront by calls to :meth:`register_contribution`.
            - ``options``: This entry can be omitted if not applicable, but
              some contributions might support options. For addressing these,
              the sub-dictionary contains a key for such contribution, and the
              value is whatever data structure the contribution accepts.
            - ``properties``: An optional entry being a list of strings that
              define the names of the properties to be calculated.

        :return: The thermodynamic model object
        """
        species = configuration["species"]
        options = configuration.get("options", {})
        contributions = {name: [self.__contributions[name],
                                options.get(name, {})]
                         for name in configuration["contributions"]}
        property_names = configuration.get("properties", None)
        return ThermoFrame(species, contributions, property_names)
