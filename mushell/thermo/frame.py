# -*- coding: utf-8 -*-

# stdlib modules
from copy import deepcopy

# external modules
from casadi import Function, SX

# internal modules
from ..utilities import flatten_dictionary, unflatten_dictionary

# TODO:
#  Overwrite __call__ and take parameter struct as input for parameters.


class ThermoFactory:
    def __init__(self):
        self.__contributions = {}

    def register_contribution(self, name, contribution_cls):
        self.__contributions[name] = contribution_cls

    def create_frame(self, configuration):
        species = configuration["species"]
        options = configuration.get("options", {})
        contributions = {name: [self.__contributions[name],
                                options.get(name, {})]
                         for name in configuration["contributions"]}
        property_names = configuration.get("properties", None)
        return ThermoFrame(species, contributions, property_names)


class ThermoFrame:
    def __init__(self, species, contributions, property_names=None):
        """This constructor establishes a thermo frame function object
        (casadi) with given species and contributions."""

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
    def function(self):
        return self.__function

    @property
    def species(self):
        return deepcopy(self.__species)

    @property
    def property_names(self):
        return self.__function.name_out()

    @parameters.setter
    def parameters(self, value):
        self.__parameter_structure = value

    @property
    def parameters(self):
        """This property is to aid the process of parametrising a model.
        It returns the structure of all required model parameters as a nested
        dictionary. Initially, all parameters are set to ``None``.
        This property must be set with actual values before the object can be
        called, or the methods ``relax`` and ``initialise`` be invoked.
        """
        return self.__parameter_structure

    @property
    def parameter_names(self):
        """Returns the names of parameters to be set. The parameters are
        dot-separated paths, containing contribution name, parameter name(s),
        and species name(s), e.g. ``Showmate heat capacity.Cp.A.H2O``.
        On the last level, all parameters are scalar and floating point.
        Discrete parameters (e.g. integers and strings) must be defined as
        configuration parameters.
        """
        # the parameter names can be queried from the contributions. This
        # yields a hierarchy (use None as values).
        # This property is the flattened variant using
        # utilities.flatten_dictionary
        return list(self.__parameter_names)

    def relax(self, state, delta_state, parameters):
        return min([cont.relax(state, delta_state, parameters.get("name", {}))
                    for name, cont in self.__contributions.items()])

    def initial_state(self, T, p, n, parameters):
        """Return a state estimate for given temperature, pressure and
        molar quantities - at given parameter set.

        This method queries all contributions top-down for an implementation
        of the initialise method. If not overwritten by any (and thus returning
        a state), the model is expected to be in Gibbs coordinates already,
        and (T, p, n) is the initial state.
        """
        for name, cont in reversed(self.__contributions.items()):
            result = cont.initial_state(T, p, n, parameters.get(name, {}))
            if result:
                return result
        return [T, p] + list(n)
