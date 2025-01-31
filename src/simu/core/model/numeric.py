"""This module implements functionality concerning the numerical handling
of the top model instance."""

from typing import Optional
from collections.abc import Callable, Sequence, Collection
from copy import deepcopy
from casadi import vertcat, SX
from pint import Unit

from ..utilities import (
    flatten_dictionary, quantity_dict_to_strings,
    parse_quantities_in_struct, nested_map)
from ..utilities.types import NestedMap, NestedMutMap, Map, MutMap
from ..utilities.quantity import Quantity, QFunction, jacobian
from ..utilities.structures import FLATTEN_SEPARATOR
from ..utilities.errors import DataFlowError

from .base import ModelProxy
from ..thermo import ThermoParameterStore
from ... import unflatten_dictionary, InitialState


# TODO:
#  - set parameters and get parameters


class NumericHandler:
    """This class implements the function object describing the top level
    model."""
    THERMO_PARAMS: str = "thermo_params"
    MODEL_PARAMS: str = "model_params"
    MODEL_PROPS: str = "model_props"
    THERMO_PROPS: str = "thermo_props"
    RESIDUALS: str = "residuals"
    STATE_VEC: str = "states"
    RES_VEC: str = "residuals"
    BOUND_VEC: str = "bounds"
    VECTORS: str = "vectors"
    JACOBIANS: str = "jacobians"

    def __init__(self, model: ModelProxy, port_properties: bool = True):
        """The option ``port_properties`` determines whether the properties
        of connected materials are also reported from a child model's
        perspective by the name of their ports."""
        self.options = {
            "port_properties": port_properties
        }
        self.model = model
        # the name vectors of vector arguments
        self.__vec_arg_names: MutMap[Sequence[str]] = {}
        self.__vec_res_names: MutMap[Sequence[str]] = {}

        # the symbolic argument structure
        self.__symargs: NestedMutMap[Quantity] = self.__collect_arguments()
        # the symbolic result structure
        self.__symres: NestedMutMap[Quantity] = self.__collect_results()
        # the numerical argument structure with initial values
        self.__arguments: MutMap[Quantity] = {}

    @property
    def function(self) -> QFunction:
        """Create a Function object based on currently available argument
        and result structures."""
        return QFunction(self.__symargs, self.__symres, "model")

    def vector_arg_names(self, key: str) -> Sequence[str]:
        """Return the names for the argument vector of given ``key``"""
        return self.__vec_arg_names[key]

    def vector_res_names(self, key: str) -> Sequence[str]:
        """Return the names for the result vector of given ``key``"""
        return self.__vec_res_names[key]

    @property
    def arguments(self) -> Map[Quantity]:
        """The function arguments as numerical values. A DataFlowError is
        thrown, if not all numerical values are known.
        A deepcopy of the structure is provided, so the returned data can be
        altered without side effects.
        """
        if not self.__arguments:
            self.__arguments = self.__collect_argument_values()
        return deepcopy(self.__arguments)

    def export_state(self) -> NestedMutMap[str]:
        """Export the internal state of the model in a hierarchical structure,
        whereas all thermodynamic states are given in :math:`T, p, n`.

        As by the philosophy of the chosen approach, only the thermodynamic
        models know how to obtain their internal state from any :math:`T, p, n`
        specification.

        The returned structure is meant to be easy to store for instance in
        yaml or json format, and easy to edit. One can use
        :func:`~simu.parse_quantities_in_struct` to convert the values of the
        data structure into :class:`~simu.Quantity` objects for programmatic
        processing.

        """
        def fetch_initial_states(model: ModelProxy) -> MutMap[Quantity]:
            """fetch material states from a specific model"""
            mat_proxy = model.materials
            return {k: m.initial_state.to_dict(m.species)
                    for k, m in mat_proxy.handler.items()
                    if k not in mat_proxy}

        thermo =  self.__fetch(self.model, fetch_initial_states, "state")
        # TODO: when non-canonical states are implemented, collect them here.

        return quantity_dict_to_strings(
            {"thermo": thermo,
             "non-canonical": {}}
        )

    def import_state(self, state: NestedMap[str],
                     allow_missing: bool=False, allow_extra: bool = False)\
            -> NestedMap[str]:
        """Imports the state data in terms of :math:`T, p, n` as exported by
        :meth:`export_state`.

        :param state: A nested mapping as returned by :meth:`export_state`,
          containing the two first-level keys ``thermo`` (for thermodynamic
          states) and ``non-canonical`` for non-canonical states.
        :param allow_missing: If true, throw a ``ValueError`` if there are model
          states in the model that are not defined in the given ``state``
          structure.
        :param allow_extra: If false, throw a ``ValueError`` if there are states
          in the given ``state`` structure that are not present in the model.
        :return: A nested mapping of same structure as ``state``, but containing
          states that are not present in the model (value = ``extra``) and
          states that were not given as part of ``state`` (value = ``missing``)
        """
        def mk_new_path(path: str, name: str) -> str:
            name = name.replace(FLATTEN_SEPARATOR, rf"\{FLATTEN_SEPARATOR}")
            return name if not path else f"{path}{FLATTEN_SEPARATOR}{name}"

        def traverse(model: ModelProxy, state_part: NestedMap[Quantity],
                     path: str):
            # process local material objects
            all_names = set()
            for name, material in model.materials.handler.items():
                new_path = mk_new_path(path, name)
                all_names.add(name)
                try:
                    new_part = state_part[name]
                except KeyError:
                    if not allow_missing:
                        raise
                    result[new_path] = "missing"
                else:
                    material.initial_state = \
                        InitialState.from_dict(new_part, material.species)

            # traverse down into model hierarchy
            for name, proxy in model.hierarchy.handler.items():
                all_names.add(name)
                new_path = mk_new_path(path, name)
                try:
                    new_part = state_part[name]
                except KeyError:
                    if not allow_missing:
                        raise
                    result[new_path] = "missing"
                else:
                    traverse(proxy, new_part, new_path)

            # detect states that are not defined in model
            for name in state_part.keys():
                new_path = mk_new_path(path, name)
                if not name in all_names:
                    if not allow_extra:
                        msg = f"{new_path} not found in model structure"
                        raise KeyError(msg)
                    result[new_path] = "extra"

        self.__arguments = {}  # force reread
        result = {}
        traverse(self.model, parse_quantities_in_struct(state["thermo"]), "")
        return unflatten_dictionary(result)


    def retain_state(self, state: Sequence[float],
                     parameters: NestedMap[Quantity]):
        """Given a numeric ``state`` vector and the current set of
        ``parameters`` as a nested mapping of quantities, store the values for
        temperature, pressure and molar quantities back into the internal
        representations of the initial thermodynamic states.
        """
        def fetch_retain_initial_state(model: ModelProxy):
            """retain initial states for a specific model"""
            nonlocal index
            for m in model.materials.handler.values():
                state_length = 2 + len(m.species)
                state_part = state[index:index + state_length]
                m.retain_initial_state(state_part, parameters)
                index += state_length

        index = 0  # mutable index
        self.__arguments = {}  # force reread
        self.__traverse(self.model, fetch_retain_initial_state)

        # todo: if there are non-canonical states, treat them now.


    def extract_parameters(self, key: str,
                           definition: NestedMap[str]) -> Quantity:
        """collect the parameter symbols addressed by the ``definition``
        argument, which defines the unit of measurement for each parameter -
        as the numerical parameter vector elements must be considered
        dimensionless.

        .. note::

            Unfortunately, we cannot easily allow unit conversion (even
            compatible units) at this point, as the parameters are already
            independent variables used to build up the symbolic graph.
            Well, it can be done by first building the original function,
            and then call the function as f(c(x)), where c(x) is the unit
            conversion.

        These parameter symbols and values are then removed from the original
        argument structure and instead added to the vector entry as
        dimensionless entities.

        :param key: The name to be used for the parameter set
        :param definition: The parameters to be extracted
        """
        def traverse(parameters: NestedMap[str],
                     symbols: NestedMutMap[Quantity],
                     arguments: NestedMutMap[Quantity]) -> \
                (Sequence[str], Sequence[SX], Sequence[float]):
            """recursively go through struct, extract and remove both values
            and symbols."""
            try:
                items = parameters.items()
            except AttributeError:
                return None, None, None

            nams, syms, args = [], [], []
            for k, value in items:
                n, s, v = traverse(value, symbols[k], arguments[k])
                if s is None:
                    if symbols[k].units != Unit(value):
                        msg = "No unit conversion possible for parameter " \
                            f"{k}: from {symbols[k].units:~} to {value}."
                        raise ValueError(msg)
                    nams.append(k)
                    syms.append(symbols[k].magnitude)
                    args.append(arguments[k].magnitude)
                    del symbols[k]
                    del arguments[k]
                else:
                    nams.extend([f"{k}/{n_i}" for n_i in n])
                    syms.extend(s)
                    args.extend(v)
            return nams, syms, args

        if key in self.__symargs[NumericHandler.VECTORS]:
            msg = f"A parameter vector of name '{key}' is already used."
            raise KeyError(msg)

        if not self.__arguments:
            self.__arguments = self.__collect_argument_values()
        nam, sym, arg = traverse(definition, self.__symargs, self.__arguments)

        result = Quantity(vertcat(*sym))
        self.__symargs[NumericHandler.VECTORS][key] = result
        values = Quantity(arg)
        self.__arguments[NumericHandler.VECTORS][key] = values
        self.__vec_arg_names[key] = nam
        return result

    def collect_properties(self, key: str,
                           definition: NestedMap[str]) -> Quantity:
        """Collect the property symbols addressed by the ``definition``
        argument, which defines the unit of measurement for each property
        to be scaled with - as the numerical property vector elements must
        be dimensionless."""
        def traverse(properties: NestedMap[str],
                     symbols: NestedMap[Quantity]) -> \
                (Sequence[str], Sequence[SX]):
            """Recursively go through property struct, collect symbols, and
            convert them into desired units."""
            try:
                items = properties.items()
            except AttributeError:
                return None, None

            nams, syms = [], []
            for k, value in items:
                n, s = traverse(value, symbols[k])
                if n is None:
                    nams.append(k)
                    syms.append(symbols[k].to(value).magnitude)
                else:
                    nams.extend([f"{k}/{n_i}" for n_i in n])
                    syms.extend(s)
            return nams, syms

        if key in self.__symres[NumericHandler.VECTORS]:
            msg = f"A property vector of name '{key}' is already used."
            raise KeyError(msg)

        nam, sym = traverse(definition, self.__symres)

        result = Quantity(vertcat(*sym))
        self.__symres[NumericHandler.VECTORS][key] = result
        self.__vec_res_names[key] = nam
        return result

    def register_jacobian(self, dependent: str, independent: str) -> str:
        """Add the given symbols to the jacobian structure. These symbols must
        be a function of the arguments, or else the function cannot be
        created. The key must be unique.
        """
        dep = self.__symres[self.VECTORS][dependent]
        ind = self.__symargs[self.VECTORS][independent]
        jac = jacobian(dep, ind)
        key = f"d_({dependent})/d_({independent})"
        self.__symres[self.JACOBIANS][key] = jac
        return key

    def __collect_arguments(self) -> NestedMutMap[Quantity]:
        """Create a function that has the following arguments, each of them as
        a flat dictionary:

            - Material States
            - Model Parameters
            - Thermodynamic Parameters

        For child models, only the free parameters are collected.
        """
        mod = self.model
        fetch = self.__fetch
        to_vector = self.__to_vector

        def fetch_material_states(model: ModelProxy) -> MutMap[Quantity]:
            """fetch material states from a specific model"""
            mat_proxy = model.materials
            return {k: m.sym_state for k, m in mat_proxy.handler.items()
                    if k not in mat_proxy}

        def fetch_parameters(model: ModelProxy) -> MutMap[Quantity]:
            """fetch model parameters from a specific model"""
            return dict(model.parameters.free)

        def fetch_store_param(model: ModelProxy) -> NestedMap[Quantity]:
            """fetch thermodynamic parameters from the stores"""
            stores = self.__fetch_thermo_stores(model)
            names = {store.name for store in stores}
            if len(names) < len(stores):
                raise ValueError("When using multiple ThermoPropertyStores, "
                                 "they have to have unique names")
            return {store.name: store.get_all_symbols() for store in stores}

        states_struct = fetch(mod, fetch_material_states, "state")
        states, state_names = to_vector(states_struct)
        self.__vec_arg_names[self.STATE_VEC] = state_names

        return {
            self.THERMO_PARAMS: fetch_store_param(mod),
            self.MODEL_PARAMS: fetch(mod, fetch_parameters, "parameter"),
            self.VECTORS: {
                self.STATE_VEC: states,
            }
        }

    def __collect_results(self) -> NestedMutMap[Quantity]:
        """The result of the function consists of

            - Model Properties
            - Thermodynamic (state) properties
            - Residuals

        All the data is to be collected from the model and all child model
        proxies.
        """
        def fetch_residuals(model: ModelProxy,
                            normed: bool = False) -> NestedMutMap[Quantity]:
            """fetch residuals from a specific model"""
            def extract(entity):
                if normed:
                    return (entity.value / entity.tolerance).to("")
                return entity.value

            # find residuals of materials
            mat_proxy = model.materials
            res = {k: m.residuals(normed) for k, m in mat_proxy.handler.items()
                   if k not in mat_proxy}
            # add residuals of model (detect name clashes)
            clash = set(res.keys()) & set(model.residuals.keys())
            if clash:
                clash = ", ".join(clash)
                msg = f"Name clash of residuals and child modules: {clash}"
                raise ValueError(msg)

            res.update({k: extract(v) for k, v in model.residuals.items()})
            return res

        def fetch_bounds(model: ModelProxy) -> MutMap[Quantity]:
            mat_proxy = model.materials
            res = {k: m.bounds for k, m in mat_proxy.handler.items()
                   if k not in mat_proxy}
            clash = set(res.keys()) & set(model.bounds.keys())
            if clash:
                clash = ", ".join(clash)
                msg = f"Name clash of bounds and child modules: {clash}"
                raise ValueError(msg)
            res.update(model.bounds)
            return res

        def fetch_mod_props(model: ModelProxy) -> MutMap[Quantity]:
            """fetch model properties from a specific model"""
            return dict(model.properties.items())

        def fetch_thermo_props(model: ModelProxy) -> MutMap[Quantity]:
            """fetch properties of materials in a specific model"""
            ports = self.options["port_properties"]
            mat_proxy = model.materials
            return {k: v for k, v in mat_proxy.handler.items()
                    if ports or k not in mat_proxy}

        mod = self.model
        fetch = self.__fetch
        to_vector = self.__to_vector

        residual_structure = fetch(mod, lambda x: fetch_residuals(x, True),
                                   "normalised residual")
        bounds_structure = fetch(mod, fetch_bounds, "bound")
        residuals, residual_names = to_vector(residual_structure)
        bounds, bound_names = to_vector(bounds_structure)
        self.__vec_res_names[self.RES_VEC] = residual_names
        self.__vec_res_names[self.BOUND_VEC] = bound_names
        return {
            self.MODEL_PROPS: fetch(mod, fetch_mod_props, "model property"),
            self.THERMO_PROPS:
                fetch(mod, fetch_thermo_props, "thermo property"),
            self.RESIDUALS: fetch(mod, lambda x: fetch_residuals(x, False),
                                  "residual"),
            self.VECTORS: {
                self.RES_VEC: residuals,
                self.BOUND_VEC: bounds
            },
            self.JACOBIANS: {}
        }

        # The following jacobian is always needed
        # TODO: no it isn't!
        #   I might create 2 functions, one for a Newton step and one
        #   much cheaper to evaluate for line search.
        # maybe define flag in constructor whether to create this one right
        # away.

        # TODO: also might add one entry which is the mean square residual

        # if residuals.magnitude.rows() and states.magnitude.rows():
        #     self.__symres[cls.DR_DX] = jacobian(residuals, states)
        # else:
        #     self.__symres[cls.DR_DX] = Quantity(SX.sym("dr_dx", 0))
        #

    def __collect_argument_values(self) -> NestedMutMap[Quantity]:
        """Fetch initial states from materials, parameter values from
        thermo parameter stores, and parameter values from parameter handlers.
        """
        def fetch_states(model: ModelProxy) -> NestedMutMap[Quantity]:
            """Fetch the initial state variables from the materials of a
            specific model"""
            result = {}
            for k, m in model.materials.handler.items():
                init = m.initial_state
                try:
                    params = m.definition.store.get_all_values()
                except KeyError:
                    msg = "Missing values for thermodynamic parameters"
                    raise DataFlowError(msg)
                state = m.definition.frame.initial_state(init, params)
                dic = {f"x_{i:03d}": Quantity(x) for i, x in enumerate(state)}
                result[k] = dic
            return result

        def fetch_store_param() -> NestedMap[Quantity]:
            """fetch thermodynamic parameter values from the stores"""
            stores = NumericHandler.__fetch_thermo_stores(self.model)
            names = {store.name for store in stores}
            if len(names) < len(stores):
                raise ValueError("When using multiple ThermoPropertyStores, "
                                 "they have to have unique names")
            return {store.name: store.get_all_values() for store in stores}

        fetch = self.__fetch
        to_vector = self.__to_vector

        states = to_vector(fetch(self.model, fetch_states, "state"))[0]
        model_param = fetch(self.model, lambda m: m.parameters.values,
                            "parameter")
        return {
            self.VECTORS: {
                self.STATE_VEC: states,
            },
            self.MODEL_PARAMS: model_param,
            self.THERMO_PARAMS: fetch_store_param()
        }

    @staticmethod
    def __to_vector(struct: NestedMap[Quantity]) -> (Quantity, Sequence[str]):
        raw = [v.magnitude for v in flatten_dictionary(struct).values()]
        names = list(flatten_dictionary(struct).keys())
        return Quantity(vertcat(*raw)), names

    @staticmethod
    def __fetch(
            root: ModelProxy,
            func: Callable[[ModelProxy], NestedMutMap[Quantity]],
            typ: str,
            path: Optional[Sequence[str]] = None) -> NestedMutMap[Quantity]:
        """Drill recursively into child models to collect all data. The result
        is a nested dictionary, such that name clashes between child models and
        parameters are not permitted and will raise a ``ValueError``.
        """
        call_self = NumericHandler.__fetch
        if path is None:
            path = []
        result: NestedMutMap[Quantity] = func(root)
        for name, proxy in root.hierarchy.handler.items():
            if name in result:
                context = ".".join(path)
                msg = f"Child model / {typ} name clash:" \
                    f"'{name}' in {context}"
                raise ValueError(msg)
            result[name] = call_self(proxy, func, typ, path + [name])
        return result

    @staticmethod
    def __traverse(
            root: ModelProxy,
            func: Callable[[ModelProxy], None]):
        """Drill recursively into child modules to perform some action."""
        call_self = NumericHandler.__traverse
        func(root)
        for name, proxy in root.hierarchy.handler.items():
            call_self(proxy, func)

    @staticmethod
    def __fetch_thermo_stores(model: ModelProxy) \
            -> Collection[ThermoParameterStore]:
        call_self = NumericHandler.__fetch_thermo_stores
        result = {m.definition.store
                  for m in model.materials.handler.values()}
        for proxy in model.hierarchy.handler.values():
            result |= call_self(proxy)
        return result
