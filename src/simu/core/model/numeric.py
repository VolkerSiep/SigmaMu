"""This module implements functionality concerning the numerical handling
of the top model instance."""

# std lib
from abc import ABC, abstractmethod
from typing import Optional, Any, Annotated
from collections.abc import Callable, Sequence, Collection
from enum import StrEnum, auto
from copy import deepcopy

# external
from casadi import vertcat, SX
from pint import Unit
from pint.registry import Quantity as QtyType
from pydantic import BaseModel, field_validator, Field, PlainValidator

# internal
from simu.core.utilities.quantity import Quantity, QFunction
from simu.core.utilities.structures import (
    flatten_dictionary, unflatten_dictionary, FLATTEN_SEPARATOR)
from simu.core.utilities.qstructures import (
    QuantityDict, quantity_dict_to_strings, parse_quantities_in_struct)
from simu.core.utilities.types import NestedMap, NestedMutMap, Map, MutMap
from simu.core.utilities.errors import DataFlowError
from simu.core.thermo.parameters import ThermoParameterStore
from simu.core.thermo.state import InitialState
from .base import ModelProxy


class PropertyFilter(ABC):
    def filter(self, properties: Map[Quantity | QuantityDict]) \
            -> Map[Quantity | QuantityDict]:
        """On filtering, this method receives the thermodynamic properties
        of a material as a mapping with keys as strings, representing the
        property name. The values are either scalar quantities or a
        quantity dictionary in case of non-scalar properties.
        """
        def filter_subkeys(key, sub_props: Quantity | QuantityDict):
            if isinstance(sub_props, Quantity):
                return sub_props
            else:
                return {
                    sub_key: value for sub_key, value in sub_props.items()
                    if self.keep_property(key, sub_key)
                }

        return {
            key: filter_subkeys(key, value) for key, value in properties.items()
            if self.keep_property(key)
        }

    @abstractmethod
    def keep_property(self, name: str, sub_key: str = None) -> bool:
        """Abstract method to decide whether a material property shall be
        included in the results of the process model.

        In case of non-scalar properties, the method is first called for the
        property itself without providing any ``sub_key``. Only if this call
        is answered with ``True``, the method is called again for each
        existing ``sub_key``.

        :param name: The name of the property
        :param sub_key: If the property is a non-scalar entity, the ``subkey``
          contains the identifier of the element, for instance the species name
          in case of mole flows or chemical potentials.
        """
        ...


class NHKeys(StrEnum):
    """Enumeration class to address sections in data structures related to
    the :class:`NumericHandler` class."""
    THERMO_PARAMS = auto()
    """Top level key in argument structure, addressing thermodynamic parameters.
    """

    MODEL_PARAMS = auto()
    """Top level key in argument structure, addressing process model parameters.
    """

    THERMO_PROPS = auto()
    """Top level key in result structure, addressing thermodynamic (or material)
    properties.
    """

    MODEL_PROPS = auto()
    """Top level key in result structure, addressing process model properties.
    """

    RESIDUALS = auto()
    """Top level key in result structure, addressing model residuals, and 
    sub-key in ``vectors`` section of result structure, containing a vector of
    dimensionless residuals, normalized by their tolerances.
    """
    STATES = auto()
    """Sub-key in ``vectors`` section of argument structure, containing a vector
    of thermodynamic state variables.
    """
    BOUNDS = auto()
    """Top level key in result structure, addressing model bounds, and 
    sub-key in ``vectors`` section of result structure, containing a lumped 
    vector of all bounds.
    """

    VECTORS = auto()
    """Top level key in both argument and result structure, pointing to
    vectorized data for efficient numerical treatment."""

    def __repr__(self):
        return f"'{self.value}'"


PQuantity = Annotated[str, PlainValidator(lambda v: Quantity(v))]


class SingleStateDump(BaseModel):
    T: PQuantity
    p: PQuantity
    n: Map[PQuantity]

    def to_dict(self) -> Map[QtyType]:
        return {"T": self.T, "p": self.p, "n": self.n}

    @field_validator("T", mode="after")
    @classmethod
    def check_temperature(cls, value: QtyType) -> QtyType:
        try:
            magnitude = value.to("K").m
        except Exception as e:
            raise ValueError(f"Invalid temperature: {value} - {e}")
        if magnitude <= 0:
            raise ValueError(f"Infeasible temperature value: {value}")
        return value

    @field_validator("p", mode="after")
    @classmethod
    def check_pressure(cls, value: QtyType) -> QtyType:
        try:
            magnitude = value.to("Pa").m
        except Exception as e:
            raise ValueError(f"Invalid Pressure: {value} - {e}")
        if magnitude <= 0:
            raise ValueError(f"Infeasible pressure value: {value}")
        return value

    @field_validator("n", mode="after")
    @classmethod
    def check_quantities(cls, value: Map[QtyType]) -> Map[QtyType]:
        for k, n_i in value.items():
            try:
                magnitude = n_i.to("mol").m
            except Exception as e:
                msg = f"Invalid Quantity for species {k}: {n_i} - {e}"
                raise ValueError(msg)
            if magnitude <= 0:
                msg = f"Infeasible quantity value for species {k}: {n_i}"
                raise ValueError(msg)
        return value


class StateDump(BaseModel):
    thermo: Map[Any]
    non_canonical: Map[Any] = Field(default=None)

    @field_validator("thermo", mode="before")
    @classmethod
    def validate_thermo(cls, value: Map[Any]) -> Map[Any]:
        def traverse(val):
            if set(val.keys()) == {"T", "p", "n"}:
                return SingleStateDump.model_validate(val)
            return {k: traverse(v) for k, v in val.items()}
        return traverse(value)

    @field_validator("non_canonical", mode="before")
    @classmethod
    def validate_non_canonical(cls, value: Map[Any] | None) -> Map[Any]:
        return {} if value is None else value


class NumericHandler:
    """This class implements the function object describing the top level
    model."""

    def __init__(self, model: ModelProxy, *,
                 property_filter: PropertyFilter = None,
                 port_properties: bool = False):
        """Create a numerical wrapper around a given model. This step is to be
        applied to any (top level) model that is to be numerically evaluated
        in any way (for solving, optimization, etc).

        :param model: The model to be wrapped. This model does not need to be
          square or well-posed. Such details are for the applied solvers to be
          fought with.
        :param property_filter: Larger models produce tens of thousands of
          properties. The house-keeping of those creates overhead internally,
          but also creates clutter for the client code. Applying a filter can
          help to limit the number of exported properties to a manageable level.
        :param port_properties: This parameter determines whether the properties
            of connected materials are also reported from a child model's
            perspective by the name of their ports. This is normally not
            interesting and thus off by default. In a generic front-end however,
            one might like to address a stream not only by its identifier in the
            containing context, but also via the port of a containing sub-model.
        """
        self.options = {
            "port_properties": port_properties
        }
        self.model = model
        self._property_filter = property_filter
        # the name vectors of vector arguments
        self.__vec_arg_names: MutMap[Sequence[str]] = {}
        self.__vec_res_names: MutMap[Sequence[str]] = {}

        # the symbolic argument structure
        self.__sym_args: NestedMutMap[Quantity] = self.__collect_arguments()
        # the symbolic result structure
        self.__sym_res: NestedMutMap[Quantity] = self.__collect_results()
        # the numerical argument structure with initial values
        self.__arguments: MutMap[Quantity] = {}

    @property
    def function(self) -> QFunction:
        """Create a Function object based on currently available argument
        and result structures."""
        return QFunction(self.__sym_args, self.__sym_res, "model")

    def vector_arg_names(self, key: str) -> Sequence[str]:
        """Return the names for the argument vector of given ``key``"""
        return self.__vec_arg_names[key]

    def vector_res_names(self, key: str) -> Sequence[str]:
        """Return the names for the result vector of given ``key``"""
        return self.__vec_res_names[key]

    @property
    def arguments(self) -> NestedMap[Quantity]:
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

        def traverse(model: ModelProxy, state_part: NestedMap[SingleStateDump],
                     path: str):
            # process local material objects
            all_names = set()
            for name, material in model.materials.handler.items():
                if name in model.materials:
                    continue
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
                        InitialState.from_dict(new_part.to_dict(),
                                               material.species)

            # traverse down into model hierarchy
            for name, proxy in model.hierarchy.handler.items():
                all_names.add(name)
                new_path = mk_new_path(path, name)
                new_part = state_part.get(name, {})
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
        state_valid = StateDump.model_validate(state)
        traverse(self.model, state_valid.thermo, "")
        return unflatten_dictionary(result)


    def retain_state(self, state: Sequence[float],
                     parameters: NestedMap[Quantity]):
        """Given a numeric ``state`` vector and the current set of
        ``parameters`` as a nested mapping of quantities, store the values for
        temperature, pressure and molar quantities back into the internal
        representations of the initial thermodynamic states.
        """
        def fetch_retain_initial_state(model: ModelProxy,
                                       states: NestedMap[float]):
            """retain initial states for a specific model"""
            mat_proxy = model.materials
            for k, m in mat_proxy.handler.items():
                if k in mat_proxy:
                    continue
                state_part = states[k].values()
                m.retain_initial_state(state_part, parameters)

        def traverse(model: ModelProxy, states: NestedMap[float]):
            fetch_retain_initial_state(model, states)
            for name, proxy in model.hierarchy.handler.items():
                if name in states:
                    traverse(proxy, states[name])

        state_struct = unflatten_dictionary(
            dict(zip(self.__vec_arg_names[NHKeys.STATES], state)))

        self.__arguments = {}  # force reread
        traverse(self.model, state_struct)

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

        if key in self.__sym_args[NHKeys.VECTORS]:
            msg = f"A parameter vector of name '{key}' is already used."
            raise KeyError(msg)

        if not self.__arguments:
            self.__arguments = self.__collect_argument_values()
        nam, sym, arg = traverse(definition, self.__sym_args, self.__arguments)

        result = Quantity(vertcat(*sym))
        self.__sym_args[NHKeys.VECTORS][key] = result
        values = Quantity(arg)
        self.__arguments[NHKeys.VECTORS][key] = values
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
            for k, item in items:
                n, s = traverse(item, symbols[k])
                if n is None:
                    nams.append(k)
                    syms.append(symbols[k].to(item).magnitude)
                else:
                    nams.extend([f"{k}/{n_i}" for n_i in n])
                    syms.extend(s)
            return nams, syms

        if key in self.__sym_res[NHKeys.VECTORS]:
            msg = f"A property vector of name '{key}' is already used."
            raise KeyError(msg)

        nam, sym = traverse(definition, self.__sym_res)

        result = Quantity(vertcat(*sym))
        self.__sym_res[NHKeys.VECTORS][key] = result
        self.__vec_res_names[key] = nam
        return result

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
        states,  state_names = to_vector(states_struct)
        self.__vec_arg_names[NHKeys.STATES] = state_names

        return {
            NHKeys.THERMO_PARAMS: fetch_store_param(mod),
            NHKeys.MODEL_PARAMS: fetch(mod, fetch_parameters, "parameter"),
            NHKeys.VECTORS: {
                NHKeys.STATES: states,
            }
        }

    def __collect_results(self) -> NestedMutMap[Quantity]:
        """The result of the function consists of

            - Model Properties
            - Thermodynamic (state) properties
            - Residuals
            - Bounds

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
            res.update(model.bounds.items())
            return res

        def fetch_mod_props(model: ModelProxy) -> MutMap[Quantity]:
            """fetch model properties from a specific model"""
            return dict(model.properties.items())

        def fetch_thermo_props(model: ModelProxy) -> MutMap[Quantity]:
            """fetch properties of materials in a specific model"""
            ports = self.options["port_properties"]
            mat_proxy = model.materials
            filter_ = self._property_filter
            f = (lambda x: x) if filter_ is None else filter_.filter
            return {k: f(v) for k, v in mat_proxy.handler.items()
                    if ports or k not in mat_proxy}

        mod = self.model
        fetch = self.__fetch
        to_vector = self.__to_vector

        residual_structure = fetch(mod, lambda x: fetch_residuals(x, True),
                                   "normalised residual")
        bounds_structure = fetch(mod, fetch_bounds, "bound")
        residuals, residual_names = to_vector(residual_structure)
        bounds, bound_names = to_vector(bounds_structure)
        self.__vec_res_names[NHKeys.RESIDUALS] = residual_names
        self.__vec_res_names[NHKeys.BOUNDS] = bound_names
        return {
            NHKeys.MODEL_PROPS:
                fetch(mod, fetch_mod_props, "model property"),
            NHKeys.THERMO_PROPS:
                fetch(mod, fetch_thermo_props, "thermo property"),
            NHKeys.RESIDUALS:
                fetch(mod, lambda x: fetch_residuals(x, False), "residual"),
            NHKeys.VECTORS: {
                NHKeys.RESIDUALS: residuals,
                NHKeys.BOUNDS: bounds
            }
        }

    def __collect_argument_values(self) -> NestedMutMap[Quantity]:
        """Fetch initial states from materials, parameter values from
        thermo parameter stores, and parameter values from parameter handlers.
        """
        def fetch_states(model: ModelProxy) -> NestedMutMap[Quantity]:
            """Fetch the initial state variables from the materials of a
            specific model"""
            result = {}
            mat_proxy = model.materials
            for k, m in mat_proxy.handler.items():
                if k in mat_proxy:
                    continue  # this is a connected port, don't collect twice

                init = m.initial_state
                frame = m.definition.frame
                param_struct = frame.parameter_structure
                try:
                    params = m.definition.store.get_values(param_struct)
                except KeyError:
                    msg = "Missing values for thermodynamic parameters"
                    raise DataFlowError(msg)
                state = frame.initial_state(init, params)
                # TODO: can I ask for proper state names from  frame?
                #  to do this, I had to get it from StateDefinition and add
                #  query functionality there.
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
            NHKeys.VECTORS: {
                NHKeys.STATES: states,
            },
            NHKeys.MODEL_PARAMS: model_param,
            NHKeys.THERMO_PARAMS: fetch_store_param()
        }

    @staticmethod
    def __to_vector(struct: NestedMap[Quantity]) -> (Quantity, Sequence[str]):
        flat = flatten_dictionary(struct)
        raw = [v.magnitude for v in flat.values()]
        return Quantity(vertcat(*raw)), list(flat.keys())

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
            res_i = call_self(proxy, func, typ, path + [name])
            if res_i:
                result[name] = res_i
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
