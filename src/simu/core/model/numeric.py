"""This module implements functionality concerning the numerical handling
of the top model instance."""

from typing import Optional
from collections.abc import Callable, Sequence, Collection
from casadi import vertcat, SX

from ..utilities import flatten_dictionary
from ..utilities.types import NestedMap, NestedMutMap, Map, MutMap
from ..utilities.quantity import Quantity, QFunction, jacobian
from ..utilities.errors import DataFlowError

from .base import ModelProxy
from ..thermo import ThermoParameterStore


class NumericHandler:
    """This class implements the function object describing the top level
    model."""
    THERMO_PARAMS: str = "thermo_params"
    MODEL_PARAMS: str = "model_params"
    STATE_VEC: str = "state_vector"
    MODEL_PROPS: str = "model_props"
    THERMO_PROPS: str = "thermo_props"
    RESIDUALS: str = "residuals"
    RES_VEC: str = "residual_vector"
    DR_DX: str = "dr_dx"

    function: QFunction

    def __init__(self, model: ModelProxy, port_properties=True):
        self.options = {
            "port_properties": port_properties
        }

        self.model = model
        self.function = self.__make_function()
        self.__arguments: Optional[MutMap[Quantity]] = None
        self.__symargs: Optional[NestedMap[Quantity]] = None
        self.__symres: Optional[NestedMap[Quantity]] = None

    def __fetch_arguments(self) -> NestedMutMap[Quantity]:
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

        cls = NumericHandler
        fetch = cls.__fetch
        to_vector = cls.__to_vector
        return {
            cls.STATE_VEC: to_vector(fetch(self.model, fetch_states, "state")),
            cls.MODEL_PARAMS: fetch(self.model, lambda m: m.parameters.values,
                                    "parameter"),
            cls.THERMO_PARAMS: fetch_store_param()
        }

    @property
    def arguments(self) -> Map[Quantity]:
        """The function arguments as numerical values. A DataFlowError is
        thrown, if not all numerical values are known."""
        if self.__arguments is None:
            self.__arguments = self.__fetch_arguments()
        return self.__arguments

    def __make_function(self):
        """Create a function that has the following arguments, each of them as
        a flat dictionary:

            - Material States
            - Model Parameters
            - Thermodynamic Parameters

        The result of the function consists likewise of

            - Model Properties
            - Thermodynamic (state) properties
            - Residuals

        All the data is to be collected from the model and all child model
        proxies. For child models, only the free parameters are
        collected.
        """

        def fetch_residuals(model: ModelProxy) -> MutMap[Quantity]:
            """fetch residuals from a specific model"""
            return {k: r.value for k, r in model.residuals.items()}

        def fetch_normed_residuals(model: ModelProxy) -> MutMap[Quantity]:
            """fetch normed residuals from a specific model"""
            return {k: (r.value / r.tolerance).to("")
                    for k, r in model.residuals.items()}

        def fetch_material_states(model: ModelProxy) -> MutMap[Quantity]:
            """fetch material states from a specific model"""
            mat_proxy = model.materials
            return {k: m.sym_state for k, m in mat_proxy.handler.items()
                    if k not in mat_proxy}

        def fetch_parameters(model: ModelProxy) -> MutMap[Quantity]:
            """fetch model parameters from a specific model"""
            return dict(model.parameters.free)

        def fetch_mod_props(model: ModelProxy) -> MutMap[Quantity]:
            """fetch model properties from a specific model"""
            return dict(model.properties)

        def fetch_thermo_props(model: ModelProxy) -> MutMap[Quantity]:
            """fetch properties of materials in a specific model"""
            ports = self.options["port_properties"]
            mat_proxy = model.materials
            return {k: v for k, v in mat_proxy.handler.items()
                    if ports or k not in mat_proxy}

        def fetch_store_param(model: ModelProxy) -> NestedMap[Quantity]:
            """fetch thermodynamic parameters from the stores"""
            stores = cls.__fetch_thermo_stores(model)
            names = {store.name for store in stores}
            if len(names) < len(stores):
                raise ValueError("When using multiple ThermoPropertyStores, "
                                 "they have to have unique names")
            return {store.name: store.get_all_symbols() for store in stores}

        mod = self.model
        cls = NumericHandler
        fetch = cls.__fetch
        to_vector = cls.__to_vector

        states = to_vector(fetch(mod, fetch_material_states, "state"))
        self.__symargs = {
            cls.THERMO_PARAMS: fetch_store_param(mod),
            cls.MODEL_PARAMS: fetch(mod, fetch_parameters, "parameter"),
            cls.STATE_VEC: states,
        }

        residuals = to_vector(fetch(mod, fetch_normed_residuals,
                                    "normalised residual"))
        self.__symres = {
            cls.MODEL_PROPS: fetch(mod, fetch_mod_props, "model property"),
            cls.THERMO_PROPS:
                fetch(mod, fetch_thermo_props, "thermo property"),
            cls.RESIDUALS: fetch(mod, fetch_residuals, "residual"),
            cls.RES_VEC: residuals
        }
        # The following jacobian is always needed
        if residuals.magnitude.rows() and states.magnitude.rows():
            self.__symres[cls.DR_DX] = jacobian(residuals, states)
        else:
            self.__symres[cls.DR_DX] = Quantity(SX.sym("dr_dx", 0))

        return QFunction(self.__symargs, self.__symres, "model")

    # TODO: implement modifying function to also deliver Jacobians
    #  Easy: dr/dx, as we always need all of both groups
    #  Challenge: How to address single thermodynamic or model parameters
    #   or properties
    #  maybe give as argument the entire structures, which then are
    #  sub-structures of the argument and result structures.

    @staticmethod
    def __to_vector(struct: NestedMap[Quantity]) -> Quantity:
        raw = [v.magnitude for v in flatten_dictionary(struct).values()]
        return Quantity(vertcat(*raw))

    @staticmethod
    def __fetch(
            root: ModelProxy,
            func: Callable[[ModelProxy], NestedMutMap[Quantity]],
            typ: str,
            path: Optional[Sequence[str]] = None) -> NestedMutMap[Quantity]:
        """Drill recursively into child models to collect all data. The result
        is a nested dictionary, such that name clashes between child models and
        parameters are not permitted and will raise a ``KeyError``.
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
                raise KeyError(msg)
            result[name] = call_self(proxy, func, typ, path + [name])
        return result

    @staticmethod
    def __fetch_thermo_stores(model: ModelProxy) \
            -> Collection[ThermoParameterStore]:
        call_self = NumericHandler.__fetch_thermo_stores
        result = {m.definition.store
                  for m in model.materials.handler.values()}
        for proxy in model.hierarchy.handler.values():
            result |= call_self(proxy)
        return result
