"""This module implements functionality concerning the numerical handling
of the top model instance."""

from typing import TYPE_CHECKING, Optional
from collections.abc import Callable, Sequence, Collection

from ..utilities.types import NestedMutMap, MutMap
from ..utilities.quantity import Quantity, QFunction
from ..utilities import flatten_dictionary

from .base import ModelProxy
from ..thermo import ThermoParameterStore


class NumericHandler:
    """This class implements the function object describing the top level
    model."""

    def __init__(self, model: ModelProxy):
        self.model = model
        self.function = self.__make_function()
        self.arguments = self.__fetch_arguments()

    def __fetch_arguments(self):
        return -1  # TODO: implement

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
        proxies. For child models, only the free parameters are to be
        collected.
        """
        def fetch(
                root: ModelProxy,
                func: Callable[[ModelProxy], NestedMutMap[Quantity]],
                typ: str,
                path: Optional[Sequence[str]] = None) -> NestedMutMap[Quantity]:
            """Drill recursively into child models to collect all free
            parameters. The result is a nested dictionary, such that name
            clashes between child models and parameters are not permitted
            and will raise a ``KeyError``.
            """
            if path is None:
                path = []
            result: NestedMutMap[Quantity] = func(root)
            for name, proxy in root.hierarchy.items():
                if name in result:
                    context = ".".join(path)
                    msg = f"Child model / {typ} name clash:" \
                        f"'{name}' in {context}"
                    raise KeyError(msg)
                result[name] = fetch(proxy, func, typ, path + [name])
            return result

        def fetch_thermo_stores(model: ModelProxy) \
                -> Collection[ThermoParameterStore]:
            result = {m.definition.store
                      for m in model.materials.handler.values()}
            for proxy in model.hierarchy.values():
                result |= fetch_thermo_stores(proxy)
            return result

        def fetch_residuals(model: ModelProxy) -> MutMap[Quantity]:
            return {k: r.value for k, r in model.residuals.items()}

        def fetch_material_states(model: ModelProxy) -> MutMap[Quantity]:
            return {k: m.sym_state for k, m in model.materials.handler.items()}

        def fetch_parameters(model: ModelProxy) -> MutMap[Quantity]:
            return dict(model.parameters.free)

        def fetch_mod_props(model: ModelProxy) -> MutMap[Quantity]:
            return model.properties

        def fetch_thermo_props(model: ModelProxy) -> MutMap[Quantity]:
            return model.materials.handler

        def fetch_store_props(model: ModelProxy) -> NestedMutMap[Quantity]:
            stores = fetch_thermo_stores(model)
            names = {store.name for store in stores}
            if len(names) < len(stores):
                raise ValueError("When using multiple ThermoPropertyStores, "
                                 "they have to have unique names")
            return {store.name: store.get_all_symbols() for store in stores}

        mod = self.model
        args = {
            "thermo_params": fetch_store_props(mod),
            "model_params": fetch(mod, fetch_parameters, "parameter"),
            "states": fetch(mod, fetch_material_states, "state"),
        }

        results = {
            "model_props": fetch(mod, fetch_mod_props, "model property"),
            "thermo_props": fetch(mod, fetch_thermo_props, "thermo property"),
            "residuals": fetch(mod, fetch_residuals, "residual")
        }

        results = flatten_dictionary(results)
        args = flatten_dictionary(args)

        return QFunction(args, results, "model")
