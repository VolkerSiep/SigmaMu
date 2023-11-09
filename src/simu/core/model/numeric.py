"""This module implements functionality concerning the numerical handling
of the top model instance."""

from typing import TYPE_CHECKING, Optional
from collections.abc import Callable, Sequence

from ..utilities.types import NestedMutMap, MutMap
from ..utilities.quantity import Quantity, QFunction
from ..utilities import flatten_dictionary

from .base import ModelProxy


class NumericHandler:
    """This class implements the function object describing the top level
    model."""

    def __init__(self, model: ModelProxy):
        self.model = model
        self.function = self.__make_function()
        self.arguments = self.__fetch_arguments()

    def __fetch_arguments(self):
        pass

    def __make_function(self):
        """Create a function that has the following arguments, each of them as
        a flat dictionary:

            - Model Parameters
            - Material States
            - Thermodynamic Parameters

        The result of the function consists likewise of

            - Properties
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
                result[name] = fetch(
                    proxy, func, typ, path + [name])
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

        mod = self.model
        args = {
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
