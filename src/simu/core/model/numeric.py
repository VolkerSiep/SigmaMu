"""This module implements functionality concerning the numerical handling
of the top model instance."""

from typing import TYPE_CHECKING, Optional
from collections.abc import Callable

from ..utilities.types import NestedMap, NestedMutMap
from ..utilities.quantity import Quantity

if TYPE_CHECKING:  # avoid circular dependencies just for typing
    from .base import ModelProxy


class NumericHandler:
    """This class implements the function object describing the top level
    model."""

    def __init__(self, model: "ModelProxy"):
        self.model = model
        self.__make_function()

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

        def fetch_in_hierarchy(
                root: "ModelProxy",
                func: Callable[["ModelProxy"], NestedMutMap[Quantity]],
                typ: str,
                path: Optional[list[str]] = None) -> NestedMutMap[Quantity]:
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
                result[name] = fetch_in_hierarchy(
                    proxy, func, typ, path + [name])
            return result

        parameters = fetch_in_hierarchy(
            self.model, lambda m: m.parameters.free, "parameter")
        args = {"parameters": parameters}

        properties = fetch_in_hierarchy(
            self.model, lambda m: m.properties.all, "property")
        results = {"properties": properties}

        # TODO:
        # make function, include also state, thermo param, residuals

        # query interface for all data:
        # all the following as quantity dictionaries:
        #   - export / import current states and parameters (thermo & model)
        #   - export current properties

        # in utilities, provide easy functions to serialise quantities or
        # dictionaries of them, e.g. "1.3 cm**2/s" should be parsed nicely for
        # yaml or json. Format with f"{a:~.16g}" to serialise.
