"""This module implements functionality related to parameter handling"""

from typing import Optional, TYPE_CHECKING
from collections.abc import KeysView

from ..utilities.quantity import (
    QuantityDict, SymbolQuantity, Quantity, NestedQuantityDict)
from ..utilities.errors import DataFlowError

if TYPE_CHECKING:  # avoid circular dependencies just for typing
    from .base import Model, ModelProxy

class ParameterHandler:
    """This class, being instantiated as the :attr:`Model.parameters`
    attribute, allows to define and access process parameters."""

    def __init__(self):
        self.__params: QuantityDict = {}
        self.__values: QuantityDict = {}

    def define(self,
               name: str,
               value: Optional[float] = None,
               unit: str = "dimless") -> None:
        """Define a parameter from within the :method:`Model.interface`
        method."""

        if name in self.__params:
            raise KeyError(f"Parameter '{name}' already defined")

        self.__params[name] = SymbolQuantity(name, unit)
        if value is not None:
            self.__values[name] = Quantity(value, unit)

    def __getitem__(self, name: str) -> Quantity:
        """Return the symbol of a defined parameter. This is to be used within
        the :method:`Model.define` method."""
        return self.__params[name]

    def names(self) -> KeysView[str]:
        """An iterator over all names of defined parameters"""
        return self.__params.keys()

    def values(self) -> QuantityDict:
        """Return the (default) values for the defined parameters that have
        such default value."""
        return dict(self.__values)  # create copy to preserve own set

    def create_proxy(self, name: str) -> "ParameterProxy":
        """Create a proxy object for configuration in hierarchy context"""
        return ParameterProxy(name, self)


    def collect_all(self, model: "Model") -> NestedQuantityDict:
        """Collect own parameters (all) and the yet free parameters of all
        child modules recursively."""
        params: NestedQuantityDict = dict(self.__params)

        def add_child_param(params: NestedQuantityDict,
                            name: str, child: "ModelProxy") -> None:
            if name in params:
                msg = f"Child model and parameter name clash: '{name}'"
                raise KeyError(msg)
            params[name] = collect(child)

        def collect(sub_model: "ModelProxy") -> NestedQuantityDict:
            params = sub_model.parameters.free
            for c_name, child in sub_model.hierarchy.items():
                add_child_param(params, c_name, child)
            return params

        for c_name, child in model.hierarchy.items():
            add_child_param(params, c_name, child)
        return params


class ParameterProxy:
    """This class is instantiated by the parent's :class:`ParameterHandler`
    to configure the parameter connections from the parent context."""

    def __init__(self, model_name: str, handler: ParameterHandler):
        self.model_name = model_name
        self.handler = handler

        self.__provided: QuantityDict = {}
        self.__free: QuantityDict = {}
        self.__values = handler.values()

    def provide(self, **kwargs: Quantity) -> None:
        """Connect a parameter from parent context to child parameter."""
        for name, quantity in kwargs.items():
            self.__assure_ok(name, quantity)
            self.__provided[name] = quantity

    def update(self, name: str, value: float, unit: str) -> None:
        """Provide new default value for child parameter"""
        quantity = Quantity(value, unit)
        self.__assure_ok(name, quantity)
        self.__values[name] = quantity

    def __assure_ok(self, name: str, quantity: Quantity) -> None:
        """Check that a new quantity can be processed and is compatible with
        the previously defined slot."""
        if name not in self.handler.names():
            msg = f"Parameter '{name}' not defined in '{self.model_name}'"
            raise KeyError(msg)
        if name in self.__provided:
            msg = f"Parameter '{name}' already provided in '{self.model_name}'"
            raise KeyError(msg)

        quantity.to(self.handler[name].units)  # complains if not compatible

    @property
    def free(self):
        """Symbols representing the parameters that have not been provided.
        These must be used to create an overall function, when parmeter values
        are provided from the outside."""
        return dict(self.__free)

    @property
    def arg(self) -> QuantityDict:
        """The function argument of the belonging model function, composed of
        the provided properties and the still free symbols."""
        arg = dict(self.__provided)
        arg.update(self.free)
        return arg

    def finalise(self) -> None:
        """Make sure there are values for all non-provided parameters with
        no value. Remove values of provided parameters"""

        # are all required connected?
        provided = set(self.__provided.keys())
        values = set(self.__values.keys())
        all_ = set(self.handler.names())
        missing = all_ - (provided | values)
        if missing:
            lst = ", ".join(map(lambda m: f"'{m}'", missing))
            msg = f"Model '{self.model_name}' has missing parameters: {lst}"
            raise DataFlowError(msg)

        # clean values of provided parameters
        self.__values = {
            name: quantity
            for name, quantity in self.__values.items() if name not in provided
        }

        # create symbols for still free variables (for later overall function)
        self.__free = {
            name: self.handler[name]
            for name in self.handler.names() if name in self.__values
        }
