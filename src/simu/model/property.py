"""This module implements functionality related to property handling"""

from collections.abc import KeysView

from ..utilities.types import QuantityDict
from ..utilities.quantity import Quantity
from ..utilities.errors import DataFlowError


class PropertyHandler:
    """This class, being instantiated as the :attr:`Model.properties`
    attribute, allows to declare and define process properties."""

    def __init__(self):
        self.__props: QuantityDict = {}
        self.__defined: QuantityDict = {}

    def declare(self, name: str, unit: str):
        """This method declares a property to be provided by the model."""
        if name in self.__defined:
            raise DataFlowError(f"Property '{name}' is already defined")
        self.__defined[name] = Quantity(unit)  # also assure the unit is valid

    def __setitem__(self, name: str, quantity: Quantity):
        """Via this operator, calculated properties can be registered as a
        model result."""
        if name not in self.__defined:
            raise KeyError(f"Property '{name}' is not declared")
        if name in self.__props:
            raise DataFlowError(f"Property '{name}' is already provided")
        self.__props[name] = quantity

    def __getitem__(self, name: str) -> Quantity:
        """Return a property as it has been registered via the ``__setitem__``
        operator (``[]``) before."""
        return self.__props[name]

    def names(self) -> KeysView[str]:
        """An iterator over all names of defined properties"""
        return self.__props.keys()

    def create_proxy(self, name: str) -> "PropertyProxy":
        """Create a proxy object for configuration in hierarchy context"""
        return PropertyProxy(name, self)


class PropertyProxy:
    """This class is instantiated by the parent's :class:`PropertyHandler`
    to allow property access to the parent context."""

    def __init__(self, name: str, handler: PropertyHandler):
        self.model_name = name
        self.handler = handler
        self.__wired_props: QuantityDict = {}

    def __getitem__(self, name: str) -> Quantity:
        """After the proxy has left the context manager, it is finalised, and
        the calculated properties can be queried using this operator."""
        return self.__wired_props[name]

    def finalise(self, wired_props: QuantityDict) -> None:
        """This method is called from the :meth:`ModelProxy.finalise` method,
        typically when leaving the context manager. Afterwards, the calculated
        properties can be used by the parent model via the ``__getitem__``
        operator (``[]``)"""
        self.__wired_props = wired_props

    @property
    def all(self) -> QuantityDict:
        """A dictionary of all properties."""
        return dict(self.__wired_props)
