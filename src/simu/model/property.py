"""This module implements functionality related to property handling"""

from collections.abc import Mapping
from typing import Iterator

from ..utilities.types import QuantityDict
from ..utilities.quantity import Quantity
from ..utilities.errors import DataFlowError


class PropertyHandler(Mapping[str, Quantity]):
    """This class, being instantiated as the :attr:`Model.properties`
    attribute, allows to declare and define process properties."""

    def __init__(self):
        self.__model_name = "N/A"
        self.__props: QuantityDict = {}
        self.__declared: QuantityDict = {}

    def __iter__(self) -> Iterator[str]:
        return iter(self.__props)

    def __len__(self) -> int:
        return len(self.__props)

    @staticmethod
    def __raise(name: str, msg: str):
        raise DataFlowError(f"Property '{name}' {msg}")

    def __setitem__(self, name: str, quantity: Quantity):
        """Via this operator, a calculated property is defined as a
        model result."""
        if name not in self.__declared:
            self.__raise(name,  "is not declared")
        if name in self.__props:
            self.__raise(name,  "is already defined")

        # check unit compatibility -> throw exception on demand
        quantity.to(self.__declared[name])

        self.__props[name] = quantity

    def __getitem__(self, name: str) -> Quantity:
        """Return a property as it has been defined via the ``__setitem__``
        operator (``[]``) before."""
        return self.__props[name]

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        return

    def declare(self, name: str, unit: str):
        """This method declares a property to be provided by the model."""
        if name in self.__declared:
            self.__raise(name,  "is already declared")
        self.__declared[name] = Quantity(unit)  # also assure the unit is valid

    def check_complete(self):
        """Check that all declared properties are defined"""
        if missing := {n for n in self.__declared if n not in self.__props}:
            msg = "The following properties are not defined: " + \
                  ", ".join(missing)
            raise DataFlowError(msg)

    def create_proxy(self) -> "PropertyProxy":
        """Create a proxy object for configuration in hierarchy context"""
        return PropertyProxy(self)


class PropertyProxy(Mapping[str, Quantity]):
    """This class is instantiated by the parent's :class:`PropertyHanler`
    to handle the property availability to the parent context."""

    def __init__(self, handler: PropertyHandler):
        self.__handler = handler
        self.__finalised = False

    def __getitem__(self, name: str) -> Quantity:
        if self.__finalised:
            return self.__handler[name]
        else:
            raise DataFlowError("Property access before model is finalised")

    def __len__(self) -> int:
        if self.__finalised:
            return len(self.__handler)
        else:
            raise DataFlowError("Property access before model is finalised")

    def __iter__(self) -> Iterator[str]:
        if self.__finalised:
            return iter(self.__handler)
        else:
            raise DataFlowError("Property access before model is finalised")

    def finalise(self):
        """Signal that the model is done with its declaration, and properties
        are now available for the target context."""
        self.__handler.check_complete()
        self.__finalised = True
