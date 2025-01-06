"""This module handles functionality concerning model hierarchy."""
from typing import TYPE_CHECKING, Type, Any
from collections.abc import Mapping, Iterator

from ..utilities.types import Map, MutMap
from ..utilities.errors import DataFlowError

if TYPE_CHECKING:  # avoid circular dependencies just for typing
    from .base import Model, ModelProxy


class HierarchyHandler(Map["ModelProxy"]):
    """This class, being instantiated as the :attr:`simu.Model.hierarchy`
    attribute, allows to define child models in a hierarchy context."""

    def __init__(self, model: "Model"):
        self.model = model
        self.__children: MutMap["ModelProxy"] = {}
        self.__declared: MutMap[Type["Model"]] = {}

    def __len__(self) -> int:
        return len(self.__children)

    def __iter__(self) -> Iterator[str]:
        return iter(self.__children)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        return

    def declare(self, name: str, model_cls: Type["Model"]) -> None:
        """Declare a sub-model in the interface, and by that

        a) Demand that it will be instantiated, and
        b) Make it available in the hierarchy proxy for browsing from parent
           context.
        """
        if name in self.__declared:
            raise KeyError(f"Child model '{name}' already declared")
        self.__declared[name] = model_cls

    def add(self, name: str, model_cls: Type["Model"],
            *args: Any, **kwargs: Any) -> "ModelProxy":
        """Add an instance of the class ``model_cls`` as child to the current
        (parent) context. A :class:`~simu.core.model.base.ModelProxy` object is
        created, registered, and returned."""
        if name in self.__children:
            raise KeyError(f"Child model '{name}' already exists")
        if (name in self.__declared and
                not issubclass(model_cls, self.__declared[name])):
            cls_name = self.__declared[name].__name__
            msg = f"Model '{name}' is not a subclass of {cls_name} as declared"
            raise ValueError(msg)

        instance = model_cls(*args, **kwargs).create_proxy(name)
        self.__children[name] = instance
        return instance

    @property
    def declared(self) -> Map[Type["Model"]]:
        """Dictionary of declared sub-models"""
        return self.__declared

    def __getitem__(self, name: str):
        """Re-obtain the proxy of named module, avoiding the need to keep a
        holding variable in the client scope code."""
        return self.__children[name]

    def create_proxy(self) -> "HierarchyProxy":
        """Create a proxy object for configuration in hierarchy context"""
        return HierarchyProxy(self)

    def check_complete(self):
        """Check that all declared child models are defined"""
        if missing := {n for n in self.__declared if n not in self.__children}:
            msg = "The following declared child models are not defined: " + \
                  ", ".join(missing)
            raise DataFlowError(msg)


class HierarchyProxy(Mapping[str, "ModelProxy"]):
    """A wrapper of the HierarchyHandler to grant access to the previously
    declared sub-models."""

    def __init__(self, handler: HierarchyHandler):
        self.handler = handler

    def __getitem__(self, name: str) -> "ModelProxy":
        if name not in self.handler.declared:
            raise KeyError(f"Child model of name 'name' not declared")
        return self.handler[name]

    def __len__(self) -> int:
        return len(self.handler.declared)

    def __iter__(self) -> Iterator[str]:
        return iter(self.handler.declared)

    def finalise(self):
        """Check that the child models declared in the interface are
        implemented"""
        self.handler.check_complete()
