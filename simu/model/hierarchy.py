"""This module handles functionality concerning model hierarchy."""
from typing import TYPE_CHECKING
from collections.abc import ItemsView

if TYPE_CHECKING:  # avoid circular dependencies just for typing
    from .base import Model, ModelProxy

ModelProxyDictionary = dict[str, "ModelProxy"]


class HierarchyHandler:
    """This class, being instantiated as the :attr:`Model.hierachy` attribute,
    allows to define child models in a hierachy context."""

    def __init__(self, model: "Model"):
        self.model = model
        self.__childs: ModelProxyDictionary = {}

    def add(self, name: str, model: "Model") -> "ModelProxy":
        """Add an instance of ``model`` as child to the current (parent)
        context. A :class:`ModelProxy` object is created, registered, and
        returned."""
        if name in self.__childs:
            raise KeyError(f"Child model '{name}' already exists")
        instance = model.create_proxy(name)
        self.__childs[name] = instance
        return instance

    def items(self) -> ItemsView[str, "ModelProxy"]:
        """Return an iterator over the dictionary of child module proxies."""
        return self.__childs.items()
