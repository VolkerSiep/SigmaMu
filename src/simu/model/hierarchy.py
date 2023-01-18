"""This module handles functionality concerning model hierarchy."""
from typing import TYPE_CHECKING
from collections.abc import ItemsView

if TYPE_CHECKING:  # avoid circular dependencies just for typing
    from .base import Model, ModelProxy

ModelProxyDictionary = dict[str, "ModelProxy"]


class HierarchyHandler:
    """This class, being instantiated as the :attr:`Model.hierarchy` attribute,
    allows to define child models in a hierarchy context."""

    def __init__(self, model: "Model"):
        self.model = model
        self.__children: ModelProxyDictionary = {}

    def add(self, name: str, model: "Model") -> "ModelProxy":
        """Add an instance of ``model`` as child to the current (parent)
        context. A :class:`ModelProxy` object is created, registered, and
        returned."""
        if name in self.__children:
            raise KeyError(f"Child model '{name}' already exists")
        instance = model.create_proxy(name)
        self.__children[name] = instance
        return instance

    def __getitem__(self, name: str):
        """Re-obtain the proxy of named module, avoiding to have to keep a
        holding variable in the client scope code."""
        return self.__children[name]

    def items(self) -> ItemsView[str, "ModelProxy"]:
        """Return an iterator over the dictionary of child module proxies."""
        return self.__children.items()
