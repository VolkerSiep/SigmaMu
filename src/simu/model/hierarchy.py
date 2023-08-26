"""This module handles functionality concerning model hierarchy."""
from typing import TYPE_CHECKING, Iterator, Type, Any
from collections.abc import Mapping

if TYPE_CHECKING:  # avoid circular dependencies just for typing
    from .base import Model, ModelProxy

ModelProxyDictionary = dict[str, "ModelProxy"]


class HierarchyHandler(Mapping):
    """This class, being instantiated as the :attr:`Model.hierarchy` attribute,
    allows to define child models in a hierarchy context."""

    def __init__(self, model: "Model"):
        self.model = model
        self.__children: ModelProxyDictionary = {}

    def __len__(self) -> int:
        return len(self.__children)

    def __iter__(self) -> Iterator[str]:
        return iter(self.__children)

    def add(self, name: str, model_cls: Type["Model"],
            *args: Any, **kwargs: Any) -> "ModelProxy":
        """Add an instance of the class ``model_cls`` as child to the current
        (parent) context. A :class:`ModelProxy` object is created, registered,
        and returned."""
        if name in self.__children:
            raise KeyError(f"Child model '{name}' already exists")
        instance = model_cls(*args, **kwargs).create_proxy(name)
        self.__children[name] = instance
        return instance

    def __getitem__(self, name: str):
        """Re-obtain the proxy of named module, avoiding to have to keep a
        holding variable in the client scope code."""
        return self.__children[name]
