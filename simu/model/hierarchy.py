"""This module facilitates model hierarchy.

    "*Hei, Ricky!*"
"""

from typing import TYPE_CHECKING
from .utils import ModelStatus

if TYPE_CHECKING:  # To avoid circular dependencies / need for abstract classes
    from .base import Model, ModelInstance

class HierarchyDefinition(dict):
    """The hierarchy definition keeps references to all child models as being a
    subclass of an ordinary dictionary.
    """

    def __init__(self):
        self.status = ModelStatus.INTERFACE

    def add(self, name: str, model: "Model") -> "ModelInstance":
        """Define an instance of the given ``model`` to be sub-model in the
        current context, named ``name``.
        """
        self.status.assure(ModelStatus.DEFINE, "defining child module")
        if name in self:
            raise KeyError(f"Sub-model of name {name} already defined")
        instance = model.instance()
        super().__setitem__(name, instance)
        return instance


    # def __setitem__(self, name: str, model: "Model"):
    #     """Define an instance of the given ``model`` to be sub-model in the
    #     current context, named ``name``.
    #     """
    #     self.status.assure(ModelStatus.DEFINE, "defining child module")
    #     if name in self:
    #         raise KeyError(f"Sub-model of name {name} already defined")
    #     instance = model.instance()
    #     super().__setitem__(name, instance)
    #     return instance

    def create_handler(self) -> "HierarchyHandler":
        self.status.assure(ModelStatus.READY, "defining child module")
        return HierarchyHandler(self)


class HierarchyHandler(dict):
    """The hierarchy handler creates instances for all child modules, when
    created from the ``HierachyDefinition`` object."""

    def __init__(self, definition: HierarchyDefinition):
        self.update(definition)
        self.status = ModelStatus.READY
