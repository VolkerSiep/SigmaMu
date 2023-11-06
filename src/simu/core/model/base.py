"""This module contains the base classes to represent (process) models."""
from abc import ABC, abstractmethod
from typing import Self, Optional

# from ..utilities import QFunction

from .parameter import ParameterHandler, ParameterProxy
from .hierarchy import HierarchyHandler, HierarchyProxy
from .property import PropertyHandler, PropertyProxy
from .material import MaterialHandler, MaterialProxy
from .residual import ResidualHandler
# from .numeric import NumericHandler


class Model(ABC):
    """This is the base class for all process models to be implemented."""

    parameters: ParameterHandler
    """A handler object that takes care of parameter configuration"""

    properties: PropertyHandler
    """A handler object that takes care of property configuration"""

    hierarchy: HierarchyHandler
    """A handler object that takes care of defining sub models"""

    materials: MaterialHandler
    """A handler object that takes care of materials"""

    residuals: ResidualHandler = None
    """A handler object that takes care of residuals"""

    def __init__(self):
        self.__proxy = None
        self.parameters = ParameterHandler(self.cls_name)
        self.properties = PropertyHandler()
        self.materials = MaterialHandler()
        self.hierarchy = HierarchyHandler(self)
        self.interface()
        self.residuals = ResidualHandler()

    @classmethod
    def top(cls, name: str = "model") -> "ModelProxy":
        """Define this model as top level model, hence instantiate,
        create proxy, and finalise it.
        """
        return cls.proxy(name).finalise()

    @classmethod
    def proxy(cls, name: str = "model") -> "ModelProxy":
        """Instantiate and create proxy of this model."""
        return cls().create_proxy(name)

    def interface(self) -> None:
        """This abstract method is to define all model parameters, material
        ports, and properties provided by this model. This makes the interface
        of the model in the hierarchical context nearly self-documenting.
        A trivial example implementation could be

        .. code-block::

            def interface(self):
                \"\"\"Here a nice documentation of the model interface\"\"\"
                self.parameters.define("length", 10.0, "m")
                self.properties.provide("area", unit="m**2")

        Above interface requires a parameter called ``length`` with a default
        value of 10 metres. It promises to calculate a property called ``area``
        which has a unit compatible to square metres.
        """

    @abstractmethod
    def define(self) -> None:
        """This abstract method is to define the model implementation,
        including the use of submodules, creation of internal materials, and
        calculation of residuals and model properties. Matching to the example
        described in the :meth:`interface` method, a trivial implementation
        could be

        .. code-block::

            def define(self):
                \"\"\"Here documentation of the internal function of the model.
                This distinction can be used to include this doc-string only
                for detailed documentation sections.\"\"\"
                length = self.parameters["length"]
                self.properties["area"] = length * length

        Here we read out the previously defined parameter ``length`` and
        calculate the property ``area``.
        """

    # following the two methods finalise and create_proxy for explicit use,
    # if the context manager is not used.

    def create_proxy(self, name: Optional[str] = "model") -> "ModelProxy":
        """Create a proxy object for configuration in hierarchy context"""
        return ModelProxy(self, name)

    @classmethod
    @property
    def cls_name(cls) -> str:
        """The name of the derived class with module path for hopefully
        unique identification"""
        module_name = cls.__module__
        cls_name = cls.__name__
        if module_name == "__main__":
            return f"{cls_name}"
        else:
            return f"{module_name}.{cls_name}"


class ModelProxy:
    """Proxy class for models, being main object to relate to when accessing
    a child model during definition of the parent model.
    """

    parameters: ParameterProxy
    """The proxy of the parameter handler, to connect and update parameters"""

    properties: PropertyProxy
    """The proxy of the property handler, making properties available"""

    hierarchy: HierarchyProxy
    """The proxy of the hierarchy handler, to parametrise child models"""

    materials: MaterialProxy
    """The proxy of the material handler, to connect material ports"""

    def __init__(self, model: Model, name: str):
        self.parameters = model.parameters.create_proxy()
        self.properties = model.properties.create_proxy()
        self.hierarchy = model.hierarchy.create_proxy()
        self.materials = model.materials.create_proxy()
        self.__model = model
        self.__set_name(name)

    def __set_name(self, name):
        """Set the name of the model for better error diagnostics"""
        self.parameters.set_name(name)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        if exc_type is not None:
            return
        self.finalise()

    def finalise(self) -> Self:
        """After the parent model has connected parameters and materials,
        this method is called to process its own modelling code"""
        self.parameters.finalise()  # parameters are final now
        self.materials.finalise()  # all ports are connected
        self.__model.define()
        self.properties.finalise()  # properties can be queried now
        self.hierarchy.finalise()  # all declared sub-models are provided
        return self

    # TODO: method that allows the numeric handler to collect the symbols
    #       and values.
