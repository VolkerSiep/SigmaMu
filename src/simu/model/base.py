"""This module contains the base classes to represent (process) models."""
from abc import ABC, abstractmethod

from ..utilities import QFunction

from .parameter import ParameterHandler
from .hierarchy import HierarchyHandler
from .property import PropertyHandler
from .numeric import NumericHandler


class Model(ABC):
    """This is the base class for all process models to be implemented."""

    parameters: ParameterHandler
    """A handler object that takes care of parameter configuration"""

    properties: PropertyHandler
    """A handler object that takes care of property configuration"""

    hierarchy: HierarchyHandler
    """A handler object that takes care of defining sub models"""

    @classmethod
    def as_top_model(cls,
                     name: str = "model"
                     ) -> tuple["ModelProxy", NumericHandler]:
        """Define this model as top level model, hence instantiate and finalise
        it.
        """
        proxy = cls().create_proxy(name).finalise()
        return proxy, proxy.numeric

    def __init__(self):
        self.parameters = ParameterHandler()
        self.properties = PropertyHandler()

        self.interface()
        self.hierarchy = HierarchyHandler(self)
        self.define()
        self.__func = self.__make_function()

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

    def create_proxy(self, name: str) -> "ModelProxy":
        """Create a proxy object for configuration in hierarchy context"""

        return ModelProxy(name, self)

    def __make_function(self) -> QFunction:
        """This method creates a function object that also includes the child
        models.
        """
        name = self.__class__.__name__
        args = {"parameters": self.parameters.all}
        props = {n: self.properties[n] for n in self.properties.names()}

        res = {"properties": props}
        return QFunction(args, res, name)

    @property
    def function(self) -> QFunction:
        """Return the complete function object of this module."""
        return self.__func


class ModelProxy:
    """Proxy class for models, being main object to relate to when accessing
    a child model during definition of the parent model.
    """

    def __init__(self, name: str, model: Model):
        """Constructor creating proxy objects for the handlers of the model,
        such that a sub-model can be configured in the parent context."""

        self.model = model
        self.parameters = model.parameters.create_proxy(name)
        self.properties = model.properties.create_proxy(name)
        self.hierarchy = model.hierarchy  # they are already proxies.

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        if exc_type is not None:
            return
        self.finalise()

    def finalise(self) -> "ModelProxy":
        """This method is typically called when the object leaves the context
        manager. It calls the child model's function to wire up the calculated
        properties and residuals."""
        # call model function on not provided parameters,
        # and calculate properties
        self.parameters.finalise()

        args = {"parameters": self.parameters.arg}
        result = self.model.function(args, squeeze_results=False)

        # There might be no properties, e.g. if model only makes residuals
        self.properties.finalise(result.get("properties", {}))

        return self

    @property
    def numeric(self) -> NumericHandler:
        """Create and return a numeric handler"""
        return NumericHandler(self)
