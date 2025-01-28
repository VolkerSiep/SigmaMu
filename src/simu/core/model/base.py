"""This module contains the base classes to represent (process) models.

"""
from abc import ABC, abstractmethod
from typing import Self

from .parameter import ParameterHandler, ParameterProxy
from .hierarchy import HierarchyHandler, HierarchyProxy
from .property import PropertyHandler, PropertyProxy
from .material import MaterialHandler, MaterialProxy
from .bound import BoundHandler, BoundProxy
from simu.core.utilities.residual import ResidualHandler, ResidualProxy


class Model(ABC):
    """This is the base class for all process models to be implemented.

    By deriving from this class, handler objects are available as class
    attributes to deal with the particular aspects. Visit their documentation
    for details.

    The model implementation is then divided into two parts: (a) the interface,
    and (b) the model implementation itself. These parts are represented by the
    two methods :meth:`interface` and :meth:`define` respectively.

    A model is then either instantiated by a parent model, whereas it is
    represented in that parent model by a
    :class:`~simu.core.model.base.ModelProxy` object, or it can be defined as
    the top level model by calling :meth:`top`.
    """

    parameters: ParameterHandler
    """The handler object that takes care of parameter configuration"""

    properties: PropertyHandler
    """The handler object that takes care of property configuration"""

    hierarchy: HierarchyHandler
    """The handler object that takes care of defining sub models"""

    materials: MaterialHandler
    """The handler object that takes care of materials"""

    residuals: ResidualHandler = None
    """The handler object that takes care of residuals"""

    bounds: BoundHandler = None
    """The handler object that takes care of domain boundaries"""

    def __init__(self):
        """The constructor is parameterless but still needs to be called by the
         subclass constructors (if implemented) to initialise the object.
         Defining a constructor for the subclasses can be useful to pass custom
         data into the model.
         """
        self.__proxy = None
        self.parameters = ParameterHandler(self.cls_name)
        self.properties = PropertyHandler()
        self.materials = MaterialHandler()
        self.hierarchy = HierarchyHandler(self)
        self.interface()
        self.residuals = ResidualHandler()
        self.bounds = BoundHandler()

    @classmethod
    def top(cls, name: str = "model") -> "ModelProxy":
        """Define this model as top level model, hence instantiate,
        create proxy, and finalise it.

        This is the recommended and fastest way to declare a model as being the
        top level model. It will create a proxy object (via :meth:`proxy`) and
        finalise the configuration of it via
        :meth:`simu.core.model.base.ModelProxy.finalise`.

        Performing the finalisation in one go after creating the proxy object
        is only possible, if the model is suitable to be top model. That is,
        it must not have material ports or parameters to be connected.

        :param name: The name of the top level model
        :type name: str
        :return: The readily configured
          :class:`~simu.core.model.base.ModelProxy` object
        """
        return cls.proxy(name).finalise()

    @classmethod
    def proxy(cls, name: str = "model") -> "ModelProxy":
        """Instantiate and create proxy of this model.

        This is a helper method called by the
        :class:`~simu.core.model.hierarchy.HierarchyHandler`
        when defining a child model. After the proxy class is generated,
        the parent model can connect materials and provide parameters. When
        that is done, the proxy object can be finalised.

        :param name: The name of the top level model
        :type name: str
        :return: The :class:`~simu.core.model.base.ModelProxy` object, ready
          for configuration from parent context. It then still needs to be
          finalised.
        """
        return cls().create_proxy(name)

    def interface(self) -> None:
        """This virtual method is to define all model parameters, material
        ports, and properties provided by this model. This makes the interface
        of the model in the hierarchical context nearly self-documenting.
        A simple example implementation could be

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
        described in the :meth:`interface` method, a simple implementation
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

    def create_proxy(self, name: str = "model") -> "ModelProxy":
        """Create a proxy object for configuration in hierarchy context. This
        is the instance variant of the :meth:`proxy` class method.
        """
        return ModelProxy(self, name)

    @classmethod
    @property
    def cls_name(cls) -> str:
        """The name of the derived class with module path for hopefully
        unique identification.

        This method is mainly used to uniquely name static parameters of
        a model.

        :return: The name of the class with the module path
        :rtype: str
        """
        module_name = cls.__module__
        cls_name = cls.__name__
        if module_name == "__main__":
            return f"{cls_name}"
        else:
            return f"{module_name}.{cls_name}"


class ModelProxy:
    """Proxy class for models, being main object to relate to when accessing
    a child model during definition of the parent model.

    As for the :class:`Model` class, this proxy version offers via its handlers
    functionality to deal with parameters, properties, hierarchy, materials,
    and residuals, but the angle of view is different:

    The ``Model`` class deals with the implementation of the model, while the
    ``ModelProxy`` class offers its access to connect to it as a client. This
    client is most likely a parent ``Model`` or, if it is a top level model,
    a :class:`~simu.NumericHandler` object.

    ``ModelProxy`` objects are created from within the :class:`Model` class,
    and do not need to be instantiated directly by `SiMu` client code.
    """

    parameters: ParameterProxy
    """The proxy of the parameter handler, to connect and update parameters"""

    properties: PropertyProxy
    """The proxy of the property handler, making properties available"""

    hierarchy: HierarchyProxy
    """The proxy of the hierarchy handler, to parametrise child models"""

    materials: MaterialProxy
    """The proxy of the material handler, to connect material ports"""

    residuals: ResidualProxy
    """A handler object that takes care of residuals. This is really just a
    non-mutable mapping of residuals, as the client code is not supposed to
    temper with the definition of the child model. The residuals are still
    browsable, as required for instance by the :class:`~simu.NumericHandler`
    object for obvious reasons."""

    bounds: BoundProxy
    """The non-mutable proxy of the bound handler, to be queried by the 
    :class:`~simu.NumericHandler`."""

    def __init__(self, model: Model, name: str):
        self.parameters = model.parameters.create_proxy()
        self.properties = model.properties.create_proxy()
        self.hierarchy = model.hierarchy.create_proxy()
        self.materials = model.materials.create_proxy()
        self.residuals = model.residuals.create_proxy()
        self.bounds = model.bounds.create_proxy()
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
        this method is called to define and finalise itself.
        The final part of this process is to validate correct configuration and
        to clean up data structures.
        """
        self.parameters.finalise()  # parameters are final now
        self.materials.finalise()  # all ports are connected
        self.__model.define()
        self.properties.finalise()  # properties can be queried now
        self.hierarchy.finalise()  # all declared sub-models are provided
        return self
