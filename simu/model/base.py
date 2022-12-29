"""This module defines the abstract :class:`Model` class."""
from abc import ABC, abstractmethod

from ..utilities import QFunction
from .utils import ModelStatus
from .property import PropertyDefinition, PropertyHandler
from .parameter import ParameterDefinition, ParameterHandler
from .numeric import NumericHandler
from .hierarchy import HierarchyDefinition, HierarchyHandler


class ModelInstance:
    """When a :class:`Model` is created, the initial configuration is common
    for all later instances, and it creates one function object to represent
    itself. A model instance object - as described by this class - is an object
    specifically used once as part of the overall (top-level) model.

    As an example, the ``Model`` object may describe a tray in a column. Then,
    one ``ModelInstance`` object is created and used for each tray.

    An end-user of ``simu``, normally does not need to relate to this class
    directly.
    """

    parameters: ParameterHandler
    """The object to handle all parameter related functionality."""

    properties: PropertyHandler
    """The object to handle all property related functionality."""

    hierarchy: HierarchyHandler
    """The object to handle all child-module related functionality."""

    numerics: NumericHandler
    """The object to handle all functionality related to performing actual
    calculations."""

    def __init__(self, model: "Model"):
        self.__status = ModelStatus.READY

        # self.material = None  # TODO: material handler
        self.parameters = model.parameters.create_handler()
        self.properties = model.properties.create_handler()
        self.hierarchy = model.hierarchy.create_handler()
        # self.residuals = None  # TODO: residual handler

        self.numerics = NumericHandler(self)

        self.__function = model.function  # by reference

    def finalise(self) -> "ModelInstance":
        """recreate the casadi nodes using the previously defined funtion.
        Use the parameters, child module properties, material objects, and
        locally defined non-canonical states to define a function,
        based on the :meth:`define` method's implementation.
        """
        self.status.assure(ModelStatus.READY, "finalising model")
        # finalise child modules
        for child in self.hierarchy.values():
            child.finalise()

        self.parameters.finalise()  # all required parameters connected?

        args = self.__collect_args()
        res = self.__function(args, squeeze_results=False)
        self.properties.finalise(res["properties"])

        self.__set_status(ModelStatus.FINALISED)
        return self

    def __collect_args(self):
        # child_properties = {
        #     name: child.properties.symbols
        #     for name, child in self.hierarchy.items()
        # }
        return {
            "parameters": self.parameters.symbols,
            # "properties": child_properties
        }

    @property
    def status(self) -> ModelStatus:
        """Returns the current status of the model instance as a
        :class:`~simu.model.utils.ModelStatus` entity."""
        return self.__status

    def __set_status(self, status: ModelStatus):
        self.__status = status
        self.parameters.status = status
        self.hierarchy.status = status


class Model(ABC):
    """The model is the central base class for process modelling.

    .. todo::

        much more documentation here!!
    """
    parameters: ParameterDefinition
    """The object to handle all parameter related functionality."""

    properties: PropertyDefinition
    """The object to handle all property related functionality."""

    hierarchy: HierarchyDefinition
    """The object to handle all child-module related functionality."""

    def __init__(self):
        # interface phase
        self.__status = ModelStatus.INTERFACE
        self.material = None  # TODO: material handler
        self.parameters = ParameterDefinition()
        self.properties = PropertyDefinition()
        self.interface()

        # definition phase
        self.hierarchy = HierarchyDefinition()
        self.residuals = None  # TODO: residual handler
        self.__set_status(ModelStatus.DEFINE)
        self.define()

        # function phase
        self.__set_status(ModelStatus.FUNCTION)
        self.__function = self.__make_function()

        # ready to be instantiated
        self.__set_status(ModelStatus.READY)

    @abstractmethod
    def interface(self) -> None:
        '''This abstract method is to define model parameters, material ports,
        and properties provided by this model. This makes the interface of
        the model in the hierarchical context nearly self-documenting.
        A trivial example implementation could be

        .. code-block::

            def interface(self):
                """Here a nice documentation of the model inteface"""
                self.parameters.define("length", 10.0, "m")
                self.properties.provide("area", unit="m**2")

        Above interface requires a parameter called ``length`` with a default
        value of 10 metres. It promises to calculate a property called ``area``
        which has a unit compatible to square metres.
        '''

    @abstractmethod
    def define(self) -> None:
        '''This abstract method is to define the model implementation,
        including the use of sub-modules, creation of internal materials, and
        calculation of residuals and model properties. Matching to the example
        described in the :meth:`interface` method, a trivial implementation
        could be

        .. code-block::

            def define(self):
                """Here documentation of the internal function of the model.
                This distinction can be used to include this doc-string only
                for detailed documenation sections."""
                length = self.parameters["length"]
                self.properties["area"] = length * length

        Here we read out the previously defined parameter ``length`` and
        calculate the property ``area``.
        '''

    @property
    def status(self) -> ModelStatus:
        """Returns the current status of the model as a
        :class:`~simu.model.utils.ModelStatus` entity."""
        return self.__status

    def __set_status(self, status: ModelStatus):
        self.__status = status
        self.parameters.status = status
        self.properties.status = status
        self.hierarchy.status = status
        # TODO: add status update to other handlers

    def instance(self) -> ModelInstance:
        """Create an instance of the model. This is normally done either from
        within the :class:`HierarchyHandler` for child models or in the
        :attr:`N` property for the top level model."""
        return ModelInstance(self)

    @property
    def N(self) -> NumericHandler:  # pylint: disable=invalid-name
        """Shortcut for ``self.instance().finalise().numerics.prepare()``

        This need letter does a couple of things required for
        using the top level model:

            - Instantiate the model to a :class:`ModelInstance` object
            - Finalise the instance object
            - query the :class:`NumericHandler` object
            - prepare the numeric handler as top level model
            - return the numeric handler

        Typical top level code is::

            class MyTopLevelModel(Model):
                def interface(self):
                    ...
                def define(self):
                    ...

            def main():
                num_obj = MyTopLevelModel().N
                ... # do some actual calculations
        """
        return self.instance().finalise().numerics.prepare()

    @property
    def function(self) -> QFunction:
        """The local casadi function that represents the calculations
        performed in the particular instance of this class."""
        return self.__function

    def __make_function(self):
        name = self.__class__.__name__
        args = {
            "parameters": self.parameters.symbols,
            "properties": {
                name: child.properties.symbols
                for name, child in self.hierarchy.items()
            }
        }
        results = {"properties": self.properties.local_symbols}
        return QFunction(args, results, name)
