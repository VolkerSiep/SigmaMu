"""This module defines the abstract :class:`Model` class."""
from abc import ABC, abstractmethod

from ..utilities import QFunction
from .utils import ModelStatus
from .property import PropertyDefinition
from .parameter import ParameterDefinition
from .numeric import NumericHandler
# from .hierarchy import HierarchyHandler


class ModelInstance:
    def __init__(self, model: "Model"):
        self.__status = ModelStatus.READY

        # self.material = None  # TODO: material handler
        self.parameters = model.parameters.create_handler()
        self.properties = model.properties.create_handler()
        # self.hierarchy = None  # TODO: model.hierarchy.create_handler()
        # self.residuals = None  # TODO: residual handler

        self.numerics = None  # NumericHandler(self)

        self.__function = model.function  # by reference

    def finalise(self):
        """recreate the casadi nodes using the previously defined funtion.
        Use the parameters, child module properties, material objects, and
        locally defined non-canonical states to define a function,
        based on the :meth:`define` method's implementation.
        """
        # finalise child modules
        # for child in self.hierarchy.values():
        #     child.finalise()

        self.parameters.finalise()  # all required parameters connected?

        args = self.__collect_args()
        res = self.__function(args, squeeze_results=False)
        self.properties.finalise(res["properties"])

        self.__set_status(ModelStatus.FINALISED)

        return self

    def create_numerics(self) -> NumericHandler:
        """Considering this model instance to be top level, this method creates
        and returns a ``NumericHandler``.
        """
        self.status.assure(ModelStatus.FINALISED, "creating numeric handler")
        return NumericHandler(self)

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
        return self.__status

    def __set_status(self, status):
        self.__status = status
        self.parameters.status = status
        # TODO: add others


class Model(ABC):
    """The model is the central base class for process modelling.

    .. todo::

        much more documentation here!!
    """

    def __init__(self):
        self.__status = ModelStatus.INTERFACE
        self.material = None  # TODO: material handler
        self.parameters = ParameterDefinition()
        self.properties = PropertyDefinition()

        self.interface()

        self.__set_status(ModelStatus.DEFINE)

        self.hierarchy = None  # TODO: HierarchyHandler()
        self.residuals = None  # TODO: residual handler

        self.define()

        self.__set_status(ModelStatus.FUNCTION)

        self.__function = self.__make_function()

        self.__set_status(ModelStatus.READY)

    @property
    def status(self) -> ModelStatus:
        """Returns the current status of the model as a
        :class:`~simu.model.utils.ModelStatus` entity."""
        return self.__status

    def __set_status(self, status: ModelStatus):
        self.__status = status
        self.parameters.status = status
        self.properties.status = status
        # TODO: add status update to other handlers

    def instance(self, finalise: bool = False) -> ModelInstance:
        """Create an instance of the model. If ``finalise`` is set to ``True``,
        also finalise the instance right away. This is useful if no further
        connections are to be done.
        """
        instance = ModelInstance(self)
        if finalise:
            instance.finalise()
        return instance

    @property
    def I(self) -> ModelInstance:  # pylint: disable=invalid-name
        """Shortcut for ``instance()``"""
        return self.instance()

    @property
    def F(self) -> ModelInstance:  # pylint: disable=invalid-name
        """Shortcut for ``instance(finalise=True)``"""
        return self.instance(finalise=True)

    @property
    def function(self) -> QFunction:
        """The local casadi function that represents the calculations
        performed in the particular instance of this class."""
        return self.__function

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

    def __make_function(self):
        name = self.__class__.__name__
        args = self.__collect_args()
        results = {"properties": self.properties.local_symbols}
        return QFunction(args, results, name)

    def __collect_args(self):
        # child_properties = {
        #     name: child.properties.symbols
        #     for name, child in self.hierarchy.items()
        # }
        return {
            "parameters": self.parameters.symbols,
            # "properties": child_properties
        }




# class MyModel(Model):

#     def interface(self):
#         # require ports
#         self.material.require("inlet", "NitricAcid")
#         # define decorators as additional arguiments here!

#         # define process parameter
#         self.variables.define("pressure_drop", 0.3, "bar")  # with default
#         self.variables.define("heat_loss", unit="kW")  # no default

#         self.variables.provide("froth_density", unit="kg/m**3")

#     def define(self):
#         # need to define casadi graph that defines output from input
#         # base class has made SX source nodes ready upfront.
#         # How do I know UOM of input?
#         #   - parameters are fine
#         #   - material variables and child module variables:
#         #     from their definitions (QFunction.res_units) ?

#         # create own internal materials
#         self.material.create("intern", "NOxGas")

#         dp = self.variables["pressure_drop"]
#         T = self.material["inlet"]["T"]

#         self.variables["froth_density"] = dp * T  # not really

#         dp_calc = self.material["s301"]["p"] - self.material["intern"]["p"]
#         self.residuals.add("dp", dp_calc - dp, "5 bar")  # last arg is tolerance basis
