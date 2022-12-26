"""This module defines the abstract :class:`Model` class."""
from abc import ABC, abstractmethod

from ..utilities import QFunction
from .property import PropertyHandler
from .parameter import ParameterDefinition
from .numeric import NumericHandler
from .hierarchy import HierarchyHandler


class Model(ABC):
    """The model is the central base class for process modelling.

    .. todo::

        much more documentation here!!
    """

    def __init__(self):
        self.material = None  # TODO: material handler
        self.parameters = ParameterDefinition()
        self.properties = PropertyHandler()

        self.interface()

        self.hierarchy = None  # TODO: HierarchyHandler()
        self.residuals = None  # TODO: residual handler

        self.define()

        self.__function = self.__make_function()

    # def __enter__(self):
    #     # prepare an object that can connect the actual symbols and ports.
    #     return ModelInstance(self)

    # def __exit__(self, exc_type, exc_value, exc_traceback):
    #     self.finalise()

    def instance(self):
        return ModelInstance(self)

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
        results = {"properties": self.properties.raw_symbols}
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


class ModelInstance:

    def __init__(self, model: Model):
        # self.material = None  # TODO: material handler
        self.parameters = model.parameters.create_handler()
        self.properties = model.properties  # TODO: need to enable multiple instances
        # self.hierarchy = None  # TODO: model.hierarchy.create_handler()
        # self.residuals = None  # TODO: residual handler

        # self.numerics = NumericHandler(self)

        self.__function = model.function  # by reference

    def finalise(self, top_level: bool = True):
        """recreate the casadi nodes using the previously defined funtion.
        Use the parameters, child module properties, material objects, and
        locally defined non-canonical states to define a function,
        based on the :meth:`define` method's implementation.
        """
        # finalise child modules
        # for child in self.hierarchy.values():
        #     child.finalise(top_level=False)

        self.parameters.finalise()  # all required parameters connected?
        self.properties.check()  # all promised properties calculated?

        args = self.__collect_args()
        res = self.__function(args, squeeze_results=False)
        self.properties.finalise(res["properties"])
        # if top_level:
        #     self.numerics.prepare()
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
