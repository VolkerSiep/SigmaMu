from abc import ABC, abstractmethod

from ..utilities import QFunction
from .property import PropertyHandler
from .parameter import ParameterHandler
from .numeric import NumericHandler


class Model(ABC):

    def __init__(self):
        self.material = None  # material handler
        self.parameters = ParameterHandler()
        self.properties = PropertyHandler()

        self.interface()

        self.hierarchy = None  # hierarchy handler
        self.residuals = None  # residual handler

        self.define()

        self.__function = self.__make_function()

        self.numerics = NumericHandler(self)

    def __enter__(self):
        # prepare an object that can connect the actual symbols and ports.
        return self

    def __exit__(self, exc_type, exc_value, exc_traceback):
        self.finalise()

    def finalise(self):
        """recreate the casadi nodes using the previously defined funtion."""
        self.parameters.check()  # all required parameters connected?
        self.properties.check()  # all promised properties calculated?

        args = {
            "parameters": self.parameters.symbols
        }  # now they are connected
        res = self.function(args, squeeze_results=False)
        self.properties.symbols = res["properties"]
        return self

    @property
    def function(self):
        return self.__function

    @abstractmethod
    def interface(self):
        ...  # todo: replace with document string

    @abstractmethod
    def define(self):
        ...  # todo: replace with document string

    def __make_function(self):
        name = self.__class__.__name__
        args = {"parameters": self.parameters.symbols}
        results = {"properties": self.properties.symbols}
        return QFunction(args, results, name)


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
