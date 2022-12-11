class Model:
    def __init__(self):
        self.material = None  # material handler
        self.variables = None  # parameter handler
        self.residuals = None  # residual handler
        self.hierarchy = None  # hierarchy handler

        self.interface()
        self.define()

    def __enter__(self):
        pass
        # create and return an object that can connect the actual
        # symbols and ports.

    def __exit__(self, exc_type, exc_value, exc_traceback):
        # create the casadi nodes using the previously defined funtion.


class MyModel(Model):

    def interface(self):
        # require ports
        self.material.require("inlet", "NitricAcid")

        # define process parameter
        self.variables.define("pressure_drop", 0.3, "bar")  # with default
        self.variables.define("heat_loss", unit="kW")  # no default

        self.variables.provide("froth_density", unit="kg/m**3")

    def define(self):
        # need to define casadi graph that defines output from input
        # base class has made SX source nodes ready upfront.
        # How do I know UOM of input?
        #   - parameters are fine
        #   - material variables and child module variables:
        #     from their definitions (QFunction.res_units) ?

        # create own internal materials
        self.material.create("intern", "NOxGas")

        dp = self.variables["pressure_drop"]
        T = self.material["inlet"]["T"]

        self.variables["froth_density"] = dp * T  # not really

        dp_calc = self.material["s301"]["p"] - self.material["intern"]["p"]
        self.residuals.add("dp", dp_calc - dp, "5 bar")  # last arg is tolerance basis




