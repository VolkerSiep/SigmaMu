from simu import Model
from .material_factory import ch4_ideal_gas

class Source(Model):
    """A model of a methane source"""

    def interface(self):
        self.parameters.define("T", 25, "degC")
        self.parameters.define("p", 1, "bar")
        self.parameters.define("V", 10, "m^3/hr")

    def define(self):
        src = self.materials.create_flow("source", ch4_ideal_gas)
        self.residuals.add("T", self.parameters["T"] - src["T"], "K")
        self.residuals.add("p", self.parameters["p"] - src["p"], "bar")
        self.residuals.add("V", self.parameters["V"] - src["V"], "m^3/h")




