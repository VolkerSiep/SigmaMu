from simu import AModel
from .material_factory import ch4_ideal_gas

class Source(AModel):
    """A model of a methane source"""

    def interface(self):
        self.pad("T", 25, "degC")
        self.pad("p", 1, "bar")
        self.pad("V", 10, "m^3/hr")

    def define(self):
        src = self.mcf("source", ch4_ideal_gas)
        self.ra("T", self.pa["T"] - src["T"], "K")
        self.ra("p", self.pa["p"] - src["p"], "bar")
        self.ra("V", self.pa["V"] - src["V"], "m^3/h")




