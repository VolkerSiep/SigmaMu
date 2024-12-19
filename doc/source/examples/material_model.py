from pprint import pprint
from simu import Model, NumericHandler, qsum
from ideal_gas_material import ch4_ideal

class Source(Model):
    """A model of methane source"""

    def interface(self):
        self.parameters.define("T", 25, "degC")
        self.parameters.define("p", 1, "bar")
        self.parameters.define("N", 1000, "kmol/hr")

    def define(self):
        flow = self.materials.create_flow("source", ch4_ideal)
        flow["N"] = qsum(flow["n"])
        self.residuals.add( "T", self.parameters["T"] - flow["T"], "K")
        self.residuals.add( "p", self.parameters["p"] - flow["p"], "bar")
        self.residuals.add( "N", self.parameters["N"] - flow["N"], "kmol/hr")

numeric = NumericHandler(Source.top())
args = numeric.arguments
pprint(args)
func = numeric.function
# pprint(func(args))

# pprint(func.arg_structure)



