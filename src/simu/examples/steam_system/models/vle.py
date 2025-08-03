"""For documentation see tutorial
*Modelling a small steam generation with a turbine*"""

from simu import AModel


class VLE(AModel):
    def interface(self):
        self.md("gas")
        self.md("liquid")

    def define(self):
        g, l = self.m["gas"], self.m["liquid"]
        self.ra("T_eq", l["T"] - g["T"], "K")
        self.ra("p_eq", l["p"] - g["p"], "bar")

        for s in g["mu"].keys() & l["mu"].keys():
            self.ra(f"mu_eq_{s}", g["mu"][s] - l["mu"][s], "kJ/mol")
