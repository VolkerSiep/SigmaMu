"""For documentation see tutorial
*Modelling a small steam generation with a turbine*"""

from simu import AModel


class SpeciesBalance(AModel):
    def __init__(self, num_in: int = 1, num_out: int = 1, tol_unit = "kmol/h"):
        self._num_in, self._num_out = num_in, num_out
        self._tol_unit = tol_unit
        super().__init__()

    def interface(self):
        for name in self.in_names():
            self.md(name)
        for name in self.out_names():
            self.md(name)

    def define(self):
        inlets = [self.m[name] for name in self.in_names()]
        outlets = [self.m[name] for name in self.out_names()]

        d_n = sum(m["n"] for m in inlets) -  sum(m["n"] for m in outlets)
        for species, dn_i in d_n.items():
            self.ra(species, dn_i, self._tol_unit)

    def in_names(self):
        for i in range(self._num_in):
            yield f"in_{i + 1:02d}"

    def out_names(self):
        for i in range(self._num_out):
            yield f"out_{i + 1:02d}"
