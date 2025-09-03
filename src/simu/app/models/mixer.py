from numpy.ma.core import equal

from simu import AModel
from .basic import SpeciesBalance, EnthalpyBalance


class Mixer(AModel):
    def __init__(self, num_out=2, equal_pressure=True,
                 tol_unit_n="kmol/h", tol_unit_h="MW", tol_unit_p="bar"):
        super().__init__()
        self._num_out = num_out
        self._tol_unit_p = tol_unit_p
        self._tol_unit_n = tol_unit_n
        self._tol_unit_h = tol_unit_h
        self._equal_pressure = equal_pressure

    def interface(self):
        self.md("feed")
        for name in self.out_names():
            self.md(name)

    def define(self):
        m = self.m
        num_out = self._num_out
        connections = {n: m[n] for n in self.out_names()} | {"in_01": m["feed"]}

        with self.ha("n-balance", SpeciesBalance,
                     num_out=num_out, tol_unit=self._tol_unit_n) as unit:
            unit.mcm(**connections)
        with self.ha("H-balance", EnthalpyBalance,
                     num_out=num_out, tol_unit=self._tol_unit_h) as unit:
            unit.mcm(**connections)

        if self._equal_pressure:
            p_in = self.m["feed"]["p"]
            for n in self.out_names():
                self.ra("p{n[3:]}", m[n]["p"] - p_in, self._tol_unit_p)

    def out_names(self):
        for i in range(self._num_out):
            yield f"out_{i + 1:02d}"
