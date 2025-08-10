"""For documentation see tutorial
*Modelling a small steam generation with a turbine*"""

from simu import AModel
from models.vle import VLE
from models.species_balance import SpeciesBalance


class PartialCondensation(AModel):
    def interface(self):
        self.md("feed")
        self.md("steam_out")
        self.md("cond_out")
        self.prd("duty", "MW")
        self.prd("vapour_fraction", "%")

    def define(self):
        f = self.m["feed"]
        s = self.m["steam_out"]
        c = self.m["cond_out"]

        with self.ha("vle", VLE) as unit:
            unit.mcm(gas=s, liquid=c)
        with self.ha("n-balance", SpeciesBalance, num_out=2) as unit:
            unit.mcm(in_01=f, out_01=s, out_02=c)

        self.pr["duty"] = f["H"] - c["H"] - s["H"]
        self.pr["vapour_fraction"] = s["N"] / f["N"]


class Turbine(AModel):
    def interface(self):
        self.md("inlet")
        self.md("steam_out")
        self.md("cond_out")
        self.pad("isentropic_efficiency", 80, "%")
        self.pad("vapour_fraction", 90, "%")
        self.prd("duty", "MW")

    def define(self):
        f = self.m["inlet"]
        o_c = self.m["cond_out"]
        o_s = self.m["steam_out"]
        i_c = self.mcf("cond_iso", o_c.definition)
        i_s = self.mcf("steam_iso", o_s.definition)

        with self.ha("isentropic", PartialCondensation) as unit:
            unit.mcm(feed=f, steam_out=i_s, cond_out=i_c)
        iso_duty = unit.pr["duty"]
        with self.ha("actual", PartialCondensation) as unit:
            unit.mcm(feed=f, steam_out=o_s, cond_out=o_c)
        duty, v_frac = unit.pr["duty"], unit.pr["vapour_fraction"]

        # set vapour fraction
        self.ra("vapour_frac", v_frac - self.pa["vapour_fraction"], "%")

        # constraints on expansion
        self.ra("discharge_pressure", o_s["p"] - i_s["p"], "bar")
        self.ra("entropy_balance", f["S"] - i_s["S"] - i_c["S"], "kW/K")

        # duties and efficiencies
        eta = self.pa["isentropic_efficiency"]
        self.ra("duty", duty - eta * iso_duty, "MW")
        self.pr["duty"] = duty
