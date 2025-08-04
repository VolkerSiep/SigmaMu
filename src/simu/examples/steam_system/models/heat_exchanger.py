"""For documentation see tutorial
*Modelling a small steam generation with a turbine*"""

from simu import AModel

from models.species_balance import SpeciesBalance
from models.vle import VLE


class Boiler(AModel):
    def interface(self):
        self.md("feed")
        self.md("steam_out")
        self.md("cond_out")

        self.pad("vap_frac", 12, "%")
        self.prd("duty", "MW")

    def define(self):
        f, c, s = self.m["feed"], self.m["cond_out"], self.m["steam_out"]
        self.ra("p_f", f["p"] - s["p"], "bar")  # zero pressure drop

        with self.ha("n_balance", SpeciesBalance, num_out=2) as unit:
            unit.mcm(in_01=f, out_01=s, out_02=c)

        with self.ha("vle", VLE) as unit:
            unit.mcm(gas=s, liquid=c)

        # vapour fraction from boiler
        self.ra("vap_frac", self.pa["vap_frac"] * f["N"] - s["N"], "kmol/h")
        self.pr["duty"] = c["H"] + s["H"] - f["H"]


class SuperHeater(AModel):
    def interface(self):
        self.md("inlet")
        self.md("outlet")
        self.prd("duty", "MW")

    def define(self):
        i, o = self.m["inlet"], self.m["outlet"]
        self.ra("p", i["p"] - o["p"], "bar")  # zero pressure drop

        with self.ha("n-balance", SpeciesBalance) as unit:
            unit.mcm(in_01=i, out_01=o)

        self.pr["duty"] = o["H"] - i["H"]


class TotalCondenser(AModel):
    def interface(self):
        self.md("steam_in")
        self.md("cond_in")
        self.md("outlet")
        self.pad("temperature", 50, "degC")
        self.prd("duty", "MW")

        self.pad("_trial_phase_size", 1, "mol/s")

    def define(self):
        i_s, i_c, o = self.m["steam_in"], self.m["cond_in"], self.m["outlet"]
        t = self.mcf("steam_trial", i_s.definition)

        with self.ha("n-balance", SpeciesBalance, num_in=2) as unit:
            unit.mcm(in_01=i_s, in_02=i_c, out_01=o)

        with self.ha("vle", VLE) as unit:
            unit.mcm(gas=t, liquid=o)
        self.ra("trial_n", t["N"] - self.pa["_trial_phase_size"], "mol/s")

        self.ra("pressure", i_s["p"] - o["p"], "bar")
        self.ra("temperature", o["T"] - self.pa["temperature"], "K")

        self.pr["duty"] = i_s["H"] + i_c["H"] - o["H"]