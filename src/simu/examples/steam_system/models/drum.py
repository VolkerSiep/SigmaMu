"""For documentation see tutorial
*Modelling a small steam generation with a turbine*"""

from simu import AModel
from models.species_balance import SpeciesBalance
from models.vle import VLE


class SteamDrum(AModel):
    def interface(self):
        self.md("bfw")
        self.md("blow_down")
        self.md("boiler_feed")
        self.md("boiler_return_cond")
        self.md("boiler_return_steam")
        self.md("steam")

        self.pad("pressure", 100, "bar")
        self.pad("blow_down", 1, "%")


    def define(self):
        bfw = self.m["bfw"]
        bd = self.m["blow_down"]
        bf = self.m["boiler_feed"]
        bc = self.m["boiler_return_cond"]
        bs = self.m["boiler_return_steam"]
        s = self.m["steam"]

        # pressure and blowdown flow spec
        self.ra("p_steam", self.pa["pressure"] - s["p"], "bar")
        self.ra("blow_down", self.pa["blow_down"] * bfw["N"] - bd["N"], "mol/s")

        # species balances
        with self.ha("n_balance", SpeciesBalance, num_in=3, num_out=3) as unit:
            unit.mcm(in_01=bfw, in_02=bc, in_03=bs,
                     out_01=bd, out_02=bf, out_03=s)

        # enthalpy balance
        res = bfw["H"] + bc["H"] + bs["H"] - bf["H"] - bd["H"] - s["H"]
        self.ra("h_balance", res, "MW")

        # phase equilibria
        with self.ha("vle", VLE) as unit:
            unit.mcm(gas=s, liquid=bf)
        #   this is a hack only valid for pure components, apologies!!
        self.ra("t_bd", bd["T"] - s["T"], "K")
        self.ra("p_bd", bd["p"] - s["p"], "bar")
