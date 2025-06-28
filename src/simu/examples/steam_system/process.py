from simu import AModel, MaterialSpec
from thermo import hp_condensate, hp_steam

h2o_spec = MaterialSpec(["H2O"])

class Boiler(AModel):
    """This boiler model assumes a condensate flow as feed, which is partially
    evaporated by applying the specified duty. No pressure drop is considered.

    As such, the unit is fully specified - however, the constraint to fix the
    obtained vaporization ratio is added, by that fixing the circulation feed
    flow to the boiler.
    """
    def interface(self):
        self.md("feed", h2o_spec)
        self.md("steam", h2o_spec)
        self.md("condensate", h2o_spec)

        self.pad("duty", 10, "MW")
        self.pad("vap_frac", 12, "%")

    def define(self):
        f = self.m["feed"]
        c = self.m["condensate"]
        s = self.m["steam"]

        # pressure from feed
        self.ra("p_f", f["p"] - s["p"], "bar")

        # material balance
        self.ra("N_bal", c["N"] + s["N"] - f["N"] , "kmol/h")

        # thermodynamic equilibrium
        self.ra("T_eq", c["T"] - s["T"], "K")
        self.ra("p_eq", c["p"] - s["p"], "bar")
        self.ra("mu_eq", c["mu"]["H2O"] - s["mu"]["H2O"], "kJ/mol")

        # given duty and evaporation
        self.ra("duty", self.pa["duty"] - (c["H"] + s["H"] - f["H"]) , "MW")

        # vapour fraction from boiler
        self.ra("vap_frac", self.pa["vap_frac"] * f["N"] - s["N"], "kmol/h")


class SteamDrum(AModel):
    def interface(self):
        self.pad("pressure", 100, "bar")
        self.pad("blow_down", 1, "%")

        self.md("bfw", h2o_spec)
        self.md("blow_down", h2o_spec)
        self.md("boiler_feed", h2o_spec)
        self.md("boiler_return_condensate", h2o_spec)
        self.md("boiler_return_steam", h2o_spec)
        self.md("steam", h2o_spec)


    def define(self):
        bfw = self.m["bfw"]
        bd = self.m["blow_down"]
        bf = self.m["boiler_feed"]
        bc = self.m["boiler_return_condensate"]
        bs = self.m["boiler_return_steam"]
        s = self.m["steam"]


        # drum pressure and equal pressure
        p_spec = self.pa["pressure"]
        self.ra("p_steam", p_spec - s["p"], "bar")
        self.ra("p_c2", p_spec - bf["p"], "bar")
        self.ra("p_c4", p_spec - bd["p"], "bar")

        # drum equilibrium
        self.ra("t_bf", bf["T"] - s["T"], "K")
        self.ra("t_bd", bd["T"] - s["T"], "K")
        self.ra("eq_drum", bf["mu"]["H2O"] - s["mu"]["H2O"], "kJ/mol")

        # blowdown
        self.ra("blow_down", self.pa["blow_down"] * bfw["N"] - bd["N"], "mol/s")

        # drum balances
        res = bfw["N"] + bc["N"] + bs["N"] - bf["N"] - bd["N"] - s["N"]
        self.ra("N_bal_drum", res, "mol/s")
        res = bfw["H"] + bc["H"] + bs["H"] - bf["H"] - bd["H"] - s["H"]
        self.ra("H_bal_drum", res, "MW")


class SteamGeneration(AModel):
    def interface(self):
        self.pad("T_bfw", 250, "degC")
        self.pad("blow_down", 1, "%")

    def define(self):
        c1 = self.mcf("c1", hp_condensate)
        c2 = self.mcf("c2", hp_condensate)
        c3 = self.mcf("c3", hp_condensate)
        s3 = self.mcf("s3", hp_steam)
        c4 = self.mcf("c4", hp_condensate)
        s5 = self.mcf("s5", hp_steam)

        with self.ha("boiler", Boiler) as boiler:
            boiler.materials.connect("feed", c2)
            boiler.materials.connect("condensate", c3)
            boiler.materials.connect("steam", s3)

        with self.ha("drum", SteamDrum) as drum:
            drum.materials.connect("bfw", c1)
            drum.materials.connect("blow_down", c4)
            drum.materials.connect("steam", s5)
            drum.materials.connect("boiler_feed", c2)
            drum.materials.connect("boiler_return_condensate", c3)
            drum.materials.connect("boiler_return_steam", s3)

        # BFW temperature and pressure (the latter from drum)
        self.ra("T_bfw", self.pa["T_bfw"] - c1["T"], "K")
        self.ra("p_c1", s5["p"] - c1["p"], "bar")

# TODO: when collecting state vector, do not collect connected states!




