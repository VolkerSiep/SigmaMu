from simu import AModel

from simu.examples.steam_system.thermo import (
    hp_condensate, hp_steam, lp_condensate, lp_steam, h2o_spec)


class VLE(AModel):
    r"""Constrain the given steam (:math:`s`) and condensate (:math:`c`) phases
    to be at thermodynamic equilibrium. 3 residuals are defined as follows:

    .. math::

        T_s = T_c \quad p_s = p_c \quad\text{and}\quad
        \mu_{s,\mathrm{H_2O}} = \mu_{c,\mathrm{H_2O}}
    """
    def interface(self):
        self.md("steam", h2o_spec)
        self.md("condensate", h2o_spec)

    def define(self):
        s, c = self.m["steam"], self.m["condensate"]
        self.ra("T_eq", c["T"] - s["T"], "K")
        self.ra("p_eq", c["p"] - s["p"], "bar")
        self.ra("mu_eq", c["mu"]["H2O"] - s["mu"]["H2O"], "kJ/mol")


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

        # self.pad("duty", 8, "MW")
        self.pad("vap_frac", 12, "%")
        self.prd("duty", "MW")


    def define(self):
        f = self.m["feed"]
        c = self.m["condensate"]
        s = self.m["steam"]

        # pressure from feed
        self.ra("p_f", f["p"] - s["p"], "bar")

        # material balance
        self.ra("N_bal", c["N"] + s["N"] - f["N"] , "kmol/h")

        # thermodynamic equilibrium
        with self.ha("vle", VLE) as unit:
            unit.materials.connect("steam", s)
            unit.materials.connect("condensate", c)

        # given duty and evaporation
        # self.ra("duty", self.pa["duty"] - (c["H"] + s["H"] - f["H"]) , "MW")
        self.pr["duty"] = c["H"] + s["H"] - f["H"]

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
        self.ra("p_c4", p_spec - bd["p"], "bar")

        # drum equilibrium
        self.ra("t_bd", bd["T"] - s["T"], "K")
        with self.ha("vle", VLE) as unit:
            unit.materials.connect("steam", s)
            unit.materials.connect("condensate", bf)

        # blowdown
        self.ra("blow_down", self.pa["blow_down"] * bfw["N"] - bd["N"], "mol/s")

        # drum balances
        res = bfw["N"] + bc["N"] + bs["N"] - bf["N"] - bd["N"] - s["N"]
        self.ra("n_balance", res, "mol/s")
        res = bfw["H"] + bc["H"] + bs["H"] - bf["H"] - bd["H"] - s["H"]
        self.ra("h_balance", res, "MW")


class SuperHeater(AModel):
    def interface(self):
        self.md("inlet", h2o_spec)
        self.md("outlet", h2o_spec)
        # self.pad("duty", 2, "MW")
        self.prd("duty", "MW")

    def define(self):
        i, o = self.m["inlet"], self.m["outlet"]
        self.ra("p", i["p"] - o["p"], "bar") # pressure
        self.ra("N_bal", o["N"]  - i["N"] , "kmol/h") # material
        # self.ra("duty", o["H"] - self.pa["duty"] - i["H"] , "MW") # enthalpy
        self.pr["duty"] = o["H"] - i["H"]

class Turbine(AModel):
    def interface(self):
        self.md("inlet", h2o_spec)
        self.md("outlet_steam", h2o_spec)
        self.md("outlet_cond", h2o_spec)
        self.pad("isentropic_efficiency", 80, "%")
        self.prd("duty", "MW")
        self.prd("vapour_fraction", "%")

    def define(self):
        f = self.m["inlet"]
        o_s = self.m["outlet_steam"]
        o_c = self.m["outlet_cond"]
        i_c = self.mcf("iso_cond", lp_condensate)
        i_s = self.mcf("iso_steam", lp_steam)

        # phase equilibria
        with self.ha("vle_o", VLE) as unit:
            unit.materials.connect("steam", o_s)
            unit.materials.connect("condensate", o_c)
        with self.ha("vle_i", VLE) as unit:
            unit.materials.connect("steam", i_s)
            unit.materials.connect("condensate", i_c)

        # conservations
        # p_out = self.pa["discharge_pressure"]
        # self.ra("outlet p", p_out - o_s["p"], "bar")
        self.ra("isentropic p", o_s["p"] - i_s["p"], "bar")
        self.ra("outlet n", f["N"] - o_s["N"] - o_c["N"], "kmol/h")
        self.ra("isentropic n", f["N"] - i_s["N"] - i_c["N"], "kmol/h")
        self.ra("isentropic S", f["S"] - i_s["S"] - i_c["S"], "kW/K")

        # duty calculation
        eta = self.pa["isentropic_efficiency"]
        self.pr["duty"] = eta * (f["H"] - i_c["H"] - i_s["H"])
        self.ra("duty", f["H"] - self.pr["duty"] - o_c["H"] - o_s["H"], "MW")
        self.pr["vapour_fraction"] = o_s["N"] / f["N"]


class TotalCondenser(AModel):
    def interface(self):
        self.md("inlet_steam", h2o_spec)
        self.md("inlet_cond", h2o_spec)
        self.md("outlet", h2o_spec)
        self.pad("temperature", 50, "degC")
        self.prd("duty", "MW")

    def define(self):
        i_s, i_c = self.m["inlet_steam"], self.m["inlet_cond"]
        o = self.m["outlet"]

        self.ra("n-balance", i_s["N"] + i_c["N"] - o["N"], "kmol/h")
        self.ra("pressure", i_s["p"] - o["p"], "bar")
        self.ra("temperature change", i_s["T"] - o["T"], "K")
        self.ra("temperature", o["T"] - self.pa["temperature"], "K")

        self.pr["duty"] = i_s["H"] + i_c["H"] - o["H"]

class SteamGeneration(AModel):
    def interface(self):
        self.pad("T_bfw", 250, "degC")
        self.pad("blow_down", 1, "%")
        self.pad("total_duty", 10, "MW")
        self.pad("vfrac_turbine", 90, "%")

    def define(self):
        c1 = self.mcf("c1", hp_condensate)
        c2 = self.mcf("c2", hp_condensate)
        c3 = self.mcf("c3", hp_condensate)
        s3 = self.mcf("s3", hp_steam)
        c4 = self.mcf("c4", hp_condensate)
        s5 = self.mcf("s5", hp_steam)
        s6 = self.mcf("s6", hp_steam)
        s7 = self.mcf("s7", lp_steam)
        c7 = self.mcf("c7", lp_condensate)
        c8 = self.mcf("c8", lp_condensate)

        with self.ha("boiler", Boiler) as unit:
            unit.materials.connect("feed", c2)
            unit.materials.connect("condensate", c3)
            unit.materials.connect("steam", s3)
        q_boiler = unit.properties["duty"]

        with self.ha("drum", SteamDrum) as unit:
            unit.materials.connect("bfw", c1)
            unit.materials.connect("blow_down", c4)
            unit.materials.connect("steam", s5)
            unit.materials.connect("boiler_feed", c2)
            unit.materials.connect("boiler_return_condensate", c3)
            unit.materials.connect("boiler_return_steam", s3)

        with self.ha("superheater", SuperHeater) as unit:
            unit.materials.connect("inlet", s5)
            unit.materials.connect("outlet", s6)
        q_superheater = unit.properties["duty"]

        with self.ha("turbine", Turbine) as unit:
            unit.materials.connect("inlet", s6)
            unit.materials.connect("outlet_steam", s7)
            unit.materials.connect("outlet_cond", c7)
        vfrac_turbine = unit.properties["vapour_fraction"]

        with self.ha("condenser", TotalCondenser) as unit:
            unit.materials.connect("inlet_steam", s7)
            unit.materials.connect("inlet_cond", c7)
            unit.materials.connect("outlet", c8)


        # BFW temperature and pressure (the latter from drum)
        self.ra("T_bfw", self.pa["T_bfw"] - c1["T"], "K")
        self.ra("p_c1", s5["p"] - c1["p"], "bar")

        # distribute duty for given vapour fraction
        res = q_superheater + q_boiler - self.pa["total_duty"]
        self.ra("total_duty", res, "MW")
        self.ra("turbine_vfrac", vfrac_turbine - self.pa["vfrac_turbine"], "%")