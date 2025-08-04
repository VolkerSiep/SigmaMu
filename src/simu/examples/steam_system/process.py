from simu import AModel

from thermo import hp_condensate, hp_steam, lp_condensate, lp_steam, h2o_spec
from models.heat_exchanger import Boiler, SuperHeater, TotalCondenser
from models.species_balance import SpeciesBalance
from models.vle import VLE


class SteamDrum(AModel):
    def interface(self):
        self.pad("pressure", 100, "bar")
        self.pad("blow_down", 1, "%")

        self.md("bfw", h2o_spec)
        self.md("blow_down", h2o_spec)
        self.md("boiler_feed", h2o_spec)
        self.md("boiler_return_cond", h2o_spec)
        self.md("boiler_return_steam", h2o_spec)
        self.md("steam", h2o_spec)


    def define(self):
        bfw = self.m["bfw"]
        bd = self.m["blow_down"]
        bf = self.m["boiler_feed"]
        bc = self.m["boiler_return_cond"]
        bs = self.m["boiler_return_steam"]
        s = self.m["steam"]

        # drum pressure and equal pressure
        p_spec = self.pa["pressure"]
        self.ra("p_steam", p_spec - s["p"], "bar")
        self.ra("p_c4", p_spec - bd["p"], "bar")

        # drum equilibrium
        self.ra("t_bd", bd["T"] - s["T"], "K")
        with self.ha("vle", VLE) as unit:
            unit.mcm(gas=s, liquid=bf)

        # blowdown
        self.ra("blow_down", self.pa["blow_down"] * bfw["N"] - bd["N"], "mol/s")

        # drum balances
        with self.ha("n_balance", SpeciesBalance, num_in=3, num_out=3) as unit:
            unit.mcm(in_01=bfw, in_02=bc, in_03=bs,
                     out_01=bd, out_02=bf, out_03=s)

        res = bfw["H"] + bc["H"] + bs["H"] - bf["H"] - bd["H"] - s["H"]
        self.ra("h_balance", res, "MW")


class Turbine(AModel):
    def interface(self):
        self.md("inlet")
        self.md("outlet_steam")
        self.md("outlet_cond")
        self.pad("isentropic_efficiency", 80, "%")
        self.prd("duty", "MW")
        self.prd("vapour_fraction", "%")

    def define(self):
        f = self.m["inlet"]
        o_c = self.m["outlet_cond"]
        o_s = self.m["outlet_steam"]
        i_c = self.mcf("iso_cond", o_c.definition)
        i_s = self.mcf("iso_steam", o_s.definition)

        with self.ha("isentropic", CompressionCondensation) as unit:
            unit.mcm(feed=f, out_steam=i_s, out_cond=i_c)
        iso_duty = unit.pr["duty"]
        with self.ha("actual", CompressionCondensation) as unit:
            unit.mcm(feed=f, out_steam=o_s, out_cond=o_c)
        duty = unit.pr["duty"]

        # conservations
        self.ra("discharge_pressure", o_s["p"] - i_s["p"], "bar")
        self.ra("entropy_balance", f["S"] - i_s["S"] - i_c["S"], "kW/K")

        # efficiency specification
        eta = self.pa["isentropic_efficiency"]
        self.ra("duty", duty - eta * iso_duty, "MW")

        self.pr["duty"] = duty  # eta * iso_duty
        self.pr["vapour_fraction"] = o_s["N"] / f["N"]


class CompressionCondensation(AModel):
    def interface(self):
        self.md("feed")
        self.md("out_steam")
        self.md("out_cond")

        self.prd("duty", "MW")

    def define(self):
        f = self.m["feed"]
        s = self.m["out_steam"]
        c = self.m["out_cond"]

        with self.ha("vle", VLE) as unit:
            unit.mcm(gas=s, liquid=c)

        with self.ha("n-balance", SpeciesBalance, num_out=2) as unit:
            unit.mcm(in_01=f, out_01=s, out_02=c)

        self.pr["duty"] = f["H"] - c["H"] - s["H"]


class SteamGeneration(AModel):
    def interface(self):
        self.pad("T_bfw", 250, "degC")
        self.pad("blow_down", 1, "%")
        self.pad("total_duty", 10, "MW")
        self.pad("vfrac_turbine", 90, "%")

    def define(self):
        # [(definition, prefix, [<stream numbers>], ...]
        streams = [(hp_steam, "s", [3, 5, 6]),
                   (hp_condensate, "c", [1, 2, 3, 4]),
                   (lp_steam, "s", [7]),
                   (lp_condensate, "c", [7, 8])]

        for definition, prefix, numbers in streams:
            for i in numbers:
                self.mcf(f"{prefix}{i}", definition)
        m = self.materials

        with self.ha("boiler", Boiler) as unit:
            unit.mcm(feed=m["c2"], cond_out=m["c3"], steam_out=m["s3"])
        q_boiler = unit.pr["duty"]

        with self.ha("drum", SteamDrum) as unit:
            unit.mcm(bfw=m["c1"], blow_down=m["c4"], steam=m["s5"],
                     boiler_feed=m["c2"], boiler_return_cond=m["c3"],
                     boiler_return_steam=m["s3"])

        with self.ha("superheater", SuperHeater) as unit:
            unit.mcm(inlet=m["s5"], outlet=m["s6"])
        q_superheater = unit.pr["duty"]

        with self.ha("turbine", Turbine) as unit:
            unit.mcm(inlet=m["s6"], outlet_steam=m["s7"], outlet_cond=m["c7"])
        vfrac_turbine = unit.pr["vapour_fraction"]

        with self.ha("condenser", TotalCondenser) as unit:
            unit.mcm(steam_in=m["s7"], cond_in=m["c7"], outlet=m["c8"])

        # BFW temperature and pressure (the latter from drum)
        self.ra("T_bfw", self.pa["T_bfw"] - m["c1"]["T"], "K")
        self.ra("p_c1", m["s5"]["p"] - m["c1"]["p"], "bar")

        # distribute duty for given vapour fraction
        res = q_superheater + q_boiler - self.pa["total_duty"]
        self.ra("total_duty", res, "MW")
        self.ra("turbine_vfrac", vfrac_turbine - self.pa["vfrac_turbine"], "%")