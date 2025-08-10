"""For documentation see tutorial
*Modelling a small steam generation with a turbine*"""

from simu import AModel

from thermo import hp_condensate, hp_steam, lp_condensate, lp_steam
from models.heat_exchanger import Boiler, SuperHeater, TotalCondenser
from models.drum import SteamDrum
from models.turbine import Turbine


class SteamGeneration(AModel):
    def interface(self):
        self.pad("T_bfw", 250, "degC")
        self.pad("blow_down", 1, "%")
        self.pad("total_duty", 10, "MW")

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
            unit.mcm(inlet=m["s6"], steam_out=m["s7"], cond_out=m["c7"])

        with self.ha("condenser", TotalCondenser) as unit:
            unit.mcm(steam_in=m["s7"], cond_in=m["c7"], outlet=m["c8"])

        # distribute duty for given vapour fraction
        res = q_superheater + q_boiler - self.pa["total_duty"]
        self.ra("total_duty", res, "MW")

        # BFW temperature and pressure (the latter from drum)
        self.ra("T_bfw", self.pa["T_bfw"] - m["c1"]["T"], "K")
        self.ra("p_c1", m["s5"]["p"] - m["c1"]["p"], "bar")
