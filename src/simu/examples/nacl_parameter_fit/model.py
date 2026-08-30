from simu import AModel
from simu.app.models.basic import PhaseEquilibrium
from thermo import materials


class PSatModel(AModel):
    def interface(self):
        self.pad("T", 100, "degC")
        self.pad("w", 0.0001, "%")
        self.pad("p_meas", 1, "bar")
        self.pad("N", 1, "mol")  # arbitrary phase size

        self.prd("q", "dimensionless")
        self.prd("p_calc", "bar")

    def define(self):
        # create materials
        liq = self.mcs("liquid", materials["liquid"])
        gas = self.mcs("gas", materials["gas"])

        # specification of the system
        pa = self.pa
        self.ra("T", liq["T"] - pa["T"], "K")
        self.ra("N_liq", liq["N"] - pa["N"], "mol")
        self.ra("N_gas", gas["N"] - pa["N"], "mol")
        self.ra("w", liq["m"]["Na+"] + liq["m"]["Cl-"] - pa["w"] * liq["M"], "g")

        # VLE (forces equal T, p and mu_h2o for both phases)
        with self.ha("vle", PhaseEquilibrium, flow=False) as vle:
            vle.mcm(phase_1=liq, phase_2=gas)

        self.pr["p_calc"] = liq["p"]  # for evaluation
        self.pr["q"] = liq["p"] / pa["p_meas"] - 1  # objective to minimize