from simu import AModel
from .thermo import materials


class PSatModel(AModel):
    def interface(self):
        self.pad("T", 25, "degC")
        self.pad("w", 5, "%")
        self.pad("p_meas", 1, "bar")
        self.pad("N", 1, "mol")

        self.prd("q", "")
        self.prd("p_calc", "bar")

    def define(self):
        # create materials
        liq = self.mcs("liquid", materials["liquid"])
        gas = self.mcs("gas", materials["gas"])

        # specification of the system
        pa = self.pa
        self.ra("T", liq["T"] - pa["T"], "K")
        self.ra("N_liq", liq["N"], pa["N"], "mol")
        self.ra("N_gas", gas["N"], pa["N"], "mol")
        self.ra("w", liq["m"]["Na+"] + liq["m"]["Cl-"] - pa["w"] * liq["M"], "g")

        # VLE
        self.ra("T_eq", liq["T"] - gas["T"], "K")
        self.ra("p_eq", liq["p"] - gas["p"], "bar")
        self.ra("h2o_eq", liq["mu"]["H2O"] - gas["mu"]["H2O"], "kJ/mol")

        self.pr["p_calc"] = liq["p"]  # for evaluation
        self.pr["q"] = liq["p"] / pa["p_meas"] - 1  # objective to minimize