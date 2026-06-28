from simu import AModel, R_GAS, exp

class MuFitModel(AModel):
    def __init__(self, ref_material, fit_material):
        self.ref_material = ref_material
        self.fit_material = fit_material
        super().__init__()

    def interface(self):
        self.pad("T", 25, "degC")
        self.pad("p", 1, "bar")
        self.pad("N", 1, "mol")
        self.prd("q_mu", "dimensionless")
        self.prd("q_h", "dimensionless")
        self.prd("f_mu", "dimensionless")
        self.prd("q_rho", "dimensionless")

    def define(self):
        ref = self.mcs("ref", self.ref_material)
        fit = self.mcs("fit", self.fit_material)

        t, p, n = (self.pa[s] for s in "TpN")

        for s in ("ref", "fit"):
            self.ra(f"T_{s}", self.m[s]["T"] - t, "K")
            self.ra(f"p_{s}", self.m[s]["p"] - p, "bar")
            self.ra(f"N_{s}", self.m[s]["N"] - n, "mol")

        self.pr["q_mu"] = (ref["mu"]["H2O"] - fit["mu"]["H2O"]) / (R_GAS * t)
        self.pr["f_mu"] = exp(self.pr["q_mu"])  # for evaluation only
        self.pr["q_h"] = (ref["H"] - fit["H"]) / (n * R_GAS * t)
        self.pr["q_rho"] = fit["rho"] / ref["rho"] - 1
