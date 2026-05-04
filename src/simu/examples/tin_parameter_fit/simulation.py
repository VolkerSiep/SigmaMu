from simu import AModel, NumericHandler, SimulationSolver, NHKeys
from simu.examples.tin_parameter_fit.thermo import tin_definition

class TinTransition(AModel):
    """Calculate the transition temperature of Tin"""
    def interface(self) -> None:
        self.pad("n", 1, "mol")
        self.pad("p", 1, "bar")
        self.pad("T_measured", 12.3, "degC")
        self.prd("dT_norm", "dimensionless")

    def define(self) -> None:
        tin = self.mcs("tin", tin_definition)
        n, mu = tin["n"], tin["mu"]
        self.ra("n_a-Sn", n["a-Sn"] - self.pa["n"], "mol")
        self.ra("n_b-Sn", n["b-Sn"] - self.pa["n"], "mol")
        self.ra("p", tin["p"] - self.pa["p"], "bar")
        self.ra("equilibrium", mu["a-Sn"] - mu["b-Sn"], "kJ/mol")
        self.pr["dT_norm"] = tin["T"] / self.pa["T_measured"] - 1

def main():
    numeric = NumericHandler(TinTransition.top())
    solver = SimulationSolver(numeric)
    report = solver.solve()
    props = report.prop_func(report.final_state)
    temperature = props[NHKeys.THERMO_PROPS]["tin"]["T"].to("degC")
    deviation = props[NHKeys.MODEL_PROPS]["dT_norm"]
    print(f"Transition temperature: {temperature:.2fP~}")
    print(f"Deviation: {deviation:.4gP~}")


if __name__ == '__main__':
    main()
