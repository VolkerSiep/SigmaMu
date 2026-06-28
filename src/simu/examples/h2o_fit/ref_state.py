from thermo import lp_condensate, condensate_new
from model import MuFitModel
from simu import NumericHandler, SimulationSolver, NHKeys, R_GAS, Quantity, InitialState

T = Quantity(25, "degC")

def main():
    lp_condensate.initial_state = InitialState.from_cbar(25, 1.0, [1.0])
    model = MuFitModel(ref_material=condensate_new, fit_material=lp_condensate)
    numeric = NumericHandler(model.create_proxy().finalise())
    solver = SimulationSolver(numeric, output=None)
    result = solver.solve().properties[NHKeys.MODEL_PROPS]
    print("d_h  = ", result["q_h"] * (R_GAS * T))
    print("d_s = ", -result["q_mu"] * R_GAS)


if __name__ == '__main__':
    main()