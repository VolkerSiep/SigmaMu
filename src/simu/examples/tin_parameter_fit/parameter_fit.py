from pathlib import Path
from yaml import safe_load
from simu import ThermoFitSolver, NumericHandler

from simu.examples.tin_parameter_fit.thermo import thermo_source
from simu.examples.tin_parameter_fit.simulation import TinTransition

THIS_DIRECTORY = Path(__file__).parent
THERMO_FIT_DEFINITION_FILE = THIS_DIRECTORY / "thermo_fit_definition.yml"

def load_definition():
    with THERMO_FIT_DEFINITION_FILE.open() as f:
        return safe_load(f)

def main():
    models = {"transition_model": NumericHandler(TinTransition.top())}
    solver = ThermoFitSolver(models, thermo_source, epsilon_q=1e-7)
    result = solver.solve(load_definition())
    print(result.final_parameters)


if __name__ == '__main__':
    main()