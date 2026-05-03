from pathlib import Path
from yaml import safe_load
from simu import ThermoFitSolver, NumericHandler

from thermo import thermo_source
from simulation import TinTransition

THERMO_FIT_DEFINITION_FILE = Path(__file__).parent / "thermo_fit_definition.yml"


def main():
    models = {"transition_model": NumericHandler(TinTransition.top())}
    with THERMO_FIT_DEFINITION_FILE.open() as f:
        definition = safe_load(f)

    solver = ThermoFitSolver(models, thermo_source)

    solver.solve(definition)


if __name__ == '__main__':
    main()