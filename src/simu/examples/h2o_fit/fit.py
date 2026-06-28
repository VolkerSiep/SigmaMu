from pathlib import Path
from yaml import safe_dump

from common import define_models, load_definition
from simu import ThermoFitSolver, quantity_dict_to_strings
from thermo import thermo_store

PARAM_FILE = Path(__file__).parent / "parameters_fit.yml"

def main():
    definition = load_definition()
    models = define_models()
    solver = ThermoFitSolver(models, thermo_store.get_source("IAPWS_fit"))
    result = solver.solve(definition)
    param = quantity_dict_to_strings(result.final_parameters)

    with PARAM_FILE.open("w") as file:
        safe_dump(param, file)


if __name__ == '__main__':
    main()