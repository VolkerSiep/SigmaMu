from pathlib import Path
from yaml import safe_load
from simu import NumericHandler, ThermoFitSolver, ThermoFitEvaluator

from simu.examples.tin_parameter_fit.thermo import thermo_source
from simu.examples.tin_parameter_fit.simulation import TinTransition

THIS_DIRECTORY = Path(__file__).parent
THERMO_FIT_DEFINITION_FILE = THIS_DIRECTORY / "thermo_fit_definition.yml"

def load_definition():
    with THERMO_FIT_DEFINITION_FILE.open() as f:
        return safe_load(f)

def fit(models, definition):
    solver = ThermoFitSolver(models, thermo_source, epsilon_q=1e-7)
    result = solver.solve(definition)
    print(result.final_parameters)
    return result.final_parameters

def evaluate(models, definition, parameters):
    evaluator = ThermoFitEvaluator(models)
    print("Original:")
    result = evaluator.solve(definition)
    print(result["by_temp"].results)
    print("\nAfter data fit:")
    result = evaluator.solve(definition, parameters)
    print(result["by_temp"].results)

def main():
    definition = load_definition()
    models = {"transition_model": NumericHandler(TinTransition.top())}
    parameters = fit(models, definition)
    evaluate(models, definition, parameters)


if __name__ == '__main__':
    main()