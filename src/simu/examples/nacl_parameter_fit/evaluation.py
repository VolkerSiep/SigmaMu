from pathlib import Path
from yaml import safe_load
from simu import ThermoFitEvaluator, NumericHandler
from model import PSatModel

THIS_DIRECTORY = Path(__file__).parent
THERMO_FIT_DEFINITION_FILE = THIS_DIRECTORY / "thermo_fit_definition.yml"
WASHBURN_H2O_FILE = THIS_DIRECTORY / "Washburn_1928_h2o.yml"


def load_definition():
    with THERMO_FIT_DEFINITION_FILE.open() as f:
        data = safe_load(f)
    with WASHBURN_H2O_FILE.open() as f:
        data["datasets"]["washburn_h2o"] = safe_load(f)
    return data

def evaluate(models, definition):
    evaluator = ThermoFitEvaluator(models)
    result = evaluator.solve(definition)
    print(result["pure_h2o"].results)

def main():
    definition = load_definition()
    models = {"p_sat": NumericHandler(PSatModel.top())}
    evaluate(models, definition)


if __name__ == '__main__':
    main()