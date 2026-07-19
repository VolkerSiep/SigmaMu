from pathlib import Path
from yaml import safe_load
from simu import NumericHandler
from model import PSatModel

THIS_DIRECTORY = Path(__file__).parent
THERMO_FIT_DEFINITION_FILE = THIS_DIRECTORY / "thermo_fit_definition.yml"
WASHBURN_H2O_FILE = THIS_DIRECTORY / "Washburn_1928_h2o.yml"
WASHBURN_ALL_FILE = THIS_DIRECTORY / "Washburn_1928_all.yml"


def load_definition():
    with THERMO_FIT_DEFINITION_FILE.open() as f:
        data = safe_load(f)
    with WASHBURN_H2O_FILE.open() as f:
        data["datasets"]["washburn_h2o"] = safe_load(f)
    with WASHBURN_ALL_FILE.open() as f:
        data["datasets"]["washburn_all"] = safe_load(f)
    return data


def define_models():
    return {"p_sat": NumericHandler(PSatModel.top())}