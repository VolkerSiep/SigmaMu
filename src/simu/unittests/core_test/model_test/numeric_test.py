from simu.core.model import NumericHandler
from .models import *


def test_parameters():
    proxy = SimpleParameterTestModel.top()
    numeric = NumericHandler(proxy)
    args = numeric.function.arg_structure
    assert args['parameters/length'] == 'm'


def test_properties():
    proxy = PropertyTestModel.top()
    numeric = NumericHandler(proxy)
    results = numeric.function.result_structure
    assert results['model_props/area'] == 'm ** 2'


def test_residuals():
    proxy = ResidualTestModel.top()
    numeric = NumericHandler(proxy)
    results = numeric.function.result_structure
    assert results['residuals/area'] == "m ** 2"


def test_material_collect_states():
    args, _ = create_material_functions()
    for k in range(4):
        assert f"states/local/x_{k:03d}" in args


def test_material_collect_props():
    _, results = create_material_functions()
    print(results)
    # TODO: make assert statement


def test_material_collect_thermo_param():
    args, _ = create_material_functions()
    print(args)
    # TODO: make assert statement


def create_material_functions():
    proxy = MaterialTestModel3.top()
    numeric = NumericHandler(proxy)
    args = numeric.function.arg_structure
    results = numeric.function.result_structure
    return args, results
