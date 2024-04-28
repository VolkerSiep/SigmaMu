from simu.core.model import NumericHandler
from simu.core.utilities import assert_reproduction
from .models import *


def test_parameters():
    proxy = SimpleParameterTestModel.top()
    numeric = NumericHandler(proxy)
    args = numeric.function.arg_structure
    assert args['model_params/length'] == 'm'


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
    start = "thermo_props/local/"
    results = [res[len(start):] for res in results if res.startswith(start)]
    assert_reproduction(results)


def test_material_collect_thermo_param():
    args, _ = create_material_functions()
    start = "thermo_params/default/"
    args = [arg[len(start):] for arg in args if arg.startswith(start)]
    assert_reproduction(args)


def test_hierarchy_collect_numerics():
    numeric = NumericHandler(HierarchyTestModel2.top())
    results = numeric.function.result_structure
    assert "model_props/square/area" in results


def test_suare_model():
    numeric = NumericHandler(SquareTestModel.top())



def create_material_functions():
    """Make a function out of a model defining materials"""
    proxy = MaterialTestModel3.top()
    numeric = NumericHandler(proxy)
    args = numeric.function.arg_structure
    results = numeric.function.result_structure
    return args, results
