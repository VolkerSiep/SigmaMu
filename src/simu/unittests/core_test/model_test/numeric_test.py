from simu.core.model import NumericHandler
from simu.core.utilities import assert_reproduction
from simu.core.thermo.parameters import StringDictThermoSource
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


def test_square_model():
    numeric = NumericHandler(SquareTestModel.top())
    ref = {"args": numeric.function.arg_structure,
           "res": numeric.function.result_structure}
    assert_reproduction(ref)


def test_square_model_args():
    model = SquareTestModel()
    material = model.no2sol
    numeric = NumericHandler(model.create_proxy().finalise())

    data = {
        'H0S0ReferenceState': {
            's_0': {'CH3-(CH2)2-CH3': '0 J/(mol*K)',
                    'CH3-CH2-CH3': '0 J/(mol*K)'},
            'dh_form': {'CH3-(CH2)2-CH3': '0 kJ/mol',
                        'CH3-CH2-CH3': '0 kJ/mol'},
            'T_ref': '25 degC',
            'p_ref': '1 atm'},
        'LinearHeatCapacity': {
            'cp_a': {'CH3-(CH2)2-CH3': '98 J/(mol*K)',
                     'CH3-CH2-CH3': '75 J/(mol*K)'},
            'cp_b': {'CH3-(CH2)2-CH3': '0 J/(mol*K*K)',
                     'CH3-CH2-CH3': '0 J/(mol*K*K)'}},
        'CriticalParameters': {
            'T_c': {'CH3-(CH2)2-CH3': '425 K', 'CH3-CH2-CH3': '370 K'},
            'p_c': {'CH3-(CH2)2-CH3': '38 bar', 'CH3-CH2-CH3': '42.5 bar'}},
        'MixingRule_A': {'T_ref': '25 degC'},
        'VolumeShift': {
            'c_i': {'CH3-(CH2)2-CH3': '0 m ** 3 / mol',
                    'CH3-CH2-CH3': '0 m ** 3 / mol'}}
    }
    material.store.add_source("default", StringDictThermoSource(data))

    args = numeric.arguments
    print(args)


def create_material_functions():
    """Make a function out of a model defining materials"""
    proxy = MaterialTestModel3.top()
    numeric = NumericHandler(proxy)
    args = numeric.function.arg_structure
    results = numeric.function.result_structure
    return args, results
