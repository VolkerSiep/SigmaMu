from simu.core.model import NumericHandler
from simu.core.utilities import assert_reproduction, flatten_dictionary
from simu.core.thermo.parameters import StringDictThermoSource
from .models import *

DATA = {
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
        'p_c': {'CH3-(CH2)2-CH3': '38 bar', 'CH3-CH2-CH3': '42.5 bar'},
        'omega': {'CH3-CH2-CH3': 0.199, 'CH3-(CH2)2-CH3': 0.153}},
    'MixingRule_A': {'T_ref': '25 degC'},
    'VolumeShift': {
        'c_i': {'CH3-(CH2)2-CH3': '0 m ** 3 / mol',
                'CH3-CH2-CH3': '0 m ** 3 / mol'}},
    'BostonMathiasAlphaFunction': {
        'eta': {'CH3-CH2-CH3': 0, 'CH3-(CH2)2-CH3': 0}}
}


def test_parameters():
    proxy = SimpleParameterTestModel.top()
    numeric = NumericHandler(proxy)
    args = numeric.function.arg_structure
    assert args['model_params']['length'] == 'm'


def test_properties():
    proxy = PropertyTestModel.top()
    numeric = NumericHandler(proxy)
    results = numeric.function.result_structure
    assert results['model_props']['area'] == 'm ** 2'


def test_residuals():
    proxy = ResidualTestModel.top()
    numeric = NumericHandler(proxy)
    results = numeric.function.result_structure
    assert results['residuals']['area'] == "m ** 2"


def test_material_collect_states():
    args, _ = create_material_functions()
    for k in range(4):
        assert f"x_{k:03d}" in args["states"]["local"]


def test_material_collect_props():
    _, results = create_material_functions()
    assert_reproduction(results["thermo_props"]["local"])


def test_material_collect_thermo_param():
    args, _ = create_material_functions()
    assert_reproduction(args["thermo_params"]["default"])


def test_hierarchy_collect_numerics():
    numeric = NumericHandler(HierarchyTestModel2.top())
    results = numeric.function.result_structure
    assert "area" in results["model_props"]["square"]


def test_square_model():
    numeric = NumericHandler(SquareTestModel.top())
    ref = {"args": numeric.function.arg_structure,
           "res": numeric.function.result_structure}
    assert_reproduction(ref)


def test_square_model_args():
    model = SquareTestModel()
    material = model.no2sol
    numeric = NumericHandler(model.create_proxy().finalise())
    material.store.add_source("default", StringDictThermoSource(DATA))
    struct = numeric.function.arg_structure
    args = numeric.arguments
    check_same_keys(struct, args)


def test_square_model_call():
    model = SquareTestModel()
    material = model.no2sol
    numeric = NumericHandler(model.create_proxy().finalise())
    material.store.add_source("default", StringDictThermoSource(DATA))
    args = numeric.arguments
    res = flatten_dictionary(numeric.function(args))
    res = {k: f"{v:.6f~}" for k, v in res.items()}
    assert_reproduction(res)


def create_material_functions():
    """Make a function out of a model defining materials"""
    proxy = MaterialTestModel3.top()
    numeric = NumericHandler(proxy)
    args = numeric.function.arg_structure
    results = numeric.function.result_structure
    return args, results


def test_collect_hierarchy_material():
    proxy = MaterialParentTestModel.top()
    for port_props in (True, False):
        numeric = NumericHandler(proxy, port_properties=port_props)
        ref = {"args": numeric.function.arg_structure,
               "res": numeric.function.result_structure}
        assert_reproduction(ref, suffix=f"{port_props}".lower())


def check_same_keys(dic1, dic2):
    """Check whether the two nested dictionaries have the same keys"""
    def is_it(d):
        try:
            d.items()
        except AttributeError:
            return False
        return True

    assert is_it(dic1) == is_it(dic2)
    if not is_it(dic1):
        return
    assert set(dic1.keys()) == set(dic2.keys())
    for key, child in dic1.items():
        check_same_keys(child, dic2[key])
