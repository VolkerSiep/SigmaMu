from numpy import squeeze
from numpy.testing import assert_allclose

from simu import NumericHandler, flatten_dictionary
from simu.core.utilities import assert_reproduction

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
    assert args[NumericHandler.MODEL_PARAMS]['length'] == 'm'


def test_properties():
    proxy = PropertyTestModel.top()
    numeric = NumericHandler(proxy)
    results = numeric.function.result_structure
    assert results[NumericHandler.MODEL_PROPS]['area'] == 'm ** 2'


def test_residuals():
    proxy = ResidualTestModel.top()
    numeric = NumericHandler(proxy)
    results = numeric.function.result_structure
    assert results[NumericHandler.RESIDUALS]['area'] == "m ** 2"


def test_material_collect_states():
    args, _ = create_material_functions()
    assert args["vectors"][NumericHandler.STATE_VEC] == ""


def test_material_collect_multiple_states():
    proxy = MaterialTestModel4.top()
    numeric = NumericHandler(proxy)
    state = numeric.arguments["vectors"][NumericHandler.STATE_VEC]
    assert len(state.magnitude.nz) == 6


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


def test_jacobian():
    model = SquareTestModel()
    material = model.no2sol
    numeric = NumericHandler(model.create_proxy().finalise())
    material.store.add_source("default", StringDictThermoSource(DATA))
    jac_id = numeric.register_jacobian(NumericHandler.RES_VEC,
                                       NumericHandler.STATE_VEC)
    args = numeric.arguments
    result = numeric.function(args)
    dr_dx = result[NumericHandler.JACOBIANS][jac_id].magnitude
    assert_reproduction(dr_dx.tolist())


def test_collect_hierarchy_material():
    proxy = MaterialParentTestModel.top()
    for port_props in (True, False):
        numeric = NumericHandler(proxy, port_properties=port_props)
        ref = {"args": numeric.function.arg_structure,
               "res": numeric.function.result_structure}
        assert_reproduction(ref, suffix=f"{port_props}".lower())


def test_extract_parameters():
    model = SquareTestModel()
    store = model.no2sol.store
    numeric = NumericHandler(model.create_proxy().finalise())
    store.add_source("default", StringDictThermoSource(DATA))
    params = {'model_params': {'N': 'mol / s', 'T': '°C',
                               'p': 'bar', 'x_c3': '%'}}
    numeric.extract_parameters("param", params)
    names = numeric.vector_arg_names("param")
    jac_id = numeric.register_jacobian(NumericHandler.RES_VEC, "param")
    args = numeric.arguments
    result = numeric.function(args)
    dr_dp = result[NumericHandler.JACOBIANS][jac_id].magnitude
    ref = {"names": names, "J": dr_dp.tolist()}
    assert_reproduction(ref)


def test_collect_properties():
    model = SquareTestModel()
    store = model.no2sol.store
    numeric = NumericHandler(model.create_proxy().finalise())
    store.add_source("default", StringDictThermoSource(DATA))
    props = {'thermo_props': {'local': {'mu': {
        'CH3-(CH2)2-CH3': 'kJ/mol', 'CH3-CH2-CH3': 'kJ/mol'}}}}
    numeric.collect_properties("mu", props)
    names = numeric.vector_res_names("mu")
    jac_id = numeric.register_jacobian("mu", NumericHandler.STATE_VEC)
    args = numeric.arguments
    result = numeric.function(args)
    dmu_dx = result[NumericHandler.JACOBIANS][jac_id].magnitude
    ref = {"names": names, "J": dmu_dx.tolist()}
    assert_reproduction(ref)


def test_export_state():
    numeric = NumericHandler(SquareTestModel.top())
    state = numeric.export_state()
    assert_reproduction(state)


def test_import_state():
    model = SquareTestModel()
    numeric = NumericHandler(model.create_proxy().finalise())
    state = {
        'thermo': {'local': {
            'T': '100 °C', 'p': '5 bar',
            'n': {'CH3-CH2-CH3': '2 mol', 'CH3-(CH2)2-CH3': '1 mol'}}},
        'non-canonical': {}}
    species = list(state["thermo"]["local"]["n"].keys())
    numeric.import_state(state)
    state = model.materials["local"].initial_state.to_dict(species)
    assert_reproduction(state)


def test_retain_initial_values():
    model = SquareTestModel()
    numeric = NumericHandler(model.create_proxy().finalise())
    material = model.materials["local"]
    material.definition.store.add_source("default", StringDictThermoSource(DATA))
    params = numeric.arguments["thermo_params"]["default"]
    state = [283.15, 2 * 0.000196732, 2, 2]
    numeric.retain_initial_values(state, params)
    pressure = material.initial_state.pressure
    assert Quantity(0.999, "MPa") < pressure < Quantity(1.001, "MPa")


def test_retain_and_args():
    model = SquareTestModel()
    numeric = NumericHandler(model.create_proxy().finalise())
    material = model.materials["local"]
    material.definition.store.add_source("default", StringDictThermoSource(DATA))
    params = numeric.arguments["thermo_params"]["default"]
    state = [283.15, 2 * 0.000196732, 2, 2]
    numeric.retain_initial_values(state, params)
    new_state  = squeeze(numeric.arguments["vectors"]["states"].magnitude)
    assert_allclose(new_state, state)


def create_material_functions():
    """Make a function out of a model defining materials"""
    proxy = MaterialTestModel3.top()
    numeric = NumericHandler(proxy)
    args = numeric.function.arg_structure
    results = numeric.function.result_structure
    return args, results


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
