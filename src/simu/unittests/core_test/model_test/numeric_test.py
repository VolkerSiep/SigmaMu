from numpy import squeeze
from numpy.testing import assert_allclose

from simu import NumericHandler, flatten_dictionary, Quantity, jacobian
from simu.examples.material_model import Source
from simu.core.utilities import assert_reproduction

from .models import *

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


def test_material_collect_states(material_model_function):
    args = material_model_function[0]
    assert args["vectors"][NumericHandler.STATE_VEC] == ""


def test_material_collect_multiple_states(material_test_model_4):
    proxy = material_test_model_4.top()
    numeric = NumericHandler(proxy)
    state = numeric.arguments["vectors"][NumericHandler.STATE_VEC]
    assert len(state.magnitude.nz) == 6


def test_material_collect_props(material_model_function):
    results = material_model_function[1]
    assert_reproduction(results["thermo_props"]["local"])


def test_material_collect_thermo_param(material_model_function):
    args = material_model_function[0]
    assert_reproduction(args["thermo_params"]["default"])


def test_hierarchy_collect_numerics():
    numeric = NumericHandler(HierarchyTestModel2.top())
    results = numeric.function.result_structure
    assert "area" in results["model_props"]["square"]


def test_square_model(square_test_model):
    numeric = NumericHandler(square_test_model.top())
    ref = {"args": numeric.function.arg_structure,
           "res": numeric.function.result_structure}
    assert_reproduction(ref)


def test_square_model_args(thermo_param, square_test_model):
    model = square_test_model()
    material = model.no2sol
    numeric = NumericHandler(model.create_proxy().finalise())
    material.store.add_source("default", thermo_param)
    struct = numeric.function.arg_structure
    args = numeric.arguments
    check_same_keys(struct, args)


def test_square_model_call(thermo_param, square_test_model):
    model = square_test_model()
    material = model.no2sol
    numeric = NumericHandler(model.create_proxy().finalise())
    material.store.add_source("default", thermo_param)
    args = numeric.arguments
    res = flatten_dictionary(numeric.function(args))
    res = {k: f"{v:.6f~}" for k, v in res.items()}
    assert_reproduction(res)


def test_jacobian(thermo_param, square_test_model):
    model = square_test_model()
    material = model.no2sol
    numeric = NumericHandler(model.create_proxy().finalise())
    material.store.add_source("default", thermo_param)
    jac_id = numeric.register_jacobian(NumericHandler.RES_VEC,
                                       NumericHandler.STATE_VEC)
    args = numeric.arguments
    result = numeric.function(args)
    dr_dx = result[NumericHandler.JACOBIANS][jac_id].magnitude
    assert_reproduction(dr_dx.tolist())


def test_collect_hierarchy_material(material_parent_test_model):
    proxy = material_parent_test_model.top()
    for port_props in (True, False):
        numeric = NumericHandler(proxy, port_properties=port_props)
        ref = {"args": numeric.function.arg_structure,
               "res": numeric.function.result_structure}
        assert_reproduction(ref, suffix=f"{port_props}".lower())


def test_extract_parameters(thermo_param, square_test_model):
    model = square_test_model()
    store = model.no2sol.store
    numeric = NumericHandler(model.create_proxy().finalise())
    store.add_source("default", thermo_param)
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


def test_collect_properties(thermo_param, square_test_model):
    model = square_test_model()
    store = model.no2sol.store
    numeric = NumericHandler(model.create_proxy().finalise())
    store.add_source("default", thermo_param)
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


def test_export_state(square_test_model):
    numeric = NumericHandler(square_test_model.top())
    state = numeric.export_state()
    assert_reproduction(state)


def test_import_state(square_test_model):
    model = square_test_model()
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


def test_retain_initial_values(thermo_param, square_test_model):
    model = square_test_model()
    numeric = NumericHandler(model.create_proxy().finalise())
    material = model.materials["local"]
    material.definition.store.add_source("default", thermo_param)
    params = numeric.arguments["thermo_params"]
    state = [283.15, 2 * 0.000196732, 2, 2]
    numeric.retain_state(state, params)
    pressure = material.initial_state.pressure
    assert Quantity(0.999, "MPa") < pressure < Quantity(1.001, "MPa")


def test_retain_and_args(thermo_param, square_test_model):
    model = square_test_model()
    numeric = NumericHandler(model.create_proxy().finalise())
    material = model.materials["local"]
    material.definition.store.add_source("default", thermo_param)
    params = numeric.arguments["thermo_params"]
    state = [283.15, 2 * 0.000196732, 2, 2]
    numeric.retain_state(state, params)
    new_state  = squeeze(numeric.arguments["vectors"]["states"].magnitude)
    assert_allclose(new_state, state)

def test_thermo_residual(model_with_residual):
    numeric = NumericHandler(model_with_residual.top())
    rs = numeric.function.result_structure
    assert rs[numeric.RESIDUALS]["liq"]["ChargeBalance"]["balance"] == "A"


def test_query_bounds():
    numeric = NumericHandler(Source.top())
    res = numeric.vector_res_names(numeric.BOUND_VEC)
    assert_reproduction(res)

def test_model_bounds():
    numeric = NumericHandler(BoundTestModel.top())
    res = numeric.vector_res_names(numeric.BOUND_VEC)
    assert_reproduction(res)


def test_bound_sensitivity():
    numeric = NumericHandler(Source.top())
    args = numeric.arguments
    names = numeric.vector_arg_names(numeric.STATE_VEC)
    state = SymbolQuantity("x", "", names)
    args[numeric.VECTORS][numeric.STATE_VEC] = state
    res = numeric.function(args, squeeze_results=False)
    res = res[numeric.VECTORS][numeric.BOUND_VEC]
    jac = jacobian(res, state).magnitude
    assert_reproduction(str(jac))


def test_vector_bound(square_test_model):
    numeric = NumericHandler(square_test_model.top())
    res = numeric.vector_res_names(numeric.BOUND_VEC)
    res = [r for r in res if r.startswith("local/IdealMix/")]
    ref = ["local/IdealMix/n/CH3-(CH2)2-CH3", "local/IdealMix/n/CH3-CH2-CH3"]
    assert res == ref


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

