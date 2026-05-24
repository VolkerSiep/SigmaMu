from numpy import squeeze
from numpy.testing import assert_allclose

from simu import (
    NumericHandler, NHKeys, flatten_dictionary, Quantity, jacobian,
    PropertyFilter)
from simu.app.numeric import ExclusionFilter
from simu.examples.material_model import Source
from simu.core.utilities.testing import assert_reproduction

from .models import *

def test_parameters():
    proxy = SimpleParameterTestModel.top()
    numeric = NumericHandler(proxy)
    args = numeric.function.arg_structure
    assert args[NHKeys.MODEL_PARAMS]['length'] == 'm'


def test_properties():
    proxy = PropertyTestModel.top()
    numeric = NumericHandler(proxy)
    results = numeric.function.result_structure
    assert results[NHKeys.MODEL_PROPS]['area'] == 'm ** 2'


def test_residuals():
    proxy = ResidualTestModel.top()
    numeric = NumericHandler(proxy)
    results = numeric.function.result_structure
    assert results[NHKeys.RESIDUALS]['area'] == "m ** 2"


def test_material_collect_states(material_model_function):
    args = material_model_function[0]
    assert args[NHKeys.VECTORS][NHKeys.STATES] == ""


def test_material_collect_multiple_states(material_test_model_4):
    proxy = material_test_model_4.top()
    numeric = NumericHandler(proxy)
    state = numeric.arguments[NHKeys.VECTORS][NHKeys.STATES]
    assert len(state.magnitude.nz) == 6


def test_material_collect_props(material_model_function):
    results = material_model_function[1]
    assert_reproduction(results[NHKeys.THERMO_PROPS]["local"])


def test_material_collect_thermo_param(material_model_function):
    args = material_model_function[0]
    assert_reproduction(args[NHKeys.THERMO_PARAMS]["default"])


def test_hierarchy_collect_numerics():
    numeric = NumericHandler(HierarchyTestModel2.top())
    results = numeric.function.result_structure
    assert "area" in results[NHKeys.MODEL_PROPS]["square"]


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


def test_collect_hierarchy_material(material_parent_test_model):
    proxy = material_parent_test_model.top()
    for port_props in (True, False):
        numeric = NumericHandler(proxy, port_properties=port_props)
        ref = {"args": numeric.function.arg_structure,
               "res": numeric.function.result_structure}
        assert_reproduction(ref, suffix=f"{port_props}".lower())


def test_filter_properties(square_test_model):
    class Filter(PropertyFilter):
        def keep_property(self, name: str, sub_key: str = None) -> bool:
            return not (name.endswith("_std") or name.endswith("_ref"))

    proxy = square_test_model.top()
    numeric = NumericHandler(proxy, property_filter=Filter())
    props = numeric.function.result_structure[NHKeys.THERMO_PROPS]["local"]
    assert_reproduction(props)


def test_filter_properties_exclusion(square_test_model):
    filter_ = ExclusionFilter({"mu_std", "S_std", "p_std", "T_ref", "p_ref"})
    proxy = square_test_model.top()
    numeric = NumericHandler(proxy, property_filter=filter_)
    props = numeric.function.result_structure[NHKeys.THERMO_PROPS]["local"]
    assert_reproduction(props)


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
    params = numeric.arguments[NHKeys.THERMO_PARAMS]
    state = [283.15, 2 * 0.000196732, 2, 2]
    numeric.retain_state(state, params)
    pressure = material.initial_state.pressure
    assert Quantity(0.999, "MPa") < pressure < Quantity(1.001, "MPa")


def test_retain_and_args(thermo_param, square_test_model):
    model = square_test_model()
    numeric = NumericHandler(model.create_proxy().finalise())
    material = model.materials["local"]
    material.definition.store.add_source("default", thermo_param)
    params = numeric.arguments[NHKeys.THERMO_PARAMS]
    state = [283.15, 2 * 0.000196732, 2, 2]
    numeric.retain_state(state, params)
    new_state  = squeeze(numeric.arguments[NHKeys.VECTORS][NHKeys.STATES].m)
    assert_allclose(new_state, state)

def test_thermo_residual(model_with_residual):
    numeric = NumericHandler(model_with_residual.top())
    rs = numeric.function.result_structure
    assert rs[NHKeys.RESIDUALS]["liq"]["ChargeBalance"]["balance"] == "A"


def test_query_bounds():
    numeric = NumericHandler(Source.top())
    res = numeric.vector_res_names(NHKeys.BOUNDS)
    assert_reproduction(res)

def test_model_bounds():
    numeric = NumericHandler(BoundTestModel.top())
    res = numeric.vector_res_names(NHKeys.BOUNDS)
    assert_reproduction(res)


def test_bound_sensitivity():
    numeric = NumericHandler(Source.top())
    args = numeric.arguments
    names = numeric.vector_arg_names(NHKeys.STATES)
    state = SymbolQuantity("x", "", names)
    args[NHKeys.VECTORS][NHKeys.STATES] = state
    res = numeric.function(args, squeeze_results=False)
    res = res[NHKeys.VECTORS][NHKeys.BOUNDS]
    jac = jacobian(res, state).magnitude
    assert_reproduction(str(jac))


def test_vector_bound(square_test_model):
    numeric = NumericHandler(square_test_model.top())
    res = numeric.vector_res_names(NHKeys.BOUNDS)
    res = [r for r in res if r.startswith("local/IdealMix/")]
    ref = ["local/IdealMix/n/CH3-(CH2)2-CH3", "local/IdealMix/n/CH3-CH2-CH3"]
    assert res == ref


def test_hierarchy_port(model_with_material_hierarchy):
    numeric = NumericHandler(model_with_material_hierarchy.top())
    x = numeric.arguments["vectors"]["states"].magnitude.nonzeros()
    assert_allclose(x, [298.15, 101325, 1])


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

