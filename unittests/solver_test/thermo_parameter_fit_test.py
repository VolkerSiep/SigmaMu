from math import isnan

from pytest import raises
from pydantic import ValidationError

from simu import NumericHandler, ThermoFitSolver
from simu.core.solver.thermofit.config import (
    DataSet, ThermoFitContribution, ThermoFitEvaluation, ThermoFitParameter,
    ThermoFitDefinition
)
from simu.core.solver.thermofit.solver import ModelContext
from simu.examples.hello_world import Square
from simu.examples.tin_parameter_fit.thermo import thermo_source
from simu.examples.tin_parameter_fit.simulation import TinTransition
from simu.examples.tin_parameter_fit.parameter_fit import load_definition


def test_instantiate_dataset(thermo_fit_configuration):
    data_set_config = thermo_fit_configuration["datasets"]["vle_1bar"]
    ds = DataSet.model_validate(data_set_config)
    assert ds.data[1][3] == 50.0


def test_dataset_invalid_uom(thermo_fit_configuration):
    data_set_config = thermo_fit_configuration["datasets"]["vle_1bar"]
    data_set_config["uom"][1] = "hansi"
    with raises(ValidationError) as e:
        DataSet.model_validate(data_set_config)
    assert "hansi" in str(e)


def test_dataset_invalid_float(thermo_fit_configuration):
    data_set_config = thermo_fit_configuration["datasets"]["vle_1bar"]
    data_set_config["data"][1][1] = "hansi"
    with raises(ValidationError) as e:
        DataSet.model_validate(data_set_config)
    assert "hansi" in str(e)


def test_dataset_wong_num_uom(thermo_fit_configuration):
    data_set_config = thermo_fit_configuration["datasets"]["vle_1bar"]
    del data_set_config["uom"][-1]
    with raises(ValidationError) as e:
        DataSet.model_validate(data_set_config)
    assert "units" in str(e)


def test_dataset_wong_num_data(thermo_fit_configuration):
    data_set_config = thermo_fit_configuration["datasets"]["vle_1bar"]
    del data_set_config["data"][3][-1]
    with raises(ValidationError) as e:
        DataSet.model_validate(data_set_config)
    assert "values in row 3" in str(e)


def test_instantiate_thermo_fit_contribution(
        thermo_fit_configuration, contribution_context_stub):
    config = thermo_fit_configuration["contributions"]["vle"]
    contribution = ThermoFitContribution.model_validate(
        config, context=contribution_context_stub
    )
    assert contribution.model_id == "vle_fit"


def test_thermo_fit_contribution_wrong_model(
        thermo_fit_configuration, contribution_context_stub):
    config = thermo_fit_configuration["contributions"]["vle"]
    config["model_id"] = "hansi"
    with raises(ValidationError) as e:
        ThermoFitContribution.model_validate(
            config, context=contribution_context_stub
        )
    assert "hansi" in str(e)


def test_thermo_fit_contribution_wrong_model_parameter(
        thermo_fit_configuration, contribution_context_stub):
    config = thermo_fit_configuration["contributions"]["vle"]
    config["data_to_model"]["x"] = "hansi"
    with raises(ValidationError) as e:
        ThermoFitContribution.model_validate(
            config, context=contribution_context_stub
        )
    assert "hansi" in str(e)


def test_thermo_fit_contribution_wrong_model_property(
        thermo_fit_configuration, contribution_context_stub):
    config = thermo_fit_configuration["contributions"]["vle"]
    config["penalties"][1] = "hansi"
    with raises(ValidationError) as e:
        ThermoFitContribution.model_validate(
            config, context=contribution_context_stub
        )
    assert "hansi" in str(e)

def test_thermo_fit_contribution_uom_penalty(
        thermo_fit_configuration, contribution_context_stub):
    config = thermo_fit_configuration["contributions"]["vle"]
    config["penalties"].append("process.p")
    with raises(ValidationError) as e:
        ThermoFitContribution.model_validate(
            config, context=contribution_context_stub
        )
    assert "process.p" in str(e)

def test_thermo_fit_evaluation(
        thermo_fit_configuration, contribution_context_stub):
    config = thermo_fit_configuration["evaluations"]["vle_p"]
    evaluation = ThermoFitEvaluation.model_validate(
        config, context=contribution_context_stub
    )
    assert "p_hat" in evaluation.properties

def test_thermo_fit_evaluation_wrong_parameter(
        thermo_fit_configuration, contribution_context_stub):
    config = thermo_fit_configuration["evaluations"]["vle_p"]
    config["properties"]["p_hat"]["name"] = "hansi"
    with raises(ValidationError) as e:
        ThermoFitEvaluation.model_validate(
            config, context=contribution_context_stub
        )
    assert "hansi" in str(e)

def test_thermo_fit_evaluation_wrong_uom(
        thermo_fit_configuration, contribution_context_stub):
    config = thermo_fit_configuration["evaluations"]["vle_p"]
    config["properties"]["p_hat"]["uom"] = "ft^2"
    with raises(ValidationError) as e:
        ThermoFitEvaluation.model_validate(
            config, context=contribution_context_stub
        )
    assert "ft^2" in str(e)

def test_thermo_fit_parameter(
        thermo_fit_configuration, contribution_context_stub):
    config = thermo_fit_configuration["parameters"]["a_b"]
    parameter = ThermoFitParameter.model_validate(
        config, context=contribution_context_stub
    )
    assert parameter.lower.magnitude == 200.0

def test_thermo_fit_parameter_wrong_sequence(
        thermo_fit_configuration, contribution_context_stub):
    config = thermo_fit_configuration["parameters"]["a_b"]
    config["lower"] = "400 K"
    with raises(ValidationError) as e:
        ThermoFitParameter.model_validate(
            config, context=contribution_context_stub
        )
    assert "400 K" in str(e)

def test_thermo_fit_parameter_different_units(
        thermo_fit_configuration, contribution_context_stub):
    config = thermo_fit_configuration["parameters"]["a_b"]
    config["lower"] = "400 m"
    with raises(ValidationError) as e:
        ThermoFitParameter.model_validate(
            config, context=contribution_context_stub
        )
    assert "400 m" in str(e)

def test_thermo_fit_parameter_wrong_name(
        thermo_fit_configuration, contribution_context_stub):
    config = thermo_fit_configuration["parameters"]["a_b"]
    config["path"] = ["hansi"]
    with raises(ValidationError) as e:
        ThermoFitParameter.model_validate(
            config, context=contribution_context_stub
        )
    assert "hansi" in str(e)

def test_thermo_fit_parameter_wrong_unit(
        thermo_fit_configuration, contribution_context_stub):
    config = thermo_fit_configuration["parameters"]["a_b"]
    config["default"] = "30 m"
    with raises(ValidationError) as e:
        ThermoFitParameter.model_validate(
            config, context=contribution_context_stub
        )
    assert "30 m" in str(e)

def test_thermo_fit_config(thermo_fit_configuration, contribution_context_stub):
    config = ThermoFitDefinition.model_validate(
        thermo_fit_configuration, context=contribution_context_stub
    )
    assert config.evaluations["vle_p"].data_to_model["x"].path == ["x"]

def test_model_contest():
    model = NumericHandler(Square.top())
    context = ModelContext(model)
    assert context.parameter_unit(["length"]) == "m"
    assert context.property_unit(["area"]) == "m ** 2"
    with raises(KeyError) as err:
        context.property_unit(["area", "Antarctica"])
    assert "Antarctica" in str(err)

def test_prepare_function_r(tin_functions):
    x = [273.15 + 8.83, 1e5, 1, 1]
    r, jac = tin_functions.f_r(x, [12.3], [44.14])
    assert r.shape == (4, 1)
    assert jac.shape == (4, 4)
    r = list(r.nonzeros())
    for i in range(1, 4):
        assert r[i] == 0
    assert 160 < r[0] < 170

def test_prepare_function_q(tin_functions):
    x = [273.15 + 8.83, 1e5, 1, 1]
    q, q_x, q_t, r_x, r_t = tin_functions.f_q(x, [12.3], [44.14])
    assert abs(q) < 0.02  # close to solution with 8.83 degC
    assert q_x.shape == (1, 4)
    assert r_t.shape == (4, 1)
    for i in range(1, 4):
        assert q_x[i] == 0.0  # q only depends on T, not p or n_i
    assert q_t == 0.0  # no direct dependency
    for i in range(1, 4):
        assert r_t[i] == 0.0  # only r[0] depends on thermo-parameter

def test_prepare_function_bx(tin_functions):
    temp = 273.15 + 8.83
    x = [temp, 1e5, 1, 1]
    dx = [-2 * temp, 0, 0, 0]
    b, a = tin_functions.f_bx(x, [12.3], [44.14], dx)
    assert b == temp
    assert a == 0.5

def test_prepare_function_bt(tin_functions):
    temp = 273.15 + 8.83
    x = [temp, 1e5, 1, 1]
    b, a = tin_functions.f_bt(x, [12.3], [44.14], [1])
    assert b == temp
    assert isnan(a)

def test_tin_parameter_fit():
    models = {"transition_model": NumericHandler(TinTransition.top())}
    solver = ThermoFitSolver(models, thermo_source)
    report = solver.solve(load_definition())
    new_param = report.thermo_source
