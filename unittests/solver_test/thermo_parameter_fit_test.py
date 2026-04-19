from pytest import raises
from pydantic import ValidationError
from simu.core.solver.thermofit import (
    DataSet, ThermoFitContribution, ThermoFitEvaluation, ThermoFitParameter
)


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
