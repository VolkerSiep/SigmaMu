from pytest import raises
from pydantic import ValidationError
from simu.core.solver.thermofit import DataSet, ThermoFitContribution


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


def test_instantiate_thermo_fit_contribution_wrong_model(
        thermo_fit_configuration, contribution_context_stub):
    config = thermo_fit_configuration["contributions"]["vle"]
    config["model_id"] = "hansi"
    with raises(ValidationError) as e:
        ThermoFitContribution.model_validate(
            config, context=contribution_context_stub
        )
    assert "hansi" in str(e)


def test_instantiate_thermo_fit_contribution_wrong_model_parameter(
        thermo_fit_configuration, contribution_context_stub):
    config = thermo_fit_configuration["contributions"]["vle"]
    config["data_to_model"]["x"] = "hansi"
    with raises(ValidationError) as e:
        ThermoFitContribution.model_validate(
            config, context=contribution_context_stub
        )
    assert "hansi" in str(e)


def test_instantiate_thermo_fit_contribution_wrong_model_property(
        thermo_fit_configuration, contribution_context_stub):
    config = thermo_fit_configuration["contributions"]["vle"]
    config["penalties"][1] = "hansi"
    with raises(ValidationError) as e:
        ThermoFitContribution.model_validate(
            config, context=contribution_context_stub
        )
    assert "hansi" in str(e)