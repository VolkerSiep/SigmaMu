from simu.core.solver.thermofit import DataSet

def test_instantiate(example_thermo_fit_configuration):
    data_set_config = example_thermo_fit_configuration["datasets"]["vle"]
    ds = DataSet.model_validate(data_set_config)
    print(ds)