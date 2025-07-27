from yaml import safe_load
from simu import ThermoParameterStore, StringDictThermoSource
from simu.app import DATA_DIR


def _populate_store() -> ThermoParameterStore:
    store = ThermoParameterStore()

    for path in (DATA_DIR / "parameters").glob("*.yml"):
        with open(path) as file:
            data = safe_load(file)
            parameter_source = StringDictThermoSource(data["data"])
        store.add_source(data["meta"]["source"], parameter_source)
    return store

predefined_parameters: ThermoParameterStore = _populate_store()
