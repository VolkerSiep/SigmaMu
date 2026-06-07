from collections.abc import Mapping
from pathlib import Path
from yaml import safe_load
from simu import (
    SpeciesDB, StringDictThermoSource, ThermoParameterStore,
    InitialState, MaterialDefinition
)
from simu.app import RegThermoFactory

CURRENT_DIR = Path(__file__).parent
CONFIG_FILE = CURRENT_DIR / "thermo_config.yml"


def _create_materials():
    # read configuration
    with CONFIG_FILE.open() as file:
        configuration = safe_load(file)

    factory = RegThermoFactory()
    species = SpeciesDB(configuration["species"])
    store = ThermoParameterStore()

    # define frames and material definitions
    result = {}
    for name, phase in configuration["phases"].items():
        s = species.get_sub_db(phase["species"])
        c = configuration["frames"][phase["frame"]]
        i = InitialState.from_std(len(s))
        f = factory.create_frame(s, c)
        result[name] = MaterialDefinition(f, i, store)

    for name, data in configuration["parameters"].items():
        source = StringDictThermoSource(data)
        store.add_source(name, source)

    return result


materials: Mapping[str, MaterialDefinition] = _create_materials()
