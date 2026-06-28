from pathlib import Path
from yaml import safe_load, safe_dump
from simu import (
    SpeciesDB, StringDictThermoSource, ThermoParameterStore,
    InitialState, MaterialDefinition
)
from simu.app import RegThermoFactory
from simu.examples.steam_system.thermo import lp_condensate, lp_steam

CURRENT_DIR = Path(__file__).parent
CONFIG_FILE = CURRENT_DIR / "thermo_config.yml"


def _create_materials(store: ThermoParameterStore):
    # read configuration
    with CONFIG_FILE.open() as file:
        configuration = safe_load(file)

    factory = RegThermoFactory()
    species = SpeciesDB(configuration["species"])

    # define frames and material definitions
    result = {}
    for name, phase in configuration["phases"].items():
        s = species.get_sub_db(phase["species"])
        c = configuration["frames"][phase["frame"]]
        result[name] = MaterialDefinition(
            factory.create_frame(s, c),
            InitialState.from_std(len(s)), store
        )

    for name, data in configuration["parameters"].items():
        source = StringDictThermoSource(data)
        store.add_source(name, source)

    missing = store.get_missing_symbols()
    if missing:
        print("Missing thermodynamic parameters:")
        print(safe_dump(missing))
        exit()

    return result

thermo_store = ThermoParameterStore()
_materials = _create_materials(thermo_store)
condensate_new = _materials["liquid"]  # to be more explicit
steam_new = _materials["gas"]