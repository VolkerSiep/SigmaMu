from pathlib import Path
from yaml import safe_load
from simu import (
    SpeciesDefinition, StringDictThermoSource, ThermoParameterStore,
    InitialState, MaterialDefinition
)
from simu.app import RegThermoFactory

_PARAMETER_FILE = Path(__file__).parent / "wagman_tin.yml"
_SOURCE_ID = "Wagman_1982"
_CONFIG = {
    "state": "GibbsState",
    "contributions": ["H0S0ReferenceState", "LinearHeatCapacity"],
}
_SPECIES = {"alpha-tin": SpeciesDefinition("Sn"),
           "beta-tin": SpeciesDefinition("Sn")}


def _create_material():
    # create thermodynamic model
    factory = RegThermoFactory()
    frame = factory.create_frame(_SPECIES, _CONFIG)
    initial_state = InitialState.from_std(2)  # 25 degC, 1 bar, [1 mol, 1 mol]

    # create thermo store and add parameter source
    store = ThermoParameterStore()
    with open(_PARAMETER_FILE) as file:
        parameters = safe_load(file)
    source = StringDictThermoSource(parameters)
    store.add_source(_SOURCE_ID, source)

    material_definition = MaterialDefinition(frame, initial_state, store)
    return material_definition, source

tin_definition, thermo_source = _create_material()

