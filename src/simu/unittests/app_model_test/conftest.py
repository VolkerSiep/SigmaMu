from yaml import safe_load
from pytest import fixture


from simu import (MaterialDefinition, SpeciesDB, StringDictThermoSource,
                  InitialState, ThermoParameterStore, ThermoFrame)
from simu.app import DATA_DIR, RegThermoFactory, ThermoStructure

SPECIES_FILE = DATA_DIR / "species.yml"
PARAM_FILE = DATA_DIR/ "parameters" / "general_example_parameters.yml"


@fixture(scope="session")
def thermo_store():
    with open(PARAM_FILE) as file:
        params = safe_load(file)
    source = params["meta"]["source"]
    params = StringDictThermoSource(params["data"])
    store = ThermoParameterStore()
    store.add_source(source, params)
    return store

@fixture(scope="session")
def species_db() -> SpeciesDB:
    with open(SPECIES_FILE) as file:
        species = safe_load(file)
    return SpeciesDB(species)

@fixture(scope="session")
def liq_and_gas_with_param(thermo_store, species_db):
    with open(SPECIES_FILE) as file:
        species = safe_load(file)
    species = species_db.get_sub_db(['C3', 'nC4'])

    factory = RegThermoFactory()
    s_gas = ThermoStructure.from_predefined("Boston-Mathias-Redlich-Kwong-Gas")
    s_liq = ThermoStructure.from_predefined("Boston-Mathias-Redlich-Kwong-Liquid")
    f_gas = factory.create_frame(species, s_gas)
    f_liq = factory.create_frame(species, s_liq)

    initial_state = InitialState.from_cbar(25, 2, [0.5, 0.5])
    m_gas = MaterialDefinition(f_gas, initial_state, thermo_store)
    m_liq = MaterialDefinition(f_liq, initial_state, thermo_store)
    return m_liq, m_gas


@fixture(scope="session")
def frame_h2o_c3_gas(species_db) -> ThermoFrame:
    factory = RegThermoFactory()
    struc = ThermoStructure.from_predefined("Boston-Mathias-Redlich-Kwong-Gas")
    return factory.create_frame(species_db.get_sub_db(["H2O", "C3"]), struc)


@fixture(scope="session")
def frame_c3_c4_gas(species_db) -> ThermoFrame:
    factory = RegThermoFactory()
    struc = ThermoStructure.from_predefined("Boston-Mathias-Redlich-Kwong-Gas")
    return factory.create_frame(species_db.get_sub_db(["C3", "nC4"]), struc)


@fixture(scope="session")
def material_h2o_c3_gas(frame_h2o_c3_gas, thermo_store) -> MaterialDefinition:
    initial_state = InitialState.from_std(2)
    return MaterialDefinition(frame_h2o_c3_gas, initial_state, thermo_store)


@fixture(scope="session")
def material_h2o_c4_gas(frame_c3_c4_gas, thermo_store) -> MaterialDefinition:
    initial_state = InitialState.from_std(2)
    return MaterialDefinition(frame_c3_c4_gas, initial_state, thermo_store)

