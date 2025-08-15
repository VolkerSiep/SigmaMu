from yaml import safe_load
from pytest import fixture


from simu import (MaterialDefinition, SpeciesDB, StringDictThermoSource,
                  InitialState, ThermoParameterStore)
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
def liq_and_gas_with_param(thermo_store):
    with open(SPECIES_FILE) as file:
        species = safe_load(file)
    species = SpeciesDB(species).get_sub_db(['C3', 'nC4'])

    factory = RegThermoFactory()
    s_gas = ThermoStructure.from_predefined("Boston-Mathias-Redlich-Kwong-Gas")
    s_liq = ThermoStructure.from_predefined("Boston-Mathias-Redlich-Kwong-Liquid")
    f_gas = factory.create_frame(species, s_gas)
    f_liq = factory.create_frame(species, s_liq)

    initial_state = InitialState.from_cbar(25, 2, [0.5, 0.5])
    m_gas = MaterialDefinition(f_gas, initial_state, thermo_store)
    m_liq = MaterialDefinition(f_liq, initial_state, thermo_store)
    return m_liq, m_gas
