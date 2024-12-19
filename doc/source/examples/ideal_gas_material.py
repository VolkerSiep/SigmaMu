from pprint import pprint
from yaml import safe_load

from simu import (
    ThermoFactory, InitialState, MaterialDefinition, ThermoParameterStore,
    StringDictThermoSource, SpeciesDB)
from simu.app.thermo import all_contributions, GibbsState

# load species database
with open("species_db.yml") as file:
    species = SpeciesDB(safe_load(file))

# load model structure database
with open("thermo_model_structures.yml") as file:
    model_structures = safe_load(file)

# load thermodynamic parameter database
with open("ideal_gas_param.yml") as file:
    parameter_source = StringDictThermoSource(safe_load(file))

factory = ThermoFactory()
factory.register_state_definition(GibbsState)
factory.register(*all_contributions)

frame = factory.create_frame(species.get_sub_db(["Methane"]),
                             model_structures["simple_ideal_gas"])

initial_state = InitialState.from_si(400, 2e5, [1.0])
store = ThermoParameterStore()
ch4_ideal = MaterialDefinition(frame, initial_state, store)
pprint(store.get_missing_symbols())

store.add_source("my_source", parameter_source)
pprint(store.get_missing_symbols())
