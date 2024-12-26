from yaml import safe_load
from pathlib import Path

from simu import (
    ThermoFactory, InitialState, MaterialDefinition, ThermoParameterStore,
    StringDictThermoSource, SpeciesDB)
from simu.app.thermo import all_contributions, all_states

CURRENT_DIR = Path(__file__).parent

class MaterialFactory:
    def __init__(self):
        # load species database
        with open(CURRENT_DIR / "species_db.yml") as file:
            self.species = SpeciesDB(safe_load(file))

        # load model structure database
        with open(CURRENT_DIR / "thermo_model_structures.yml") as file:
            self.model_structures = safe_load(file)

        # load thermodynamic parameter database
        with open(CURRENT_DIR / "ideal_gas_param.yml") as file:
            parameter_source = StringDictThermoSource(safe_load(file))

        self.store = ThermoParameterStore()
        self.store.add_source("my_source", parameter_source)

        factory = ThermoFactory()
        for state in all_states:
            factory.register_state_definition(state)
        factory.register(*all_contributions)
        self.factory = factory

    def create(self, model_name, species):
        frame = self.factory.create_frame(
            self.species.get_sub_db(species),
            self.model_structures[model_name])

        initial_state = InitialState.from_si(400, 2e5, [1.0] * len(species))
        return MaterialDefinition(frame, initial_state, self.store)


material_factory = MaterialFactory()

ch4_ideal_gas = material_factory.create("simple_ideal_gas", ["Methane"])