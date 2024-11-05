from simu import ThermoFactory, SpeciesDefinition
from simu.app.thermo import all_contributions, GibbsState

factory = ThermoFactory()
factory.register_state_definition(GibbsState)
factory.register(*all_contributions)

config = {
    "state": "GibbsState",
    "contributions": [
        "H0S0ReferenceState",
        "LinearHeatCapacity",
        "StandardState",
        "IdealMix",
        "GibbsIdealGas"
    ],
}
species = {"Methane": SpeciesDefinition("CH4")}
frame = factory.create_frame(species, config)

from pprint import pprint
pprint(frame.parameter_structure, width=90)
print()
pprint(frame.property_structure)
