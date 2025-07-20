from simu.app import RegThermoFactory

factory = RegThermoFactory()

from simu import SpeciesDefinition

config = {
    "state": "GibbsState",
    "contributions": [
        "H0S0ReferenceState",
        "LinearHeatCapacity",
        "IdealMix",
        "GibbsIdealGas"
    ],
}
species = {"Methane": SpeciesDefinition("CH4")}
frame = factory.create_frame(species, config)

from simu import parse_quantities_in_struct

parameters = parse_quantities_in_struct({
    'H0S0ReferenceState': {
        'T_ref': '25 degC',
        'dh_form': {'Methane': '-74.87 kJ/mol'},
        'p_ref': '1 bar',
        's_0': {'Methane': '188.66 J/K/mol'}},
    'LinearHeatCapacity': {
        'cp_a': {'Methane': '35.69 J/K/mol'},
        'cp_b': {'Methane': '50 mJ/K**2/mol'}}
})

result = frame([400, 1e5, 1.0], parameters)

from simu import InitialState

tpn = InitialState.from_si(400, 2e5, [1.0])
state = frame.initial_state(tpn, parameters)
result2 = frame(state, parameters)
