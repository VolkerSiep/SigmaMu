from simu import (InitialState, MaterialDefinition, ThermoParameterStore,
                  SpeciesDefinition, ThermoFrame)
from simu.app import RegThermoFactory, ThermoStructure, predefined_parameters


def _create_frames() -> tuple[ThermoFrame, ThermoFrame]:
    factory = RegThermoFactory()
    species_def = {"H2O": SpeciesDefinition("H2O")}

    struct = ThermoStructure.from_predefined("IAPWS-Liquid")
    cond = factory.create_frame(species_def, struct + "GenericProperties")

    struct = ThermoStructure.from_predefined("IAPWS-Gas")
    steam = factory.create_frame(species_def, struct + "GenericProperties")
    return cond, steam


def _create_material(frame: ThermoFrame, parameters: ThermoParameterStore,
                     t: float, p: float, n: float) -> MaterialDefinition:
    initial_state = InitialState.from_si(t, p, [n])
    return MaterialDefinition(frame, initial_state, parameters)


_f_cond, _f_steam = _create_frames()
_parameters = predefined_parameters()

hp_steam = _create_material(_f_steam, _parameters, 600, 100e5, 100)
hp_condensate = _create_material(_f_cond, _parameters, 600, 100e5, 100)
lp_steam = _create_material(_f_steam, _parameters, 373, 1e5, 100)
lp_condensate = _create_material(_f_cond, _parameters, 373, 1e5, 100)
