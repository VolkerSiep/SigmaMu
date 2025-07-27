from simu import (
    InitialState, MaterialDefinition, SpeciesDefinition,
    MaterialSpec, ThermoFrame)
from simu.app import RegThermoFactory, ThermoStructure, predefined_parameters


def _create_frames() -> (ThermoFrame, ThermoFrame):
    factory = RegThermoFactory()
    species_def = {"H2O": SpeciesDefinition("H2O")}

    struct = ThermoStructure.from_predefined("IAPWS-Liquid")
    cond = factory.create_frame(species_def, struct + "GenericProperties")

    struct = ThermoStructure.from_predefined("IAPWS-Gas")
    steam = factory.create_frame(species_def, struct + "GenericProperties")
    return cond, steam


def _create_material(frame: ThermoFrame,
                     t: float, p: float, n: float) -> MaterialDefinition:
    initial_state = InitialState.from_si(t, p, [n])
    return MaterialDefinition(frame, initial_state, predefined_parameters)


_f_cond, _f_steam = _create_frames()

hp_steam = _create_material(_f_steam, 600, 100e5, 100)
hp_condensate = _create_material(_f_cond, 600, 100e5, 100)
lp_steam = _create_material(_f_steam, 373, 1e5, 100)
lp_condensate = _create_material(_f_cond, 373, 1e5, 100)

h2o_spec = MaterialSpec(["H2O"])

