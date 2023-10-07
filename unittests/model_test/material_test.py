"""Unit tests for material related objects"""

from pytest import raises

from simu.thermo import InitialState
from simu.thermo.species import SpeciesDefinition
from simu.thermo.factory import ExampleThermoFactory
from simu.thermo.parameters import ThermoParameterStore
from simu.thermo.material import MaterialDefinition
from simu.utilities import assert_reproduction

# class BigNAugmentor(Augmentor):
#     def define(self, material):
#         material["N"] = sum(material["n"])

RK_LIQ = "Boston-Mathias-Redlich-Kwong-Liquid"


def test_create_material_definition():
    factory = ExampleThermoFactory()
    species = {"H2O": SpeciesDefinition("H2O")}
    frame = factory.create_frame(species, RK_LIQ)
    store = ThermoParameterStore()
    initial_state = InitialState.from_std(1)
    _ = MaterialDefinition(frame, initial_state, store)


def test_create_material_definition_wrong_init():
    factory = ExampleThermoFactory()
    species = {"H2O": SpeciesDefinition("H2O")}
    frame = factory.create_frame(species, RK_LIQ)
    store = ThermoParameterStore()
    initial_state = InitialState.from_std(2)
    with raises(ValueError) as err:
        _ = MaterialDefinition(frame, initial_state, store)
    assert "Incompatible initial state" in str(err.value)


def test_create_material():
    factory = ExampleThermoFactory()
    species = {"H2O": SpeciesDefinition("H2O")}
    frame = factory.create_frame(species, RK_LIQ)
    store = ThermoParameterStore()
    initial_state = InitialState.from_std(1)
    material_def = MaterialDefinition(frame, initial_state, store)
    material = material_def.create_state()
    res = {name: f"{value.units:~}" for name, value in material.items()}
    assert_reproduction(res)
