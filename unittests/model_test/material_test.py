"""Unit tests for material related objects"""

from pytest import raises

from simu.materials import ExampleThermoFactory, ThermoPropertyStore
from simu.model.material import MaterialDefinition
from simu.utilities.errors import DimensionalityError


# class BigNAugmentor(Augmentor):
#     def define(self, material):
#         material["N"] = sum(material["n"])


def test_create_frame():
    factory = ExampleThermoFactory()
    assert "Water-RK-Liquid" in factory.configuration_names
    frame = factory.create_frame("Water-RK-Liquid")
    print(frame.parameter_structure)


def test_get_thermo_properties():
    factory = ExampleThermoFactory()
    frame = factory.create_frame("Water-RK-Liquid")
    store = ThermoPropertyStore()
    symbols = store.get_symbols(frame.parameter_structure)
    assert symbols["CriticalParameters"]["T_c"]["H2O"].units == "kelvin"


def test_get_thermo_properties_twice():
    factory = ExampleThermoFactory()
    frame = factory.create_frame("Water-RK-Liquid")
    store = ThermoPropertyStore()
    symbols1 = store.get_symbols(frame.parameter_structure)
    symbols2 = store.get_symbols(frame.parameter_structure)
    qty1 = symbols1["CriticalParameters"]["T_c"]["H2O"]
    qty2 = symbols2["CriticalParameters"]["T_c"]["H2O"]
    assert qty1 is qty2


def test_get_thermo_properties_twice_wrong():
    store = ThermoPropertyStore()
    struct_1 = {"A": {"B": "K"}}
    struct_2 = {"A": {"B": "m"}}
    store.get_symbols(struct_1)
    with raises(DimensionalityError):
        store.get_symbols(struct_2)


# def test_create_material_definition():
#     factory = ExampleThermoFactory()
#     mat_def = MaterialDefinition()
