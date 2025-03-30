"""Unit tests for material related objects"""

from pytest import raises

from simu import (
    InitialState, ThermoFactory, ThermoParameterStore, MaterialDefinition,
    MaterialSpec, SpeciesDefinition, SpeciesDB, flatten_dictionary, Quantity)
from simu.core.model.material import MaterialHandler
from simu.core.thermo.material import MaterialLab
from simu.core.utilities import assert_reproduction
from simu.app import all_contributions, GibbsState
from simu.app.thermo.factories import ExampleThermoFactory


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


def test_create_material(material_h2o_rk_liq):
    flat = flatten_dictionary(material_h2o_rk_liq)
    res = {name: f"{value.units:~}" for name, value in flat.items()}
    assert_reproduction(res)


def test_retain_initial_state_material(material_def_h2o_with_param):
    material = material_def_h2o_with_param.create_flow()
    store = material_def_h2o_with_param.store
    param = {store.name: store.get_all_values()}
    x = [300, 2e5, 10]
    material.retain_initial_state(x, param)
    ini = material.initial_state
    assert ini.temperature == Quantity(300.0, "K")
    assert ini.pressure == Quantity(2, "bar")
    assert ini.mol_vector == Quantity([10], "mol")


def test_species_db():
    _ = SpeciesDB({"Water": "H2O", "Ethanol": "C2H5OH"})


def test_material_lab():
    species= SpeciesDB({"Water": "H2O", "Ethanol": "C2H5OH"})
    store = ThermoParameterStore()
    factory = ThermoFactory()
    factory.register_state_definition(GibbsState)
    factory.register(*all_contributions)
    lab = MaterialLab(factory, species, store)

    structure = {
        "state": "GibbsState",
        "contributions": [
            "H0S0ReferenceState",
            "LinearHeatCapacity",
            "StandardState",
            "IdealMix",
            "GibbsIdealGas"]
    }
    initial_state = InitialState.from_std(2)
    definition = lab.define_material(species, initial_state, structure)
    material = definition.create_flow()
    assert_reproduction(list(material.keys()))


def test_handler_def_port():
    handler = MaterialHandler()
    handler.define_port("inlet")


def test_handler_create_proxy():
    handler = MaterialHandler()
    handler.define_port("inlet")
    proxy = handler.create_proxy()
    assert "inlet" in proxy.free_ports()


def test_handler_connect_port(material_h2o_rk_liq):
    handler = MaterialHandler()
    handler.define_port("inlet")
    proxy = handler.create_proxy()
    proxy.connect("inlet", material_h2o_rk_liq)


def test_handler_connect_no_port(material_h2o_rk_liq):
    handler = MaterialHandler()
    proxy = handler.create_proxy()
    with raises(KeyError):
        proxy.connect("inlet", material_h2o_rk_liq)


def test_handler_connect_incompatible_port(material_h2o_rk_liq):
    handler = MaterialHandler()
    handler.define_port("inlet", MaterialSpec(flow=False))
    proxy = handler.create_proxy()
    with raises(ValueError):
        proxy.connect("inlet", material_h2o_rk_liq)


def test_handler_create_material(material_h2o_rk_liq):
    handler = MaterialHandler()
    handler.create_flow("inlet", material_h2o_rk_liq.definition)


def test_handler_create_material_despite_port(material_h2o_rk_liq):
    handler = MaterialHandler()
    handler.define_port("inlet")
    with raises(KeyError):
        handler.create_flow("inlet", material_h2o_rk_liq.definition)
