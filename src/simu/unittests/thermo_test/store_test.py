from _pytest.python_api import raises
from pint import DimensionalityError, UndefinedUnitError

from simu import ThermoParameterStore, StringDictThermoSource, SpeciesDefinition
from simu.app.thermo.factories import ThermoStructure, RegThermoFactory
from simu.core.utilities.testing import assert_reproduction


RK_LIQ = "Boston-Mathias-Redlich-Kwong-Liquid"


def test_create_frame(rk_h2o_frame):
    assert_reproduction(rk_h2o_frame.parameter_structure)


def test_get_thermo_properties(rk_h2o_frame):
    store = ThermoParameterStore()
    symbols = store.get_symbols(rk_h2o_frame.parameter_structure)
    assert symbols["CriticalParameters"]["T_c"]["H2O"].units == "kelvin"


def test_get_thermo_properties_twice(rk_h2o_frame):
    store = ThermoParameterStore()
    symbols1 = store.get_symbols(rk_h2o_frame.parameter_structure)
    symbols2 = store.get_symbols(rk_h2o_frame.parameter_structure)
    qty1 = symbols1["CriticalParameters"]["T_c"]["H2O"]
    qty2 = symbols2["CriticalParameters"]["T_c"]["H2O"]
    assert qty1 is qty2


def test_get_thermo_properties_twice_wrong():
    store = ThermoParameterStore()
    struct_1 = {"A": {"B": "K"}}
    struct_2 = {"A": {"B": "m"}}
    store.get_symbols(struct_1)
    with raises(DimensionalityError):
        store.get_symbols(struct_2)


def test_thermo_source_get_item():
    struct = {"T": {"H2O": "100 K", "NH3": "200 K"},
              "p": {"H2O": "10 bar", "NH3": "20 atm"}}
    source = StringDictThermoSource(struct)
    t_h2o = source["T", "H2O"]
    assert f"{t_h2o:~}" == "100 K"


def test_thermo_source_try_nonleaf():
    struct = {"T": {"H2O": {"A": "1 K", "B": "2 K"}, "NH3": "200 K"},
              "p": {"H2O": "10 bar", "NH3": "20 atm"}}
    source = StringDictThermoSource(struct)
    with raises(KeyError):
        _ = source["T", "H2O"]


def test_thermo_source_try_non_quantity():
    struct = {"T": {"H2O": "Bite me!", "NH3": "200 K"},
              "p": {"H2O": "10 bar", "NH3": "20 atm"}}
    with raises(UndefinedUnitError):
        _ = StringDictThermoSource(struct)


def test_get_thermo_property_values_sources():
    store = ThermoParameterStore()
    struct = {"T": {"H2O": "K"}, "p": {"H2O": "bar"}}
    store.get_symbols(struct)

    struct = {"T": {"H2O": "100 K", "NH3": "200 K"},
              "p": {"H2O": "10 bar", "NH3": "20 atm"}}
    store.add_source("Dagbladet", StringDictThermoSource(struct))

    values = store.get_all_values()
    assert f"{values['T']['H2O']:~}" == "100 K"
    sources = store.get_sources()
    assert sources["T"]["H2O"] == "Dagbladet"


def test_get_thermo_missing():
    store = ThermoParameterStore()
    struct = {"T": {"CH4": "K", "H2O": "K"}}
    store.get_symbols(struct)

    struct = {"T": {"H2O": "100 K", "NH3": "200 K"},
              "p": {"H2O": "10 bar", "NH3": "20 atm"}}
    store.add_source("Dagbladet", StringDictThermoSource(struct))

    missing = store.get_missing_symbols()
    assert missing["T"]["CH4"].strip() == "K"


def test_get_missing_dimensionless():
    store = ThermoParameterStore()
    struct = {"T": {"CH4": "", "H2O": ""}}
    store.get_symbols(struct)
    missing = store.get_missing_symbols()
    assert len(missing["T"]) == 2


def test_get_thermo_property_values_two_sources():
    store = ThermoParameterStore()
    struct = {"T": {"H2O": "K"}, "p": {"H2O": "bar"}}
    store.get_symbols(struct)

    struct = {"T": {"H2O": "100 K", "NH3": "200 K"}}
    store.add_source("Dagbladet", StringDictThermoSource(struct))
    struct = {"p": {"H2O": "10 bar", "NH3": "20 atm"}}
    store.add_source("VG", StringDictThermoSource(struct))

    _ = store.get_all_values()
    sources = store.get_sources()
    assert sources["T"]["H2O"] == "Dagbladet"
    assert sources["p"]["H2O"] == "VG"


def test_get_same_thermo_property_values_two_sources():
    store = ThermoParameterStore()
    struct = {"T": {"H2O": "K"}, "p": {"H2O": "bar"}}
    store.get_symbols(struct)

    struct = {"T": {"H2O": "100 K", "NH3": "200 K"},
              "p": {"H2O": "10 bar", "NH3": "20 atm"}}
    store.add_source("Dagbladet", StringDictThermoSource(struct))
    struct = {"p": {"H2O": "10 bar", "NH3": "20 atm"}}
    store.add_source("VG", StringDictThermoSource(struct))

    _ = store.get_all_values()
    sources = store.get_sources()
    assert sources["T"]["H2O"] == "Dagbladet"
    assert sources["p"]["H2O"] == "VG"
