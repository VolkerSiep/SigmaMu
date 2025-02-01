# -*- coding: utf-8 -*-
"""Test module for basic contributions"""

from pytest import raises

# internal modules
from simu.core.thermo import InitialState, SpeciesDefinition
from simu.app.thermo.contributions import (
    GibbsIdealGas, H0S0ReferenceState, HelmholtzIdealGas, IdealMix,
    LinearHeatCapacity, ConstantGibbsVolume, MolecularWeight, ChargeBalance)
from simu.core.utilities import (
    ParameterDictionary, Quantity, SymbolQuantity, assert_reproduction,
    base_unit)


# auxiliary functions
def sym(name, units):
    """Return a scalar symbol of given name and units"""
    return SymbolQuantity(name, base_unit(units))


def vec(name, size, units):
    """Return a vector symbol of given name, units, and size"""
    return SymbolQuantity(name, base_unit(units), size)


def test_h0s0_reference_state(species_definitions_ab):
    """Test definition of H0S0ReferenceState contribution"""

    res = {"T": sym("T", "K"), "n": vec("n", 2, "mol")}
    bounds = {}
    par = ParameterDictionary()

    cont = H0S0ReferenceState(species_definitions_ab, {})
    cont.define(res)
    to_reproduce = {
        "res": {i: str(res[i])
                for i in "S mu".split()},
        "par_names": list(cont.parameters.keys())
    }
    assert_reproduction(to_reproduce)


def test_linear_heat_capacity(species_definitions_ab):
    """Test definition of LinearHeatCapacity contribution"""
    res = {
        "T": sym("T", "K"),
        "T_ref": sym("T_ref", "K"),
        "n": vec("n", 2, "mol"),
        "S": sym("S_ref", "J/K"),
        "mu": vec("mu_ref", 2, "J/mol")
    }
    bounds = {}
    cont = LinearHeatCapacity(species_definitions_ab, {})
    par = ParameterDictionary()
    cont.define(res)
    result = {i: str(res[i]).split(", ") for i in "S mu".split()}
    assert_reproduction(result)


def test_ideal_mix(species_definitions_ab):
    """Test definition of IdealMix contribution"""
    res = {
        "T": sym("T", "K"),
        "n": vec("n", 2, "mol"),
        "S": sym("S_std", "J/K"),
        "mu": vec("mu_std", 2, "J/mol")
    }
    bounds = {}
    cont = IdealMix(species_definitions_ab, {})
    cont.define(res)
    result = {i: str(res[i]).split(", ") for i in "S mu".split()}
    assert_reproduction(result)


def test_gibbs_ideal_gas(species_definitions_ab):
    """Test definition of GibbsIdealGas contribution"""
    res = {
        "T": sym("T", "K"),
        "p": sym("p", "bar"),
        "n": vec('n', 2, "mol"),
        "p_ref": sym("p_ref", "bar"),
        "S": sym("S_im", "J/K"),
        "mu": vec("mu_im", 2, "J/mol")
    }
    cont = GibbsIdealGas(species_definitions_ab, {})
    cont.define(res)
    result = {i: str(res[i]).split(", ") for i in "S V mu".split()}
    assert_reproduction(result)

    # unit correct?
    assert res["S"].is_compatible_with("J/K")
    assert res["V"].is_compatible_with("m**3")
    assert res["mu"].is_compatible_with("J/mol")


def test_helmholtz_ideal_gas(species_definitions_ab):
    """Test definition of GibbsIdealGas contribution"""
    res = {
        "T": sym("T", "K"),
        "V": sym("V", "m ** 3"),
        "n": vec('n', 2, "mol"),
        "p_ref": sym("p_ref", "bar"),
        "S": sym("S_im", "J/K"),
        "mu": vec("mu_im", 2, "J/mol")
    }
    bounds = {}
    cont = HelmholtzIdealGas(species_definitions_ab, {})
    cont.define(res)
    result = {i: str(res[i]).split(", ") for i in "S p mu".split()}
    assert_reproduction(result)


def test_helmholtz_ideal_gas_initialise(species_definitions_ab):
    """Test initialisation via Helmholtz ideal gas contribution"""
    cont = HelmholtzIdealGas(species_definitions_ab, {})
    # normally, we would need to provide numeric quantities as results,
    #  but these are not used for ideal gas initialisation.
    initial_state = InitialState(temperature=Quantity("25 degC"),
                                 pressure=Quantity("1 bar"),
                                 mol_vector=Quantity([1, 1], "mol"))

    state = cont.initial_state(initial_state, {})
    ref_volume = 2 * 8.31446261815324 * (273.15 + 25) / 1e5
    assert abs(state[1] / ref_volume - 1) < 1e-7


def test_constant_gibbs_volume(species_definitions_ab):
    """Test definition of constant gibbs volume contribution"""
    res = {
        "p": sym("p", "Pa"),
        "p_ref": sym("p_ref", "Pa"),
        "n": vec("n", 2, "mol"),
        "mu": vec("mu_std", 2, "J/mol")
    }
    bounds = {}
    cont = ConstantGibbsVolume(species_definitions_ab, {})
    cont.define(res)
    result = {i: str(res[i]).split(", ") for i in "V mu".split()}
    assert_reproduction(result)


def test_molecular_weight(species_definitions_ab):
    """Test definition of molecular weights"""
    res, bounds = {}, {}
    cont = MolecularWeight(species_definitions_ab, {})
    cont.define(res)
    mw = res["mw"].to("g/mol").magnitude
    assert_reproduction(str(mw))


def test_charge_balance():
    species = {"A": SpeciesDefinition("Na:1+"),
               "B": SpeciesDefinition("Cl:1-")}
    res, bounds = {"n": vec("n", 2, "mol")}, {}
    cont = ChargeBalance(species, {})
    cont.define(res)
    res_str = str(cont.residuals["balance"].value.m)
    assert res_str == "(n_0-n_1)"

def test_charge_balance_only_positive():
    species = {"A": SpeciesDefinition("Na:1+"),
               "B": SpeciesDefinition("K:1+")}
    res, bounds = {"n": vec("n", 2, "mol")}, {}
    cont = ChargeBalance(species, {})
    with raises(ValueError):
        cont.define(res)

def test_no_charge_balance(species_definitions_ab):
    res, bounds = {"n": vec("n", 2, "mol")}, {}
    cont = ChargeBalance(species_definitions_ab, {})
    cont.define(res)
    assert not cont.residuals
