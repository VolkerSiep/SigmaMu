# -*- coding: utf-8 -*-
"""Test module for ideal contributions"""

# internal modules
from simu.thermo import (GibbsIdealGas, H0S0ReferenceState, HelmholtzIdealGas,
                         IdealMix, LinearHeatCapacity)
from simu.utilities import (ParameterDictionary, Quantity, SymbolQuantity,
                            assert_reproduction, base_unit)


# auxiliary functions
def sym(name, units):
    return SymbolQuantity(name, base_unit(units))


def vec(name, size, units):
    return SymbolQuantity(name, base_unit(units), size)


def test_h0s0_reference_state():
    """Test definition of H0S0ReferenceState contribution"""

    res = {"T": sym("T", "K"), "n": vec("n", 2, "mol")}
    par = ParameterDictionary()
    cont = H0S0ReferenceState(["A", "B"], {})
    cont.define(res, par)
    to_reproduce = {
        "res": {i: str(res[i])
                for i in "S mu".split()},
        "par_names": list(par.keys())
    }
    assert_reproduction(to_reproduce)


def test_linear_heat_capacity():
    """Test definition of LinearHeatCapacity contribution"""
    res = {
        "T": sym("T", "K"),
        "T_ref": sym("T_ref", "K"),
        "n": vec("n", 2, "mol"),
        "S": sym("S_ref", "J/K"),
        "mu": vec("mu_ref", 2, "J/mol")
    }
    cont = LinearHeatCapacity(["A", "B"], {})
    par = ParameterDictionary()
    cont.define(res, par)
    result = {i: str(res[i]).split(", ") for i in "S mu".split()}
    assert_reproduction(result)


def test_ideal_mix():
    """Test definition of IdealMix contribution"""
    res = {
        "T": sym("T", "K"),
        "n": vec("n", 2, "mol"),
        "S": sym("S_std", "J/K"),
        "mu": vec("mu_std", 2, "J/mol")
    }
    cont = IdealMix(["A", "B"], {})
    cont.define(res, {})
    result = {i: str(res[i]).split(", ") for i in "S mu".split()}
    assert_reproduction(result)


def test_gibbs_ideal_gas():
    """Test definition of GibbsIdealGas contribution"""
    res = {
        "T": sym("T", "K"),
        "p": sym("p", "bar"),
        "n": vec('n', 2, "mol"),
        "p_ref": sym("p_ref", "bar"),
        "S": sym("S_im", "J/K"),
        "mu": vec("mu_im", 2, "J/mol")
    }
    cont = GibbsIdealGas(["A", "B"], {})
    cont.define(res, {})
    result = {i: str(res[i]).split(", ") for i in "S V mu".split()}
    assert_reproduction(result)

    # unit correct?
    assert res["S"].is_compatible_with("J/K")
    assert res["V"].is_compatible_with("m**3")
    assert res["mu"].is_compatible_with("J/mol")


def test_helmholtz_ideal_gas():
    """Test definition of GibbsIdealGas contribution"""
    res = {
        "T": sym("T", "K"),
        "V": sym("V", "m ** 3"),
        "n": vec('n', 2, "mol"),
        "p_ref": sym("p_ref", "bar"),
        "S": sym("S_im", "J/K"),
        "mu": vec("mu_im", 2, "J/mol")
    }
    cont = HelmholtzIdealGas(["A", "B"], {})
    cont.define(res, {})
    result = {i: str(res[i]).split(", ") for i in "S p mu".split()}
    assert_reproduction(result)


def test_helmholtz_ideal_gas_initialise():
    """Test initialisation via Helmholtz ideal gas contribution"""
    cont = HelmholtzIdealGas(["A", "B"], {})
    # normally, we would need to provide numeric quantities as results,
    #  but these are not used for ideal gas initialisation.
    state = cont.initial_state(Quantity("25 degC"), Quantity("1 bar"),
                               Quantity([1, 1], "mol"), {})
    V_ref = 2 * 8.31446261815324 * (273.15 + 25) / 1e5
    assert abs(state[1] / V_ref - 1) < 1e-7