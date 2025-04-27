from numpy import linspace
from numpy.testing import assert_allclose

from simu import R_GAS
from simu.app.thermo.contributions.iapws.standard import (
    ReducedStateIAPWS, StandardStateIAPWS, IdealGasIAPWS)
from simu.core.utilities.testing import assert_reproduction
from simu.core.utilities.quantity import jacobian, QFunction, Quantity

from .utils import sym, vec

def test_reduced_state_iapws(species_definitions_ab):
    res = {"T": sym("T", "K"), "V": sym("V", "m^3"), "n": vec("n", 2, "mol"),
           "mw": vec("mw", 2, "kg/mol")}
    cont = ReducedStateIAPWS(species_definitions_ab)
    cont.define(res)
    reference = {
        "tau": f"{res["_tau"]:~}",
        "rho": f"{res["_rho"]:~}"
    }
    assert_reproduction(reference)

def test_standard_state_iapws(species_definitions_h2o):
    res = {"T": sym("T", "K"), "n": vec("n", 1, "mol"),
           "_tau": sym("tau", "")}
    cont = StandardStateIAPWS(species_definitions_h2o)
    cont.define(res)

    assert f"{res["mu"].units:~}" == "J / mol"
    assert f"{res["S"].units:~}" == "J / K"

def test_standard_state_iapws_deri(species_definitions_h2o):
    res_0 = {"T": sym("T", "K"), "V": sym("V", "m^3"), "n": vec("n", 1, "mol"),
             "mw": vec("mw", 1, "kg/mol")}
    res = dict(res_0)

    r_state = ReducedStateIAPWS(species_definitions_h2o)
    s_state = StandardStateIAPWS(species_definitions_h2o)
    r_state.define(res)
    s_state.define(res)
    param = r_state.parameters | s_state.parameters
    residual = jacobian(res["mu"], res["T"]) + jacobian(res["S"], res["n"])
    f = QFunction({"param": param, "state": res_0}, {"r": residual})

    # I put some T, V, n and then permutate unity for all n_i to check equality
    args = {
        "param": {
            "rho_c": {"H2O": Quantity("322 kg/m**3")},
            "T_c": {"H2O": Quantity("647.096 K")}},
        "state": {
            "T": Quantity("300.0 K"),
            "V": Quantity("0.024 m**3"),
            "n": Quantity([2.0], "mol"),
            "mw": Quantity([0.018], "kg/mol")}
    }

    def evaluate(k):
        for i in range(1, 9):
            p = 1 if i == k else 0
            args["param"][f"n_{i}"] = {"H2O": Quantity(p)}
            if i > 3:
                args["param"][f"g_{i}"] = {"H2O": Quantity(2.345)}
        return float(f(args)["r"].magnitude)

    result = [evaluate(k) for k in range(1, 9)]
    assert_allclose(result, [0.0] * len(result), atol=1e-13)


def test_ideal_gas_iapws(species_definitions_ab):
    res = {"T": sym("T", "K"), "V": sym("V", "m**3"), "n": vec("n", 2, "mol"),
           "_rho": vec("rho", 2, ""), "mu": vec("mu0", 2, "J/mol"),
           "S": sym("S0", "J/K")}
    cont = IdealGasIAPWS(species_definitions_ab)
    cont.define(res)
    reference = {
        "mu": f"{res["mu"]:~}",
        "S": f"{res["S"]:~}",
        "p": f"{res["p"]:~}"
    }
    assert_reproduction(reference)


def test_iapws_call_ideal_gas(iapws_ideal_gas_model):
    frame, param = iapws_ideal_gas_model
    state = [373.15, 1e-2, 1.0]
    res = frame(state, param)["props"]
    z = res["p"] * res["V"] / (R_GAS * res["T"] * res["n"])
    assert z.to("dimensionless").magnitude == 1.0


def test_iapws_ideal_gas_enthalpy(iapws_ideal_gas_model):
    def enthalpy(temp):
        state = [temp, 100.0, 1.0]  # very low pressure
        res = frame(state, param)["props"]
        return res["mu"] * res["n"] + res["T"] * res["S"]

    frame, param = iapws_ideal_gas_model
    t0 = 273.15
    temperatures = linspace(20, 200, num=10)
    h = [float(enthalpy(t + t0).to("kJ").magnitude) for t in temperatures]
    h_ref = [45.734, 46.406, 47.081, 47.758, 48.437, 49.120, 49.805, 50.495,
             51.188, 51.885]  # from NIST at low pressure
    assert_allclose(h, h_ref, rtol=1e-4)


def test_iapws_ideal_gas_entropy(iapws_ideal_gas_model):
    def entropy(temp):
        state = [temp, 100.0, 1.0]  # very low pressure
        res = frame(state, param)["props"]
        return res["S"]

    frame, param = iapws_ideal_gas_model
    t0 = 273.15
    temperatures = linspace(20, 200, num=10)
    s = [float(entropy(t + t0).to("J/K").magnitude) for t in temperatures]
    s_ref = [194.13, 195.80, 197.37, 198.86, 200.28, 201.62, 202.91, 204.15,
             205.34, 206.48]  # from NIST at same volume
    assert_allclose(s, s_ref, rtol=1e-4)
