# -*- coding: utf-8 -*-
"""Test module for cubic EOS contributions"""

# external modules
from numpy import linspace
from pytest import mark, raises
from matplotlib import pyplot

# internal modules
from simu.core.utilities import (
    assert_reproduction, base_unit, ParameterDictionary, SymbolQuantity,
    jacobian, QFunction, Quantity as Q, base_magnitude, qsum, unit_registry)
from simu.core.utilities.constants import R_GAS
from simu.core.thermo import InitialState
from simu.app.thermo.contributions import (
    CriticalParameters, LinearMixingRule, RedlichKwongEOSLiquid,
    RedlichKwongEOSGas, NonSymmetricMixingRule, RedlichKwongAFunction,
    RedlichKwongBFunction, RedlichKwongMFactor, VolumeShift)
from simu.app.thermo.contributions.cubic.rk import RedlichKwongEOS


# auxiliary functions
def sym(name: str, units: str) -> SymbolQuantity:
    """Define a scalar symbolic quantity"""
    return SymbolQuantity(name, base_unit(units))


def vec(name: str, size: int, units: str) -> SymbolQuantity:
    "define a vector symbolic quantity"
    return SymbolQuantity(name, base_unit(units), size)


def test_critical_parameters(species_definitions_ab):
    """Test definition of CriticalParameters contribution"""
    res = {}
    par = ParameterDictionary()
    cont = CriticalParameters(species_definitions_ab, {})
    cont.define(res)
    assert_reproduction(res)


def test_volume_shift(species_definitions_ab):
    """Test definition of VolumeShift contribution"""
    res = {}
    par = ParameterDictionary()
    cont = VolumeShift(species_definitions_ab, {})
    cont.define(res)
    assert_reproduction([res, cont.parameters])


def test_linear_mixing_rule(species_definitions_ab):
    """Test definition of LinearMixingRule contribution"""
    res = {"T": sym("T", "K"), "n": vec("n", 2, "mol"),
           "c_i": vec("c_i", 2, "m**3/mol")}
    par = ParameterDictionary()
    opt = {"target": "c"}
    cont = LinearMixingRule(species_definitions_ab, opt)
    cont.define(res)
    assert_reproduction(res["c"])


def test_redlich_kwong_eos(species_definitions_ab):
    """Test definition of RedlichKwongEOS contribution"""
    res = {
        "T": sym("T", "K"),
        "V": sym("V", "m**3"),
        "n": vec("n", 2, "mol"),
        "S": sym('S', "J/K"),
        "p": sym('p', "Pa"),
        "mu": vec('mu', 2, "J/mol")
    }
    res.update({
        "_ceos_a":
        sym("A0", "Pa*m**6") + res["T"] * sym('dAdT', "Pa*m**6/K"),
        "_ceos_b":
        sym("B0", "m**3") + res["T"] * sym('dBdT', "m**3/K"),
        "_ceos_c":
        sym("C0", "m**3"),
        "_state":
        vec("x", 4, "dimless")
    })
    cont = RedlichKwongEOSLiquid(species_definitions_ab, {})
    cont.define(res)
    keys = "S p mu _ceos_a_T _ceos_b_T".split()
    result = {k: res[k] for k in keys}
    assert_reproduction(result)


def test_abstract_class_init(species_definitions_ab):
    """Check that abstract RedlichKwongEOS class cannot be instantiated"""
    with raises(TypeError) as exception_info:
        # pylint: disable=abstract-class-instantiated
        RedlichKwongEOS(species_definitions_ab, {})
    assert "abstract" in str(exception_info.value)


def test_non_symmetric_mixing_rule(species_definitions_abc):
    """Test definition of NonSymmetricMixingRule contribution"""
    res = {
        "T": sym("T", "K"),
        "n": vec("n", 3, "mol"),
        "a_i": vec("a_i", 3, "Pa*m**6/mol")
    }
    options = {
        "k_1": [["A", "B"], ["A", "C"]],
        "k_2": [["A", "B"]],
        "l_1": [["B", "A"], ["C", "B"]],
        "target": "a"
    }
    cont = NonSymmetricMixingRule(species_definitions_abc, options)
    par = ParameterDictionary()
    cont.define(res)
    assert_reproduction(res["a"])


def test_non_symmetric_mixing_rule_no_interaction(species_definitions_abc):
    """Test definition of NonSymmetricMixingRule contribution"""
    res = {
        "T": sym("T", "K"),
        "n": vec("n", 3, "mol"),
        "a_i": vec("a_i", 3, "Pa*m**6/mol")
    }
    options = {
        "k_1": [],
        "k_2": [],
        "l_1": [],
        "target": "a"
    }
    cont = NonSymmetricMixingRule(species_definitions_abc, options)
    par = ParameterDictionary()
    cont.define(res)


def test_redlich_kwong_a_function(species_definitions_ab):
    """Test definition of RedlichKwongAFunction contribution"""
    res = {
        "_alpha": vec('alpha', 2, "dimless"),
        "_T_c": vec('T_c', 2, "K"),
        "_p_c": vec('p_c', 2, "bar")
    }
    cont = RedlichKwongAFunction(species_definitions_ab, {})
    cont.define(res)
    assert_reproduction(res["_ceos_a_i"])


def test_redlich_kwong_b_function(species_definitions_ab):
    """Test definition of RedlichKwongBFunction contribution"""
    res = {"_T_c": vec('T_c', 2, "K"), "_p_c": vec('p_c', 2, "bar")}
    cont = RedlichKwongBFunction(species_definitions_ab, {})
    cont.define(res)
    assert_reproduction(res["_ceos_b_i"])


def test_rk_m_factor(species_definitions_ab):
    """Test definition of RedlichKwongMFactor contribution"""
    res = {"_omega": vec('w', 2, "dimless")}
    cont = RedlichKwongMFactor(species_definitions_ab, {})
    cont.define(res)
    assert_reproduction(res["_m_factor"])


def test_boston_mathias_alpha_function(boston_mathias_alpha_function):
    """Test definition of BostonMathiasAlphaFunction contribution"""
    res, par = boston_mathias_alpha_function
    assert_reproduction(res["_alpha"][0])


def test_boston_mathias_alpha_function_smoothness(boston_mathias_alpha_function):
    """Check smoothness of alpha function at critical temperature, where
    the expression switches to the super-critical extrapolation"""
    res, par = boston_mathias_alpha_function

    args = {k: res[k] for k in ["T", "_T_c", "_m_factor"]}
    args["_eta"] = par.get_vector_quantity("eta")

    result = {"_alpha": res["_alpha"]}
    result["_a_t"] = jacobian(result["_alpha"], args["T"])
    result["_a_tt"] = jacobian(result["_a_t"], args["T"])

    f = QFunction(args, result)

    # now the numbers
    args = {
        "_T_c": Q([300, 400], "K"),
        "_m_factor": Q([0.6, 0.6]),
        "_eta": Q([0.12, 0.06])
    }

    def props(eps):
        args["T"] = args["_T_c"][0] + eps
        res = f(args)
        return {k: v[0] for k, v in res.items()}

    eps = Q(1e-10, "K")
    sub = props(-1.0 * eps)  # sub-critical
    sup = props(eps)  # super-critical

    assert abs(sub["_alpha"] - 1) < 1e-7, "sub-critical alpha unequal unity"
    assert abs(sup["_alpha"] - 1) < 1e-7, "super-critical alpha unequal unity"
    res = abs(sup["_a_t"] / sub["_a_t"] - 1)
    assert res < 1e-7, "first derivative not smooth"
    res = abs(sup["_a_tt"] / sub["_a_tt"] - 1)
    assert res < 1e-7, "second derivative not smooth"


def test_initialise_rk(species_definitions_ab):
    """Test volume initialisation of RK-model"""
    T, p = Q("100 degC"), Q("1 bar")
    n = Q([0.5, 0.5], "mol")
    # try to imitate water
    res = {
        "_ceos_a": Q("15 Pa * m**6"),
        "_ceos_b": Q("25 ml"),
        "_ceos_c": Q("10 ml")
    }
    liq = RedlichKwongEOSLiquid(species_definitions_ab, {})
    gas = RedlichKwongEOSGas(species_definitions_ab, {})
    ini_state = InitialState(temperature=T, pressure=p, mol_vector=n)
    v_liq = liq.initial_state(ini_state, res)[1]
    v_gas = gas.initial_state(ini_state, res)[1]
    assert abs(v_liq - 1.526e-5) < 1e-8  # value 1.5... is validated
    assert 0.02 < v_gas < 0.03


@mark.parametrize("cls", [RedlichKwongEOSGas, RedlichKwongEOSLiquid])
def test_initialise_rk2(cls, species_definitions_ab):
    """test initialisation of rk gas and liquid"""

    # define upstream expected results
    res = {
        "T": sym("T", "K"),
        "V": sym("V", "m**3"),
        "n": vec("n", 2, "mol"),
        "_ceos_a": sym("A", "bar*m**6"),
        "_ceos_b": sym("B", "m**3"),
        "_ceos_c": sym("C", "m**3"),
        "_state": vec('x', 4, "dimless")
    }
    ideal = {
        "S": sym("S", "J/K"),
        "mu": vec("mu", 2, "J/mol"),
        "p": sym("p", "Pa")
    }

    res.update(ideal)
    cont = cls(species_definitions_ab, {})
    cont.define(res)  # now the ideal part in res is overwritten

    # now with number quantities
    T, p, n = Q("100 degC"), Q("1 bar"), Q([0.5, 0.5], "mol")
    # try to imitate water
    res_num = {
        "T": T,
        "n": n,
        "_ceos_a": Q("15 Pa * m**6"),
        "_ceos_b": Q("25 ml"),
        "_ceos_c": Q("10 ml"),
    }
    ideal_num = {"S": Q("0 J/K"), "p": Q("0 Pa"), "mu": Q([0, 0], "J/mol")}
    ini_state = InitialState(temperature=T, pressure=p, mol_vector=n)
    state = cont.initial_state(ini_state, res_num)

    # is the rest of the state (except volume) reproduced?
    assert base_magnitude(T) == state[0]
    assert base_magnitude(n).tolist() == state[2:]
    res_num["V"] = Q(state[1], "m**3")

    # calculate contribution values with initial state to reproduce pressure
    args = {k: res[k] for k in res_num}
    args.update(ideal)
    func = QFunction(args, {"p_res": res["p"]})

    args = {**res_num, **ideal_num}
    p2 = func(args)["p_res"] + sum(n) * R_GAS * T / args["V"]
    assert abs(p2 - p) < Q("1 Pa")


# *** auxiliary routines (TODO: plot should be an example instead!)

def plot_pv(res):
    """auxiliary method to plot pv-graph and linear/quadratic approximation"""
    T, V = res["T"], res["V"]

    def p(V):
        A, B, C = [res[i] for i in "_ceos_a _ceos_b _ceos_c".split()]
        A /= 33.7
        VC = V + C
        return qsum(res["n"]) * R_GAS * T / (VC - B) - A / VC / (VC + B)

    # only plot if running this file interactively
    volumes = linspace(45, 52, num=100) * Q("1 ml")
    pressures = p(volumes)
    lin_p = p(V) + (volumes - V) * res["_dp_dV"]
    ddp_dv2 = res["_ddp_dV_dx"][1] * Q("m**-3")
    quad_p = lin_p + (volumes - V)**2 * ddp_dv2 / 2

    unit_registry.setup_matplotlib(True)  # allow plotting with units
    pyplot.plot(volumes, pressures.to("bar"), 1)
    pyplot.plot(volumes, lin_p.to("bar"), ":")
    pyplot.plot(volumes, quad_p.to("bar"), "--")
    pyplot.xlim([volumes[0], volumes[-1]])
    # pyplot.ylim(Q([0, 15], "bar"))
    pyplot.xlabel("V [m3]")
    pyplot.ylabel("p [bar]")
    pyplot.grid()
    pyplot.show()


if __name__ == "__main__":
    from pytest import main
    from sys import argv
    # only this file, very verbose and print stdout when started from here.
    main([__file__, "-v", "-v", "-s", "-rP"] + argv[1:])
