# -*- coding: utf-8 -*-
"""Test module for cubic EOS contributions"""

# external modules
from numpy import linspace
from pytest import mark, raises
import pylab

# internal modules
from simu.core.utilities import (
    assert_reproduction, base_unit, ParameterDictionary, SymbolQuantity,
    jacobian, QFunction, Quantity as Q, base_magnitude, sum1, unit_registry)
from simu.core.utilities.constants import R_GAS
from simu.core.thermo import InitialState
from simu.app.thermo.contributions import (
    CriticalParameters, LinearMixingRule, RedlichKwongEOSLiquid,
    RedlichKwongEOSGas, NonSymmetricMixingRule, RedlichKwongAFunction,
    RedlichKwongBFunction, RedlichKwongMFactor, BostonMathiasAlphaFunction,
    VolumeShift)
from simu.app.thermo.contributions.cubic.rk import RedlichKwongEOS


# auxiliary functions
def sym(name: str, units: str) -> SymbolQuantity:
    """Define a scalar symbolic quantity"""
    return SymbolQuantity(name, base_unit(units))


def vec(name: str, size: int, units: str) -> SymbolQuantity:
    "define a vector symbolic quantity"
    return SymbolQuantity(name, base_unit(units), size)


def test_critical_parameters():
    """Test definition of CriticalParameters contribution"""
    res = {}
    par = ParameterDictionary()
    cont = CriticalParameters(["A", "B"], {})
    cont.define(res, par)
    assert_reproduction(res)


def test_volume_shift():
    """Test definition of VolumeShift contribution"""
    res = {}
    par = ParameterDictionary()
    cont = VolumeShift(["A", "B"], {})
    cont.define(res, par)
    assert_reproduction([res, par])


def test_linear_mixing_rule():
    """Test definition of LinearMixingRule contribution"""
    res = {"T": sym("T", "K"), "n": vec("n", 2, "mol"),
           "c_i": vec("c_i", 2, "m**3/mol")}
    par = ParameterDictionary()
    opt = {"target": "c"}
    cont = LinearMixingRule(["A", "B"], opt)
    cont.define(res, par)
    assert_reproduction(res["c"])


def test_redlich_kwong_eos():
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
    cont = RedlichKwongEOSLiquid(["A", "B"], {})
    cont.define(res, {})
    keys = "S p mu _ceos_a_T _ceos_b_T _VBC _dp_dV".split()
    result = {k: res[k] for k in keys}
    assert_reproduction(result)


def test_abstract_class_init():
    """Check that abstract RedlichKwongEOS class cannot be instantiated"""
    with raises(TypeError) as excinfo:
        # pylint: disable=abstract-class-instantiated
        RedlichKwongEOS(["A", "B"], {})
    assert "abstract" in str(excinfo.value)


def test_non_symmmetric_mixing_rule():
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
    cont = NonSymmetricMixingRule(["A", "B", "C"], options)
    par = ParameterDictionary()
    cont.define(res, par)
    assert_reproduction(res["a"])


def test_non_symmmetric_mixing_rule_no_interaction():
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
    cont = NonSymmetricMixingRule(["A", "B", "C"], options)
    par = ParameterDictionary()
    cont.define(res, par)


def test_redlich_kwong_a_function():
    """Test definition of RedlichKwongAFunction contribution"""
    res = {
        "_alpha": vec('alpha', 2, "dimless"),
        "_T_c": vec('T_c', 2, "K"),
        "_p_c": vec('p_c', 2, "bar")
    }
    cont = RedlichKwongAFunction(["A", "B"], {})
    cont.define(res, {})
    assert_reproduction(res["_ceos_a_i"])


def test_redlich_kwong_b_function():
    """Test definition of RedlichKwongBFunction contribution"""
    res = {"_T_c": vec('T_c', 2, "K"), "_p_c": vec('p_c', 2, "bar")}
    cont = RedlichKwongBFunction(["A", "B"], {})
    cont.define(res, {})
    assert_reproduction(res["_ceos_b_i"])


def test_rk_m_factor():
    """Test definition of RedlichKwongMFactor contribution"""
    res = {"_omega": vec('w', 2, "dimless")}
    cont = RedlichKwongMFactor(["A", "B"], {})
    cont.define(res, {})
    assert_reproduction(res["_m_factor"])


def test_boston_mathias_alpha_function():
    """Test definition of BostonMathiasAlphaFunction contribution"""
    res, par = create_boston_mathias_alpha_function()
    assert_reproduction(res["_alpha"][0])


def test_boston_mathias_alpha_function_smoothness():
    """Check smoothness of alpha function at critical temperature, where
    the expressions switches to the super-critical extrapolation"""
    res, par = create_boston_mathias_alpha_function()

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


def test_initialise_rk():
    """Test volume initialisation of RK-model"""
    T, p = Q("100 degC"), Q("1 bar")
    n = Q([0.5, 0.5], "mol")
    # try to imitate water
    res = {
        "_ceos_a": Q("15 Pa * m**6"),
        "_ceos_b": Q("25 ml"),
        "_ceos_c": Q("10 ml")
    }
    liq = RedlichKwongEOSLiquid(["A", "B"], {})
    gas = RedlichKwongEOSGas(["A", "B"], {})
    ini_state = InitialState(temperature=T, pressure=p, mol_vector=n)
    v_liq = liq.initial_state(ini_state, res)[1]
    v_gas = gas.initial_state(ini_state, res)[1]
    assert abs(v_liq - 1.526e-5) < 1e-8  # value 1.5... is validated
    assert 0.02 < v_gas < 0.03


@mark.parametrize("cls", [RedlichKwongEOSGas, RedlichKwongEOSLiquid])
def test_initialise_rk2(cls):
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
    cont = cls(["A", "B"], {})
    cont.define(res, {})  # now the ideal part in res is overwritten

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
    assert base_magnitude(n).tolist() == state[2].tolist()
    res_num["V"] = Q(state[1], "m**3")

    # calculate contribution values with initial state to reproduce pressure
    args = {k: res[k] for k in res_num}
    args.update(ideal)
    func = QFunction(args, {"p_res": res["p"]})

    args = {**res_num, **ideal_num}
    p2 = func(args)["p_res"] + sum(n) * R_GAS * T / args["V"]
    assert abs(p2 - p) < Q("1 Pa")


def test_relax_rk():
    """test relaxation method for RK model"""
    v_magnitude = 1.5260379390336834e-05
    V = Q(v_magnitude, "m**3")
    res = {
        "_state": Q([370.0, v_magnitude, 0.5, 0.5]),
        "V": V,
        "p": Q("1 bar"),
        "_VBC": V + Q("-15 ml"),
        "_dVBC_dx": Q([0, 1e6, 10 - 25, 10 - 25], "ml"),
        "_dp_dV": Q("-0.1 bar/ml"),
        "_ddp_dV_dx": Q([0, 1e4, 0, 0], "bar/ml"),
        "_dp_dx": Q([0, 0, 0, 0], "bar")
    }
    cont = RedlichKwongEOSLiquid(["A", "B"], {})
    delta = [0.1, -v_magnitude, 0.1, 0.1]  # this one is float
    beta = cont.relax(res, delta)
    assert beta < 1


def test_relax_rk_advanced():
    """test relaxation method for RK models, limited by dp/DV:

    - provoke parameters to have dp/dV limiting
    - relax only in V-direction
    - provide p_V and p_V_x[1] approximately correctly
    - relax from V = 4.6e-5 by dV = +0.4e-5"""

    state = [370.0, 4.6e-05, 0.5, 0.5]
    T = Q(state[0], "K")
    V = Q(state[1], "m**3")
    n = Q(state[2:], "mol")
    res = {
        "_state": Q(state),
        "T": T,
        "V": V,
        "n": n,
        "p": Q(11, "bar"),
        "_ceos_a": Q("15 Pa * m**6"),
        "_ceos_b": Q("25 ml"),
        "_ceos_c": Q("10 ml"),
        "_dVBC_dx": Q([0, 1e6, 10 - 25, 10 - 25], "ml"),
        "_dp_dV": Q("-2.5 bar/ml"),
        "_ddp_dV_dx": Q([0, 8.5e5, 0, 0], "bar/ml"),
        "_dp_dx": Q([0, 0, 0, 0], "bar")
    }
    res["_VBC"] = V + res["_ceos_c"] - res["_ceos_b"]
    # in above data, p, p_V and p_V_x are set approximately to reproduce
    # taylor approximation in graph. For plotting, p is evaluated exactly.


    # TODO: make an example folder with these things.
    # if user_agree("Plot pV-plot for relaxation dp/dV < 0"):
    #     plot_pv(res)

    cont = RedlichKwongEOSLiquid(["A", "B"], {})
    beta = cont.relax(res, [0.0, 4e-6, 0.0, 0.0])
    assert 0.7 < beta < 0.8  # from plot, this is where dp/dV turns positive

# *** auxiliary routines

def create_boston_mathias_alpha_function():
    res = {
        "_m_factor": vec("m", 2, "dimless"),
        "_T_c": vec("T_c", 2, "K"),
        "T": sym("T", "K")
    }
    par = ParameterDictionary()
    cont = BostonMathiasAlphaFunction(["A", "B"], {})
    cont.define(res, par)
    return res, par


def plot_pv(res):
    """auxiliary method to plot pv-graph and linear/quadratic approximation"""
    T, V = res["T"], res["V"]

    def p(V):
        A, B, C = [res[i] for i in "_ceos_a _ceos_b _ceos_c".split()]
        A /= 33.7
        VC = V + C
        return sum1(res["n"]) * R_GAS * T / (VC - B) - A / VC / (VC + B)

    # only plot if running this file interactively
    volumes = linspace(45, 52, num=100) * Q("1 ml")
    pressures = p(volumes)
    lin_p = p(V) + (volumes - V) * res["_dp_dV"]
    ddp_dv2 = res["_ddp_dV_dx"][1] * Q("m**-3")
    quad_p = lin_p + (volumes - V)**2 * ddp_dv2 / 2

    unit_registry.setup_matplotlib(True)  # allow plotting with units
    pylab.plot(volumes, pressures.to("bar"), 1)
    pylab.plot(volumes, lin_p.to("bar"), ":")
    pylab.plot(volumes, quad_p.to("bar"), "--")
    pylab.xlim([volumes[0], volumes[-1]])
    # pylab.ylim(Q([0, 15], "bar"))
    pylab.xlabel("V [m3]")
    pylab.ylabel("p [bar]")
    pylab.grid()
    pylab.show()


if __name__ == "__main__":
    from pytest import main
    from sys import argv
    # only this file, very verbose and print stdout when started from here.
    main([__file__, "-v", "-v", "-s", "-rP"] + argv[1:])
