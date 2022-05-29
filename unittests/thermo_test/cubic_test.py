# -*- coding: utf-8 -*-

# external modules
from casadi import DM, SX, Function, jacobian, vertcat
from pytest import raises, mark

# internal modules
from simu.thermo.cubic.rk import RedlichKwongEOSLiquid, RedlichKwongEOSGas
from simu.constants import R_GAS
from simu.utilities import assert_reproduction, user_agree


def test_critical_parameters():
    from simu.thermo.cubic import CriticalParameters

    res = {}
    par = {"T_c": {"A": SX.sym('T_c.A'),
                   "B": SX.sym('T_c.B')},
           "p_c": {"A": SX.sym('p_c.A'),
                   "B": SX.sym('p_c.B')},
           "omega": {"A": SX.sym('omega.A'),
                      "B": SX.sym('omega.B')}}
    cont = CriticalParameters(["A", "B"], {})
    cont.define(res, par)
    result = {key: str(value) for key, value in res.items()}
    assert_reproduction(result)


def test_linear_mixing_rule():
    from simu.thermo.cubic import LinearMixingRule

    res = {"T": SX.sym('T'), "n": SX.sym('n', 2)}
    par = {"c_i": {"A": SX.sym('c_i.A'),
                   "B": SX.sym('c_i.B')}}
    opt = {"target": "c", "src_mode": LinearMixingRule.PARAMETER}
    cont = LinearMixingRule(["A", "B"], opt)
    cont.define(res, par)
    result = str(res["c"])
    assert_reproduction(result)


def test_RedlichKwongEOS():
    from simu.thermo.cubic.rk import RedlichKwongEOSLiquid
    res = {"T": SX.sym('T'), "V": SX.sym('V'), "n": SX.sym('n', 2),
           "S": SX.sym('S'), "p": SX.sym('p'), "mu": SX.sym('mu', 2)}
    res["ceos_a"] = SX.sym('A0') + res["T"] * SX.sym('dAdT')
    res["ceos_b"] = SX.sym('B0') + res["T"] * SX.sym('dBdT')
    res["ceos_c"] = SX.sym('C0')
    res["state"] = SX.sym('x', 4)
    cont = RedlichKwongEOSLiquid(["A", "B"], {})
    cont.define(res, {})
    result = {i: str(res[i]) for i in "S p mu".split()}
    assert_reproduction(result)


def test_RedlicKwongAbstract():
    from simu.thermo.cubic.rk import RedlichKwongEOS
    with raises(TypeError) as excinfo:
        RedlichKwongEOS(["A", "B"], {})
    assert "abstract" in str(excinfo.value)


def test_NonSymmmetricMixingRule():
    from simu.thermo.cubic import NonSymmetricMixingRule
    res = {"T": SX.sym('T'), "n": SX.sym('n', 3),
           "a_i": SX.sym('a_i', 3)}
    options = {
      "k_1": [["A", "B"], ["A", "C"]],
      "k_2": [["A", "B"]],
      "l_1": [["B", "A"], ["C", "B"]],
      "target": "a"}
    cont = NonSymmetricMixingRule(["A", "B", "C"], options)
    par = {'T_ref': SX.sym("T_ref"),
           'k_1': {'A': {'B': SX.sym("k_1_AB"),
                         'C': SX.sym("k_1_AC")}},
           'k_2': {'A': {'B': SX.sym("k_2_AB")}},
           'l_1': {'B': {'A': SX.sym("k_1_AB")},
                   'C': {'B': SX.sym("k_1_AB")}}}
    cont.define(res, par)
    result = str(res["a"])
    assert_reproduction(result)


def test_redlich_kwong_a_function():
    from simu.thermo.cubic.rk import RedlichKwongAFunction
    res = {"alpha": SX.sym('alpha', 2),
           "T_c": SX.sym('T_c', 2), "p_c": SX.sym('p_c', 2)}
    cont = RedlichKwongAFunction(["A", "B"], {})
    cont.define(res, {})
    result = str(res["ceos_a_i"])
    assert_reproduction(result)


def test_redlich_kwong_b_function():
    from simu.thermo.cubic.rk import RedlichKwongBFunction
    res = {"T_c": SX.sym('T_c', 2), "p_c": SX.sym('p_c', 2)}
    cont = RedlichKwongBFunction(["A", "B"], {})
    cont.define(res, {})
    result = str(res["ceos_b_i"])
    assert_reproduction(result)


def test_rk_m_factor():
    from simu.thermo.cubic.rk import RedlichKwongMFactor
    res = {"omega": SX.sym('w', 2)}
    cont = RedlichKwongMFactor(["A", "B"], {})
    cont.define(res, {})
    result = str(res["m_factor"])
    assert_reproduction(result)


def test_BostonMathiasAlphaFunction():
    from simu.thermo.cubic import BostonMathiasAlphaFunction
    res = {"m_factor": SX.sym('m', 2), "T_c": SX.sym('T_c', 2),
           "T": SX.sym('T')}
    par = {"eta": {"A": SX.sym('eta.A'),
                   "B": SX.sym('eta.B')}}
    cont = BostonMathiasAlphaFunction(["A", "B"], {})
    cont.define(res, par)
    result = str(res["alpha"][0])
    assert_reproduction(result)
    return res, par


def test_BostonMathiasAlphaFunction_smoothness():
    """Check smoothness of alpha function at critical temperature, where
    the expressions switches to the super-critical extrapolation"""
    res, par = test_BostonMathiasAlphaFunction()
    T, T_c, m, alpha = res["T"], res["T_c"], res["m_factor"], res["alpha"]
    eta = vertcat(par["eta"]["A"], par["eta"]["B"])

    dadt = jacobian(alpha, T)
    d2adt2 = jacobian(dadt, T)

    f = Function("alpha", [T, T_c, m, eta], [alpha, dadt, d2adt2])

    def props(eps):
        T, Tc, m, eta = 300, [300, 400], [0.6, 0.6], [0.12, 0.06]
        res = f(T, Tc, m, eta)
        return [float(r[0]) for r in res]

    eps = 1e-10
    a_sub, dadt_sub, d2adt2_sub = props(-eps)  # sub-critical
    a_sup, dadt_sup, d2adt2_sup = props(eps)  # super-critical

    assert abs(a_sub - 1) < 1e-7, "sub-critical alpha unequal unity"
    assert abs(a_sup - 1) < 1e-7, "super-critical alpha unequal unity"
    assert abs(dadt_sub - dadt_sup) < 1e-7, "first derivative not smooth"
    assert abs(d2adt2_sub - d2adt2_sup) < 1e-7, "second derivative not smooth"


def test_initialise_rk():
    T, p, n = 370.0, 1e5, [0.5, 0.5]
    # try to imitate water
    res = {"ceos_a": 15, "ceos_b": 2.5e-5, "ceos_c": 1e-5}
    liq = RedlichKwongEOSLiquid(["A", "B"], {})
    gas = RedlichKwongEOSGas(["A", "B"], {})
    v_liq = liq.initial_state(T, p, n, res)[1]
    v_gas = gas.initial_state(T, p, n, res)[1]
    assert abs(v_liq - 1.526e-5) < 1e-8  # value 1.5... is validated
    assert 0.02 < v_gas < 0.03


@mark.parametrize("Class", [RedlichKwongEOSGas, RedlichKwongEOSLiquid])
def test_initialise_rk2(Class):
    from simu.constants import R_GAS
    # define upstream expected results
    res = {"T": SX.sym('T'), "V": SX.sym('V'), "n": SX.sym('n', 2),
           "ceos_a": SX.sym('A'), "ceos_b": SX.sym('B'), "ceos_c": SX.sym('C'),
           "S": 0, "mu": 0, "p": 0, "state": SX.sym('x', 4)}
    cont = Class(["A", "B"], {})
    cont.define(res, {})

    T, p, n = 370.0, 1e5, [0.5, 0.5]
    res_float = {"T": T, "n": n,
                 "ceos_a": 15, "ceos_b": 2.5e-5, "ceos_c": 1e-5}
    state = cont.initial_state(T, p, n, res_float)

    # is the rest of the state (except volume) reproduced?
    for ref, val in zip([T, None] + n, state):
        assert ref is None or val == ref
    res_float["V"] = state[1]

    # calculate contribution values with initial state to reproduce pressure
    f = Function("f", [res[k] for k in res_float], [res["p"]])
    p2 = f(*res_float.values()) + R_GAS * T / state[1]
    assert abs(p2 - p) < 1


def test_relax_rk():
    T, V, n = 370.0, 1.5260379390336834e-05, [0.5, 0.5]
    res= {"T": T, "V": V, "n": n, "p": 1e5,
          "ceos_a": 15, "ceos_b": 2.5e-5, "ceos_c": 1e-5,
          "VBC": V + 1e-5 - 2.5e-5,
          "dVBC_dx": DM([0, 1, 1e-5 - 2.5e-5, 1e-5 - 2.5e-5]),
          "dp_dV": -1e10, "ddp_dV_dx": DM([0, 1e15, 0, 0]),
          "dp_dx": DM([0, 0, 0, 0])}
    cont = RedlichKwongEOSLiquid(["A", "B"], {})
    beta = cont.relax(res, DM([0.1, -V, 0.1, 0.1]))
    assert beta < 1


def test_relax_rk_advanced():
    # - provoke parameters to have dp/dV limiting
    # - relax only in V-direction
    # - provide p_V and p_V_x[1] approximately correctly
    # - relax from V = 4.6e-5 by dV = +0.4e-5

    T, V, n = 370.0, 4.6e-5, [0.5, 0.5]
    res= {"T": DM(T), "V": DM(V), "n": DM(n), "p": DM(1.1e6),
          "ceos_a": DM(15), "ceos_b": DM(2.5e-5), "ceos_c": DM(1e-5),
          "VBC": DM(V + 1e-5 - 2.5e-5),
          "dVBC_dx": DM([0, 1, 1e-5 - 2.5e-5, 1e-5 - 2.5e-5]),
          "dp_dV": DM(-2.5e11), "ddp_dV_dx": DM([0, 8.5e16, 0, 0]),
          "dp_dx": DM([0, 0, 0, 0])}
    # in above data, p, p_V and p_V_x are set approximately to reproduce
    # taylor approximation in graph. For plotting, p is evaluated exactly.
    if user_agree("Plot pV-plot for relaxation dp/dV < 0"):
        plot_pv(res)

    cont = RedlichKwongEOSLiquid(["A", "B"], {})
    beta = cont.relax(res, DM([0.0, 4e-6, 0.0, 0.0]))
    assert 0.7 < beta < 0.8  # from plot, this is where dp/dV turns positive


# *** auxiliary routines

def plot_pv(res):
    T, V  = res["T"], res["V"]
    def p(V):
        A, B, C = [res[i] for i in "ceos_a ceos_b ceos_c".split()]
        A /= 33.7
        VC = V + C
        return R_GAS * T / (VC - B)  - A / VC / (VC + B)

    # only plot if running this file interactively
    from numpy import linspace
    import pylab

    volumes = linspace(4.5e-5, 5.2e-5, num=100)
    pressures = p(volumes)
    lin_p = p(V) + (volumes - V) * res["dp_dV"]
    quad_p = lin_p + 0.5 * (volumes - V) ** 2 * res["ddp_dV_dx"][1]
    pylab.plot(volumes, pressures.full(), 1)
    pylab.plot(volumes, lin_p.full(), ":")
    pylab.plot(volumes, quad_p.full(), "--")
    pylab.xlim([volumes[0], volumes[-1]])
    pylab.ylim([0, 1.5e6])
    pylab.xlabel("V [m3]")
    pylab.ylabel("p [Pa]")
    pylab.grid()
    pylab.show()


if __name__ == "__main__":
    from pytest import main
    from sys import argv
    # only this file, very verbose and print stdout when started from here.
    main([__file__, "-v", "-v", "-rP"] + argv[1:])
