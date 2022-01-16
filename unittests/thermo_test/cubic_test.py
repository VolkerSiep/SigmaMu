# -*- coding: utf-8 -*-

# stdlib modules
from sys import path
from pathlib import Path

# external modules
from casadi import DM, SX, Function, jacobian, vertcat
from pytest import raises, mark

# internal modules
from mushell.thermo.cubic.rk import RedlichKwongEOSLiquid, RedlichKwongEOSGas

# reproductiontest
path.append(str(Path(__file__).absolute().parents[1]))
from reproductiontest import assert_reproduction

def test_critical_parameters():
    from mushell.thermo.cubic import CriticalParameters

    res = {}
    par = {"T_C": {"A": SX.sym('T_C.A'),
                   "B": SX.sym('T_C.B')},
           "P_C": {"A": SX.sym('P_C.A'),
                   "B": SX.sym('P_C.B')},
           "OMEGA": {"A": SX.sym('OMEGA.A'),
                      "B": SX.sym('OMEGA.B')}}
    cont = CriticalParameters(["A", "B"], {})
    cont.define(res, par)
    result = {key: str(value) for key, value in res.items()}
    assert_reproduction(result)


def test_linear_peneloux_volume_shift():
    from mushell.thermo.cubic import LinearPenelouxVolumeShift

    res = {"T": SX.sym('T'), "n": SX.sym('n', 2)}
    par = {"c_i": {"A": SX.sym('c_i.A'),
                   "B": SX.sym('c_i.B')}}
    cont = LinearPenelouxVolumeShift(["A", "B"], {})
    cont.define(res, par)
    result = str(res["CEOS_C"])
    assert_reproduction(result)


def test_RedlichKwongEOS():
    from mushell.thermo.cubic.rk import RedlichKwongEOSLiquid
    res = {"T": SX.sym('T'), "V": SX.sym('V'), "n": SX.sym('n', 2),
           "S": SX.sym('S'), "p": SX.sym('p'), "mu": SX.sym('mu', 2)}
    res["RK_A"] = SX.sym('A0') + res["T"] * SX.sym('dAdT')
    res["RK_B"] = SX.sym('B0') + res["T"] * SX.sym('dBdT')
    res["CEOS_C"] = SX.sym('C0')
    res["state"] = SX.sym('x', 4)
    cont = RedlichKwongEOSLiquid(["A", "B"], {})
    cont.define(res, {})
    result = {i: str(res[i]) for i in "S p mu".split()}
    assert_reproduction(result)


def test_RedlicKwongAbstract():
    from mushell.thermo.cubic.rk import RedlichKwongEOS
    with raises(TypeError) as excinfo:
        RedlichKwongEOS(["A", "B"], {})
    assert "abstract" in str(excinfo.value)


def test_NonSymmmetricMixingRule():
    from mushell.thermo.cubic import NonSymmmetricMixingRule
    res = {"T": SX.sym('T'), "n": SX.sym('n', 3),
           "RK_A_I": SX.sym('a_i', 3)}
    options = {
      "k_1": [["A", "B"], ["A", "C"]],
      "k_2": [["A", "B"]],
      "l_1": [["B", "A"], ["C", "B"]]}
    cont = NonSymmmetricMixingRule(["A", "B", "C"], options)
    par = {'T_ref': SX.sym("T_ref"),
           'k_1': {'A': {'B': SX.sym("k_1_AB"),
                         'C': SX.sym("k_1_AC")}},
           'k_2': {'A': {'B': SX.sym("k_2_AB")}},
           'l_1': {'B': {'A': SX.sym("k_1_AB")},
                   'C': {'B': SX.sym("k_1_AB")}}}
    cont.define(res, par)
    result = str(res["RK_A"])
    assert_reproduction(result)


def test_redlich_kwong_a_function():
    from mushell.thermo.cubic.rk import RedlichKwongAFunction
    res = {"ALPHA_I": SX.sym('alpha', 2),
           "T_C": SX.sym('T_c', 2), "P_C": SX.sym('p_c', 2)}
    cont = RedlichKwongAFunction(["A", "B"], {})
    cont.define(res, {})
    result = str(res["RK_A_I"])
    assert_reproduction(result)


def test_redlich_kwong_b_function():
    from mushell.thermo.cubic.rk import RedlichKwongBFunction
    res = {"T_C": SX.sym('T_c', 2), "P_C": SX.sym('p_c', 2)}
    cont = RedlichKwongBFunction(["A", "B"], {})
    cont.define(res, {})
    result = str(res["RK_B_I"])
    assert_reproduction(result)


def test_rk_m_factor():
    from mushell.thermo.cubic.rk import RedlichKwongMFactor
    res = {"OMEGA": SX.sym('w', 2)}
    cont = RedlichKwongMFactor(["A", "B"], {})
    cont.define(res, {})
    result = str(res["MFAC"])
    assert_reproduction(result)


def test_BostonMathiasAlphaFunction():
    from mushell.thermo.cubic import BostonMathiasAlphaFunction
    res = {"MFAC": SX.sym('m', 2), "T_C": SX.sym('T_c', 2),
           "T": SX.sym('T')}
    par = {"ETA": {"A": SX.sym('ETA.A'),
                   "B": SX.sym('ETA.B')}}
    cont = BostonMathiasAlphaFunction(["A", "B"], {})
    cont.define(res, par)
    result = str(res["ALPHA"][0])
    assert_reproduction(result)
    return res, par


def test_BostonMathiasAlphaFunction_smoothness():
    """Check smoothness of alpha function at critical temperature, where
    the expressions switches to the super-critical extrapolation"""
    res, par = test_BostonMathiasAlphaFunction()
    T, T_c, m, alpha = res["T"], res["T_C"], res["MFAC"], res["ALPHA"]
    eta = vertcat(par["ETA"]["A"], par["ETA"]["B"])

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
    res = {"RK_A": 15, "RK_B": 2.5e-5, "CEOS_C": 1e-5}
    liq = RedlichKwongEOSLiquid(["A", "B"], {})
    gas = RedlichKwongEOSGas(["A", "B"], {})
    v_liq = liq.initial_state(T, p, n, res)[1]
    v_gas = gas.initial_state(T, p, n, res)[1]
    assert v_liq < 2e-5
    assert 0.02 < v_gas < 0.03


@mark.parametrize("Class", [RedlichKwongEOSGas, RedlichKwongEOSLiquid])
def test_initialise_rk2(Class):
    from mushell.constants import R_GAS
    # define upstream expected results
    res = {"T": SX.sym('T'), "V": SX.sym('V'), "n": SX.sym('n', 2),
           "RK_A": SX.sym('A'), "RK_B": SX.sym('B'), "CEOS_C": SX.sym('C'),
           "S": 0, "mu": 0, "p": 0, "state": SX.sym('x', 4)}
    cont = Class(["A", "B"], {})
    cont.define(res, {})

    T, p, n = 370.0, 1e5, [0.5, 0.5]
    res_float = {"T": T, "n": n, "RK_A": 15, "RK_B": 2.5e-5, "CEOS_C": 1e-5}
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
    res= {"T": T, "V": V, "n": n,
          "RK_A": 15, "RK_B": 2.5e-5, "CEOS_C": 1e-5,
          "VBC": V + 1e-5 - 2.5e-5,
          "VBC_x": DM([0, 1, 1e-5 - 2.5e-5, 1e-5 - 2.5e-5])}
    cont = RedlichKwongEOSLiquid(["A", "B"], {})
    beta = cont.relax(res, DM([0.1, -V, 0.1, 0.1]))  # check V - B + C > 0
    assert beta < 1
    print(beta)
    # TODO:
    #  - plot p(V) and check plauslibility
    #  - Also implement that V > 0
    #  - then dp/dV < 0
    #  - check all that with the graph (plot log(p) over log(V))

if __name__ == "__main__":
    from pytest import main
    # only this file, very verbose and print stdout when started from here.
    main([__file__, "-v", "-v", "-rP"])
