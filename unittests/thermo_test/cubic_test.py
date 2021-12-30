# -*- coding: utf-8 -*-

# stdlib modules
from sys import path
from pathlib import Path

# external modules
from casadi import SX
from pytest import raises

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


if __name__ == "__main__":
    from pytest import main
    # only this file, very verbose and print stdout when started from here.
    main([__file__, "-v", "-v", "-rP"])
