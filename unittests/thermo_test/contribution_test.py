# -*- coding: utf-8 -*-

# stdlib modules
from sys import path
from pathlib import Path

# external modules
from casadi import SX

# reproductiontest
path.append(str(Path(__file__).absolute().parents[1]))
from reproductiontest import assert_reproduction

def test_H0S0ReferenceState():
    from mushell.thermo import H0S0ReferenceState
    res = {"T": SX.sym('T'), "n": SX.sym('n', 2)}
    par = {"dh_form": {"A": SX.sym('dh_form.A'),
                       "B": SX.sym('dh_form.B')},
           "s_0": {"A": SX.sym('s_0.A'),
                   "B": SX.sym('s_0.B')},
           "T_ref": SX.sym('T_ref'),
           "p_ref": SX.sym('p_ref')}
    cont = H0S0ReferenceState(["A", "B"], {})
    cont.define(res, par)
    result = {i: str(res[i]) for i in "S mu".split()}
    assert_reproduction(result)

def test_LinearHeatCapacity():
    from mushell.thermo import LinearHeatCapacity
    res = {"T": SX.sym('T'), "n": SX.sym('n', 2),
           "S": SX.sym('S_ref'), "mu": SX.sym('mu_ref'),
           "T_ref": SX.sym('T_ref')}
    par = {"cp_a": {"A": SX.sym('cp_a.A'),
                    "B": SX.sym('cp_a.B')},
           "cp_b": {"A": SX.sym('cp_b.A'),
                    "B": SX.sym('cp_b.B')}}
    cont = LinearHeatCapacity(["A", "B"], {})
    cont.define(res, par)
    result = {i: str(res[i]).split(", ") for i in "S mu".split()}
    assert_reproduction(result)

def test_IdealMix():
    from mushell.thermo import IdealMix
    res = {"T": SX.sym('T'), "n": SX.sym('n', 2),
           "S": SX.sym('S_std'), "mu": SX.sym('mu_std')}
    cont = IdealMix(["A", "B"], {})
    cont.define(res, {})
    result = {i: str(res[i]).split(", ") for i in "S mu".split()}
    assert_reproduction(result)

def test_GibbsIdealGas():
    from mushell.thermo import GibbsIdealGas
    res = {"T": SX.sym('T'), "p": SX.sym('p'), "n": SX.sym('n', 2),
           "p_ref": SX.sym('p_ref'),
           "S": SX.sym('S_im'), "mu": SX.sym('mu_im')}
    cont = GibbsIdealGas(["A", "B"], {})
    cont.define(res, {})
    result = {i: str(res[i]).split(", ") for i in "S V mu".split()}
    assert_reproduction(result)


if __name__ == "__main__":
    from pytest import main
    # only this file, very verbose and print stdout when started from here.
    main([__file__, "-v", "-v", "-rP"])
