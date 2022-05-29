# -*- coding: utf-8 -*-

# external modules
from casadi import SX

# internal modules
from simu.utilities import assert_reproduction


def test_derivative():
    from simu.thermo import Derivative
    T = SX.sym("T")
    C = SX.sym("C", 3)
    poly = C[0] + (C[1] + C[2] * T) * T
    res = {"T": T, "poly": poly}
    opt = {"y": "poly", "x": "T"}
    deri = Derivative(["A", "B"], opt)
    deri.define(res, {})
    assert_reproduction(str(res["dpoly_dT"]))



if __name__ == "__main__":
    from pytest import main
    from sys import argv
    # only this file, very verbose and print stdout when started from here.
    main([__file__, "-v", "-v", "-rP"] + argv[1:])
