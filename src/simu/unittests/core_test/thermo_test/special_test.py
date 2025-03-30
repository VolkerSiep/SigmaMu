# -*- coding: utf-8 -*-

# internal modules
from simu.core.utilities import assert_reproduction, SymbolQuantity, base_unit
from simu.app.thermo.contributions.special import Derivative


# auxiliary functions
def sym(name: str, units: str) -> SymbolQuantity:
    """Define a scalar symbolic quantity"""
    return SymbolQuantity(name, base_unit(units))


def vec(name: str, size: int, units: str) -> SymbolQuantity:
    """define a vector symbolic quantity"""
    return SymbolQuantity(name, base_unit(units), size)


def test_derivative():
    """Test definition of Derivative contribution"""
    tau = sym("tau", "dimless")
    C = vec("C", 3, "dimless")
    poly = C[0] + (C[1] + C[2] * tau) * tau
    res = {"T": tau, "poly": poly}
    opt = {"y": "poly", "x": "T"}
    deri = Derivative(["A", "B"], opt)
    deri.define(res)
    assert_reproduction(str(res["dpoly_dT"]))
