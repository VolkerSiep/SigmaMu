from simu.app.thermo.contributions.iapws.standard import (
    ReducedStateIAPWS, StandardStateIAPWS, IdealGasIAPWS)
from simu.core.utilities.testing import assert_reproduction
# from simu.core.thermo import InitialState, SpeciesDefinition

from .utils import sym, vec

def test_reduced_state_iapws(species_definitions_ab):
    res = {"T": sym("T", "K"), "V": sym("V", "m^3"), "n": vec("n", 2, "mol"),
           "mw": vec("mw", 2, "kg/mol")}
    cont = ReducedStateIAPWS(species_definitions_ab, {})
    cont.define(res)
    reference = {
        "tau": f"{res["_tau"]:~}",
        "rho": f"{res["_rho"]:~}"
    }
    assert_reproduction(reference)

# TODO:
#   -  test standard state
#   -  test ideal gas state
#   -  parameterize standard state and check enthalpy of steam at low pressure
