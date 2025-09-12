from casadi import DM

from simu import Quantity
from simu.app.thermo.contributions.electrolytes.pitzer import (
    ElectrolyteBasics, PitzerDebyeHueckel, PitzerBinaryInteraction)
from simu.core.utilities.testing import assert_reproduction

from .utils import vec, sym

def test_basics(species_definitions_elec):
    res = {"n": vec("n", 3, "mol")}
    cont = ElectrolyteBasics(species_definitions_elec)
    cont.define(res)
    del res["n"]
    res = {k: f"{v:~}" for k, v in res.items()}
    assert_reproduction(res)

def test_pdh(species_definitions_elec):
    res = {"T": sym("T", "K"), "n": vec("n", 3, "mol"),
           "I": sym("I", "dimless"), "_delta_i_s": vec("d_is", 3, "dimless"),
           "charge": vec("c", 3, "e/mol"), "m_solvent": sym("m_s", "kg"),
           "mw_solvent": sym("M_s", "g/mol"),
           "molality": vec("molality", 3, "mol/kg")}
    inp_keys = set(res.keys())
    res.update(mu=vec("mu", 3, "kJ/mol"), S=sym("S", "J/K"))
    cont = PitzerDebyeHueckel(species_definitions_elec)
    cont.define(res)
    res = {k: f"{v:~}" for k, v in res.items() if not k in inp_keys}
    assert_reproduction(res)

def test_pitzer_binary(species_definitions_elec):
    d_si = DM.zeros(3)
    d_si[0] = 1
    res = {"T": sym("T", "K"), "n": vec("n", 3, "mol"),
           "I": sym("I", "dimless"), "_delta_i_s": Quantity(d_si, "dimless"),
           "charge": vec("c", 3, "e/mol"), "m_solvent": sym("m_s", "kg"),
           "mw_solvent": sym("M_s", "g/mol"),
           "molality": vec("m", 3, "mol/kg")}
    inp_keys = set(res.keys())
    res.update(mu=vec("mu", 3, "kJ/mol"), S=sym("S", "J/K"))
    opts = {f"beta_{k}{m + 1}": [["Na+", "SO42-"]]
            for k in (0, 1) for m in range(5)}
    cont = PitzerBinaryInteraction(species_definitions_elec, opts)
    cont.define(res)
    res = {k: f"{v:~}" for k, v in res.items() if not k in inp_keys}
    assert_reproduction(res)