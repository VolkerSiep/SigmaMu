from casadi import DM

from simu import Quantity, SymbolQuantity
from simu.app.thermo.contributions.electrolytes.pitzer import (
    ElectrolyteBasics, PitzerDebyeHueckel, PitzerBinaryInteraction,
    ExcessBasePitzer, PitzerTernaryInteraction)
from simu.core.utilities.testing import assert_reproduction

from .utils import vec, sym


def test_basics(species_definitions_electrolyte):
    res = {"n": vec("n", 3, "mol")}
    cont = ElectrolyteBasics(species_definitions_electrolyte)
    cont.define(res)
    del res["n"]
    res = {k: f"{v:~}" for k, v in res.items()}
    assert_reproduction(res)


def test_excess_base_pitzer(species_definitions_electrolyte,
                            res_input_electrolyte):
    class SimpleChi(ExcessBasePitzer):
        def define_chi(self, res):
            return {
                "chi": SymbolQuantity("chi", "dimless"),
                "chi_t": SymbolQuantity("chi_t", "1/K"),
                "chi_i": SymbolQuantity("chi_i", "dimless"),
                "chi_m": SymbolQuantity("chi_m", "dimless", self.species)}

    res, inp_keys = res_input_electrolyte
    cont = SimpleChi(species_definitions_electrolyte)
    cont.define(res)
    res = {k: f"{v:~}" for k, v in res.items() if not k in inp_keys}
    assert_reproduction(res)


def test_pdh(species_definitions_electrolyte, res_input_electrolyte):
    res, inp_keys = res_input_electrolyte
    cont = PitzerDebyeHueckel(species_definitions_electrolyte)
    cont.define(res)
    res = {k: f"{v:~}" for k, v in res.items() if not k in inp_keys}
    assert_reproduction(res)


def test_pitzer_binary(species_definitions_electrolyte, res_input_electrolyte):
    res, inp_keys = res_input_electrolyte
    opts = {f"beta_{k}{m + 1}": [["Na+", "SO42-"]]
            for k in (0, 1) for m in range(5)}
    cont = PitzerBinaryInteraction(species_definitions_electrolyte, opts)
    cont.define(res)
    res = {k: f"{v:~}" for k, v in res.items() if not k in inp_keys}
    assert_reproduction(res)


def test_pitzer_ternary(species_definitions_electrolyte, res_input_electrolyte):
    res, inp_keys = res_input_electrolyte
    # in real, the solute should not be part of the interaction.
    opts = {f"gamma_{m + 1}": [["H2O", "Na+", "SO42-"]] for m in range(5)}
    cont = PitzerTernaryInteraction(species_definitions_electrolyte, opts)
    cont.define(res)
    res = {k: f"{v:~}" for k, v in res.items() if not k in inp_keys}
    assert_reproduction(res)
