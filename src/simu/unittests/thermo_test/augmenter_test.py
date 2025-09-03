from pprint import pprint

from simu import SpeciesDefinition
from simu.app.thermo.contributions.augmenters.general import (
    GenericProperties, Elemental)
from simu.core.utilities.testing import assert_reproduction

from .utils import sym, vec

def test_generic_properties(species_definitions_ab):
    res = {"T": sym("T", "K"), "p": sym("p", "Pa"), "n": vec("n", 2, "mol"),
           "S": sym("S", "J/K"), "V": sym("V", "m**3"),
           "mu": vec("mu", 2, "J/mol"), "mw": vec("mw", 2, "kg/mol")}
    cont = GenericProperties(species_definitions_ab, {})
    cont.define(res)
    props = {n: f"{res[n]:~}" for n in GenericProperties.provides}
    assert_reproduction(props)


def test_elemental():
    species = {n: SpeciesDefinition(n) for n in ["Na2SO4", "NaHSO4", "Na2S2O3"]}
    res = {"n": vec("n", 3, "mol"), "mw": vec("mw", 3, "kg/mol")}
    cont = Elemental(species, {})
    cont.define(res)

    # all vectors same?
    for key, elem in cont.vectors.items():
        assert elem == ["H", "Na", "O", "S"], key
    assert cont.vectors.keys() == {"m_e", "n_e", "w_e", "x_e"}

    # reproduce calculation
    props = {n: f"{res[n]:~}" for n in Elemental.provides}
    assert_reproduction(props)
