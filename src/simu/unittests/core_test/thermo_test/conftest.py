from pytest import fixture
from simu import (SpeciesDefinition, SymbolQuantity, base_unit,
                  ParameterDictionary)
from simu.app import BostonMathiasAlphaFunction
from simu.core.utilities.types import Map


@fixture(scope="session")
def species_definitions_ab() -> Map[SpeciesDefinition]:
    """A simple example species definition map with 2 species, A and B"""
    return {"A": SpeciesDefinition("N2"),
            "B": SpeciesDefinition("O2")}


@fixture(scope="session")
def species_definitions_abc() -> Map[SpeciesDefinition]:
    """A simple example species definition map with 3 species, A, B and C"""
    return {"A": SpeciesDefinition("N2"),
            "B": SpeciesDefinition("O2"),
            "C": SpeciesDefinition("Ar")}


@fixture(scope="session")
def boston_mathias_alpha_function(species_definitions_ab):
    def sym(name: str, units: str) -> SymbolQuantity:
        return SymbolQuantity(name, base_unit(units))

    def vec(name: str, size: int, units: str) -> SymbolQuantity:
        return SymbolQuantity(name, base_unit(units), size)

    res = {
        "_m_factor": vec("m", 2, "dimless"),
        "_T_c": vec("T_c", 2, "K"),
        "T": sym("T", "K")
    }
    cont = BostonMathiasAlphaFunction(species_definitions_ab, {})
    cont.define(res)
    return res, cont.parameters


