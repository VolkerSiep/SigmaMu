from yaml import safe_load
from pytest import fixture

from simu import (SpeciesDefinition, SymbolQuantity, base_unit, ThermoFactory,
                  parse_quantities_in_struct)
from simu.core.utilities.types import Map
from simu.app import (
    all_contributions, BostonMathiasAlphaFunction, HelmholtzState)
from simu.app.data import DATA_DIR


@fixture(scope="session")
def species_definitions_h2o() -> Map[SpeciesDefinition]:
    """A simple example species definition map with 1 species, H2O"""
    return {"H2O": SpeciesDefinition("H2O")}


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

@fixture(scope="session")
def frame_factory():
    """Create a ThermoFactory and register standard state contributions"""
    fac = ThermoFactory()
    fac.register(*all_contributions)
    fac.register_state_definition(HelmholtzState)
    return fac


@fixture(scope="session")
def simple_frame(frame_factory):
    """Create a ThermoFrame based on just standard state contributions"""
    config = {
        "species": ["N2", "O2"],
        "state": "HelmholtzState",
        "contributions": [
            "H0S0ReferenceState", "LinearHeatCapacity", "StandardState",
            "IdealMix", "HelmholtzIdealGas"
        ],
    }
    species = {"N2": SpeciesDefinition("N2"), "O2": SpeciesDefinition("O2")}
    return frame_factory.create_frame(species, config)


@fixture(scope="session")
def iapws_ideal_gas_model(species_definitions_h2o):
    fac = ThermoFactory()
    fac.register(*all_contributions)
    fac.register_state_definition(HelmholtzState)
    config = {
        "species": ["H2O"],
        "state": "HelmholtzState",
        "contributions": [
            "MolecularWeight", "ReducedStateIAPWS",
            "StandardStateIAPWS", "IdealGasIAPWS"
        ]
    }
    frame = fac.create_frame(species_definitions_h2o, config)
    with open(DATA_DIR / "parameters_iapws.yml") as file:
        params = parse_quantities_in_struct(safe_load(file)["data"])
    params = {k: v for k, v in params.items() if not k.startswith("Residual")}
    return frame, params

@fixture(scope="session")
def iapws_model(species_definitions_h2o):
    fac = ThermoFactory()
    fac.register(*all_contributions)
    fac.register_state_definition(HelmholtzState)
    config = {
        "species": ["H2O"],
        "state": "HelmholtzState",
        "contributions": [
            "MolecularWeight", "ReducedStateIAPWS", "StandardStateIAPWS",
            "IdealGasIAPWS",
            "Residual1IAPWS",
            "Residual2IAPWS",
            "Residual3IAPWS",
            "Residual4IAPWS"
        ]
    }
    frame = fac.create_frame(species_definitions_h2o, config)
    with open(DATA_DIR / "parameters_iapws.yml") as file:
        params = parse_quantities_in_struct(safe_load(file)["data"])
    # just for a test
    params2 = {k: v for k, v in params.items()
              if not (
                k.startswith("Residual1") or
                k.startswith("Residual2") or
                k.startswith("Residual3") or
                k.startswith("Residual4")
    )}
    return frame, params
