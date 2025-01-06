# -*- coding: utf-8 -*-
"""Test module for governing thermo objects"""

# stdlib modules
from sys import argv
from pathlib import Path

# external modules
from pytest import main
from logging import DEBUG
from yaml import safe_load

# internal modules
from simu import ThermoFactory, InitialState, SpeciesDefinition
from simu.core.thermo import ThermoFactory
from simu.app.thermo.contributions import (
    H0S0ReferenceState, LinearHeatCapacity, StandardState, IdealMix,
    HelmholtzIdealGas)
from simu.app.thermo import HelmholtzState
from simu.core.utilities import (
    assert_reproduction, parse_quantities_in_struct, Quantity as Qty)


filename = Path(__file__).resolve().parent / "example_parameters.yml"
with open(filename, encoding="utf-8") as file:
    example_parameters = parse_quantities_in_struct(safe_load(file))


def test_create_thermo_factory():
    """just create a ThermoFactory"""
    ThermoFactory()


def test_register_contributions():
    """Create a ThermoFactory and register some contribtions"""
    create_frame_factory()


def test_create_frame(caplog):
    """create a ThermoFrame object"""
    fac = create_frame_factory()
    config = {
        "species": ["N2", "O2", "Ar", "CO2", "H2O"],
        "state": "HelmholtzState",
        "contributions": [
            "H0S0ReferenceState", "LinearHeatCapacity", "StandardState"
        ],
    }
    species = {f: SpeciesDefinition(f) for f in "N2 O2 Ar CO2 H2O".split()}
    with caplog.at_level(DEBUG, logger="simu"):
        fac.create_frame(species, config)
    msg = "\n".join([r.message for r in caplog.records])
    for contrib in config["contributions"]:
        assert contrib in msg


def test_parameter_structure():
    """Retrieve an (empty) parameter structure from created frame"""
    frame = create_simple_frame()
    assert_reproduction(dict(frame.parameter_structure))


def test_property_structure():
    """Retrieve the names of defined properties from created frame"""
    frame = create_simple_frame()
    assert_reproduction(frame.property_structure)


def test_call_frame_flow():
    """Call a created frame with numerical values"""
    result = call_frame(flow=True)[2]
    result = {k: v for k, v in result.items() if k in {"S", "mu"}}
    assert_reproduction(result)


def test_call_frame_state():
    """Call a created frame with numerical values"""
    result = call_frame(flow=False)[2]
    result = {k: v for k, v in result.items() if k in {"S", "mu"}}
    assert_reproduction(result)


def test_initial_state():
    """Test whether initialisation of a Helmholtz ideal gas contribution
    gives the correct volume"""
    frame = create_simple_frame()
    initial_state = InitialState(temperature=Qty("25 degC"),
                                 pressure=Qty("1 bar"),
                                 mol_vector=Qty([1, 1], "mol"))
    x = frame.initial_state(initial_state, example_parameters)
    assert_reproduction(x[1])


# *** helper functions

def call_frame(flow: bool = True):
    """Call a frame object with a state and return all with the result"""
    frame = create_simple_frame()
    state = Qty([398.15, 0.0448, 1, 1])  # = T, V, *n
    result = frame(state, example_parameters, flow=flow)
    return frame, state, result["props"]


def create_frame_factory():
    """Create a ThermoFactory and register standard state contributions"""
    fac = ThermoFactory()
    fac.register(H0S0ReferenceState, LinearHeatCapacity, StandardState,
                 IdealMix, HelmholtzIdealGas)
    fac.register_state_definition(HelmholtzState)
    return fac


def create_simple_frame():
    """Create a ThermoFrame based on just standard state contributions"""
    fac = create_frame_factory()
    config = {
        "species": ["N2", "O2"],
        "state": "HelmholtzState",
        "contributions": [
            "H0S0ReferenceState", "LinearHeatCapacity", "StandardState",
            "IdealMix", "HelmholtzIdealGas"
        ],
    }
    species = {"N2": SpeciesDefinition("N2"), "O2": SpeciesDefinition("O2")}
    return fac.create_frame(species, config)


if __name__ == "__main__":
    # only this file, very verbose and print stdout when started from here.
    main([__file__, "-v", "-v", "-s", "-rP"] + argv[1:])
