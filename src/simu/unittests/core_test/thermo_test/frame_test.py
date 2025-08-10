"""Test module for governing thermo objects"""

# stdlib modules
from sys import argv
from pathlib import Path

# external modules
from pytest import main
from logging import DEBUG
from yaml import safe_load

# internal modules
from simu import (
    ThermoFactory, InitialState, SpeciesDefinition,
    parse_quantities_in_struct, Quantity as Q)
from simu.core.utilities.testing import assert_reproduction


filename = Path(__file__).resolve().parent / "example_parameters.yml"
with open(filename, encoding="utf-8") as file:
    example_parameters = parse_quantities_in_struct(safe_load(file))


def test_create_thermo_factory():
    """just create a ThermoFactory"""
    ThermoFactory()


def test_register_contributions(frame_factory):
    """Create a ThermoFactory and register some contributions"""
    pass  # just testing the fixture


def test_create_frame(caplog, frame_factory):
    """create a ThermoFrame object"""
    config = {
        "species": ["N2", "O2", "Ar", "CO2", "H2O"],
        "state": "HelmholtzState",
        "contributions": [
            "H0S0ReferenceState", "LinearHeatCapacity", "StandardState"
        ],
    }
    species = {f: SpeciesDefinition(f) for f in "N2 O2 Ar CO2 H2O".split()}
    with caplog.at_level(DEBUG, logger="simu"):
        frame_factory.create_frame(species, config)
    msg = "\n".join([r.message for r in caplog.records])
    for contrib in config["contributions"]:
        assert contrib in msg


def test_parameter_structure(simple_frame):
    """Retrieve an (empty) parameter structure from created frame"""
    assert_reproduction(dict(simple_frame.parameter_structure))


def test_property_structure(simple_frame):
    """Retrieve the names of defined properties from created frame"""
    assert_reproduction(simple_frame.property_structure)


def test_call_frame_flow(simple_frame):
    """Call a created frame with numerical values"""
    result = call_frame(simple_frame, flow=True)[2]
    result = {k: v for k, v in result.items() if k in {"S", "mu"}}
    assert_reproduction(result)


def test_call_frame_state(simple_frame):
    """Call a created frame with numerical values"""
    result = call_frame(simple_frame, flow=False)[2]
    result = {k: v for k, v in result.items() if k in {"S", "mu"}}
    assert_reproduction(result)


def test_initial_state(simple_frame):
    """Test whether initialisation of a Helmholtz ideal gas contribution
    gives the correct volume"""
    initial_state = InitialState(temperature=Q("25 degC"),
                                 pressure=Q("1 bar"),
                                 mol_vector=Q([1, 1], "mol"))
    x = simple_frame.initial_state(initial_state, example_parameters)
    assert_reproduction(x[1])


# *** helper functions

def call_frame(frame, flow: bool = True):
    """Call a frame object with a state and return all with the result"""
    state = Q([398.15, 0.0448, 1, 1])  # = T, V, *n
    result = frame(state, example_parameters, flow=flow)
    return frame, state, result["props"]


if __name__ == "__main__":
    # only this file, very verbose and print stdout when started from here.
    main([__file__, "-v", "-v", "-s", "-rP"] + argv[1:])
