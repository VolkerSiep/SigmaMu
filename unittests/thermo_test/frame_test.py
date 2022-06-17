# -*- coding: utf-8 -*-

"""Test module for governing thermo objects"""

# stdlib modules
from sys import argv

# external modules
from pytest import main

# internal modules
from simu.thermo import (H0S0ReferenceState, HelmholtzState,
                         LinearHeatCapacity, StandardState, ThermoFactory)
from simu.utilities import assert_reproduction


def test_create_thermo_factory():
    """just create a ThermoFactory"""
    return ThermoFactory()


def test_register_contributions():
    """Create a ThermoFactory and register some contribtions"""
    create_frame_factory()


def test_create_frame():
    """create a ThermoFrame object"""
    fac = create_frame_factory()
    config = {
        "species": ["N2", "O2", "Ar", "CO2", "H2O"],
        "state": "HelmholtzState",
        "contributions": [
            "H0S0ReferenceState",
            "LinearHeatCapacity",
            "StandardState"
            ],
        }
    return fac.create_frame(config)


def test_parameter_structure():
    """Retrieve an (empty) parameter structure from created frame"""
    frame = create_simple_frame()
    assert_reproduction(dict(frame.parameters))


def test_parameter_names():
    """Retrieve parameter names from created frame"""
    frame = create_simple_frame()
    assert_reproduction(list(frame.parameters.view(flat=True).keys()))


def test_property_names():
    """Retrieve the names of defined properties from created frame"""
    frame = create_simple_frame()
    assert_reproduction(frame.property_names)


def test_call_frame():
    """Call a created frame with numerical values"""
    call_frame()


def test_relax():
    """Simple test of relaxation method"""
    frame, state, result = call_frame()
    delta_state = [-500, 10, -2, 1]
    assert frame.relax(result, delta_state) == state[0] / 500


# helper functions
def call_frame():
    """Call a frame object with a state and return all with the result"""
    frame = create_simple_frame()
    params = {'H0S0ReferenceState': {
                'dh_form': {'N2': 0.0, 'O2': 0.0},
                's_0': {'N2': 0.0, 'O2': 0.0},
                'T_ref': 298.15,
                'p_ref': 101325.0
                },
              'LinearHeatCapacity': {
                'cp_a': {'N2': 29.12379083, 'O2': 20.786},
                'cp_b': {'N2': 5.28694e-4, 'O2': 0.0}
                }
              }
    frame.parameters.set_struct_values(params)
    state = [398.15, 101325.0, 1, 1]
    result = frame(state)
    return frame, state, result


def create_frame_factory():
    """Create a ThermoFactory and register standard state contributions"""
    fac = test_create_thermo_factory()
    fac.register(H0S0ReferenceState, LinearHeatCapacity, StandardState)
    fac.register_state_definition(HelmholtzState)
    return fac


def create_simple_frame():
    """Create a ThermoFrame based on just standard state contributions"""
    fac = create_frame_factory()
    config = {
        "species": ["N2", "O2"],
        "state": "HelmholtzState",
        "contributions": [
            "H0S0ReferenceState",
            "LinearHeatCapacity",
            "StandardState"
            ],
        }
    return fac.create_frame(config)


if __name__ == "__main__":
    # only this file, very verbose and print stdout when started from here.
    main([__file__, "-v", "-v", "-s", "-rP"] + argv[1:])

