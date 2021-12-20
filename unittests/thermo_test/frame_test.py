# -*- coding: utf-8 -*-

# stdlib modules
from sys import path
from pathlib import Path

# reproductiontest
path.append(str(Path(__file__).absolute().parents[1]))
from reproductiontest import assert_reproduction


def test_create_thermo_factory():
    from mushell.thermo import ThermoFactory
    return ThermoFactory()

def test_register_contributions():
    create_frame_factory()

def test_create_frame():
    fac = create_frame_factory()
    config = {
        "species": ["N2", "O2", "Ar", "CO2", "H2O"],
        "contributions": [
            "state#Helmholtz",
            "reference_state#H0S0",
            "heat_capacity#linear",
            "standard_state"
            ],
        }
    return fac.create_frame(config)

def test_parameter_structure():
    frame = create_simple_frame()
    assert_reproduction(frame.parameters)

def test_parameter_names():
    frame = create_simple_frame()
    assert_reproduction(frame.parameter_names)

def test_property_names():
    frame = create_simple_frame()
    assert_reproduction(frame.property_names)

def test_call_frame():
    call_frame()

def test_relax():
    frame, state, result = call_frame()
    delta_state = [-500, 10, -2, 1]
    assert frame.relax(result, delta_state) == state[0] / 500


# helper functions
def call_frame():
    """Call a frame object with a state and return all with the result"""
    frame = create_simple_frame()
    params = {'H0S0 reference state': {
                'dh_form': {'N2': 0.0, 'O2': 0.0},
                's_0': {'N2': 0.0, 'O2': 0.0},
                'T_ref': 298.15,
                'p_ref': 101325.0
                },
              'Linear heat capacity': {
                'cp_a': {'N2': 29.12379083, 'O2': 20.786},
                'cp_b': {'N2': 5.28694e-4, 'O2': 0.0}
                }
             }
    frame.parameters = params
    state = [398.15, 101325.0, 1, 1]
    result = frame(state)
    return frame, state, result

def create_frame_factory():
    """Create a ThermoFactory and register standard state contributions"""
    from mushell.thermo import (HelmholtzState, H0S0ReferenceState,
                                LinearHeatCapacity, StandardState)
    fac = test_create_thermo_factory()
    cont = [HelmholtzState, H0S0ReferenceState, LinearHeatCapacity,
            StandardState]
    for con in cont:
        fac.register_contribution(con)
    print(fac.contribution_names)
    return fac

def create_simple_frame():
    """Create a ThermoFrame based on just standard state contributions"""
    fac = create_frame_factory()
    config = {
        "species": ["N2", "O2"],
        "contributions": [
            "state#Helmholtz",
            "reference_state#H0S0",
            "heat_capacity#linear",
            "standard_state"
            ],
        }
    return fac.create_frame(config)

if __name__ == "__main__":
    from pytest import main
    # only this file, very verbose and print stdout when started from here.
    main([__file__, "-v", "-v", "-rP"])
