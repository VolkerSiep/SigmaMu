# -*- coding: utf-8 -*-
"""Test module for governing thermo objects"""

# stdlib modules
from sys import argv

# external modules
from pytest import main

# internal modules
from simu.thermo import (H0S0ReferenceState, HelmholtzState, ThermoFactory,
                         LinearHeatCapacity, StandardState, IdealMix,
                         HelmholtzIdealGas)
from simu.utilities import assert_reproduction, Quantity as Q


def test_create_thermo_factory():
    """just create a ThermoFactory"""
    ThermoFactory()


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
            "H0S0ReferenceState", "LinearHeatCapacity", "StandardState"
        ],
    }
    fac.create_frame(config)


def test_parameter_structure():
    """Retrieve an (empty) parameter structure from created frame"""
    frame = create_simple_frame()
    assert_reproduction(dict(frame.parameter_structure))


def test_property_structure():
    """Retrieve the names of defined properties from created frame"""
    frame = create_simple_frame()
    assert_reproduction(frame.property_structure)


def test_call_frame():
    """Call a created frame with numerical values"""
    result = call_frame()[2]
    result = {k: v for k, v in result.items() if k in {"S", "mu"}}
    assert_reproduction(result)


def test_relax():
    """Simple test of relaxation method"""
    frame, state, result = call_frame()
    delta_state = [-500, 10, -0.5, 1]
    beta = frame.relax(result, delta_state)
    assert beta == state[0] / 500


def test_initial_state():
    """Test whether initialisation of a Helmholtz ideal gas contributioon
    gives the correct volume"""
    frame = create_simple_frame()
    T, p, n = Q("25 degC"), Q("1 bar"), Q([1, 1], "mol")
    x = frame.initial_state(T, p, n, example_parameters)
    assert_reproduction(x[1])


# helper functions / data

example_parameters = {
    'H0S0ReferenceState': {
        'dh_form': {
            'N2': Q("0 J/mol"),
            'O2': Q("0 J/mol")
        },
        's_0': {
            'N2': Q("0 J/mol/K"),
            'O2': Q("0 J/mol/K")
        },
        'T_ref': Q("25 degC"),
        'p_ref': Q("1 atm")
    },
    'LinearHeatCapacity': {
        'cp_a': {
            'N2': Q("29.12379083 J/mol/K"),
            'O2': Q("20.786 J/mol/K")
        },
        'cp_b': {
            'N2': Q("5.28694e-4 J/mol/K**2"),
            'O2': Q("0.0 J/mol/K**2")
        }
    }
}


def call_frame():
    """Call a frame object with a state and return all with the result"""
    frame = create_simple_frame()
    state = Q([398.15, 0.0448, 1, 1])  # = T, V, *n
    result = frame(state, example_parameters)
    return frame, state, result


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
    return fac.create_frame(config)


if __name__ == "__main__":
    # only this file, very verbose and print stdout when started from here.
    main([__file__, "-v", "-v", "-s", "-rP"] + argv[1:])
