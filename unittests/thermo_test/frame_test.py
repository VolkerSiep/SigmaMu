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
    from mushell.thermo import (HelmholtzState, H0S0ReferenceState,
                                LinearHeatCapacity, StandardState)
    fac = test_create_thermo_factory()
    cont = {"Helmholtz state": HelmholtzState,
            "H0S0 reference state": H0S0ReferenceState,
            "Linear heat capacity": LinearHeatCapacity,
            "Standard state": StandardState}
    for name, con in cont.items():
        fac.register_contribution(name, con)
    return fac

def test_create_frame():
    fac = test_register_contributions()
    config = {
        "species": ["N2", "O2", "Ar", "CO2", "H2O"],
        "contributions": [
            "Helmholtz state",
            "H0S0 reference state",
            "Linear heat capacity",
            "Standard state",
            ],
        }
    return fac.create_frame(config)

def test_parameter_structure():
    frame = test_create_frame()
    assert_reproduction(frame.parameters)

def test_parameter_names():
    frame = test_create_frame()
    assert_reproduction(frame.parameter_names)

def test_property_names():
    frame = test_create_frame()
    assert_reproduction(frame.property_names)


# TODO: test to actually call the function

if __name__ == "__main__":
    from pytest import main
    from sys import argv
    # only this file, very verbose and print stdout when started from here.
    argv.extend([__file__, "-v", "-v", "-rP"])
    main()
