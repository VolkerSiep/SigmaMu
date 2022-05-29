#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Create an instance of a full RK-EOS with standard state for an actual system
and use it for something interesting, e.g. draw a simple phase diagram
"""

from yaml import load, SafeLoader

from simu.thermo import (
    ThermoFactory, HelmholtzState, H0S0ReferenceState, LinearHeatCapacity,
    StandardState, IdealMix, HelmholtzIdealGas)
from simu.thermo.cubic import (
    NonSymmetricMixingRule, LinearMixingRule,
    BostonMathiasAlphaFunction, CriticalParameters,
    RedlichKwongEOSLiquid, RedlichKwongEOSGas, RedlichKwongAFunction,
    RedlichKwongBFunction, RedlichKwongMFactor)


def create_factory():
    fac = ThermoFactory()
    fac.register(
        H0S0ReferenceState, LinearHeatCapacity, StandardState, IdealMix,
        HelmholtzIdealGas, NonSymmetricMixingRule, LinearMixingRule,
        BostonMathiasAlphaFunction, RedlichKwongEOSLiquid, RedlichKwongEOSGas,
        RedlichKwongAFunction, RedlichKwongBFunction, RedlichKwongMFactor,
        CriticalParameters)
    fac.register_state_definition(HelmholtzState)
    return fac

def create_frame(factory, name):
    with open("frame_definitions.yml") as file:
        config = load(file, SafeLoader)
    return factory.create_frame(config[name])

def main():
    fac = create_factory()
    gas = create_frame(fac, "BostonMathiasRedlichKwongGas")
    liq = create_frame(fac, "BostonMathiasRedlichKwongLiquid")

    with open("parameters.yml") as file:
        parameters = load(file, SafeLoader)

    gas.parameters = liq.parameters = parameters

    x_liq = liq.initial_state(298.15, 10e5, [0.5, 0.5])
    x_gas = gas.initial_state(298.15, 10e5, [0.5, 0.5])

    for name, prop in zip(liq.property_names, liq(x_liq)):
        print(name, prop)


# TODO:
#  - Update documentation
#  - to avoid need for custom derivatives for standard use, shall I
#    define a contribution that defines derivatives of properties
#    with respect to the state or parts of the state.
#
#    Well, I need mu_T = dmu_dT, mu_V = dmu_dV, S_T = dS_dT, p_V = dp_dV
#      and then the derivative of S, p, mu and above w.r.t. x.
#    How to handle that some derivatives are already there? -
#      Naming convention or just overwrite?

#     y_x
#     y\x
#     y@x
#     y|x
#     y!x
#     dy_dx  ddy_dx_dT  ddy_dx_dx  <- that's it!



if __name__ == "__main__":
    main()
