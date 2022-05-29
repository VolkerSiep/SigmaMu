#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Create an instance of a full RK-EOS with standard state for an actual system
and use it for something interesting, e.g. draw a simple phase diagram
"""

from simu.thermo import (
    ThermoFactory, HelmholtzState, H0S0ReferenceState, LinearHeatCapacity,
    StandardState, IdealMix, HelmholtzIdealGas)
from simu.thermo.cubic import (
    NonSymmetricMixingRule, LinearMixingRule,
    BostonMathiasAlphaFunction, CriticalParameters)
from simu.thermo.cubic.rk import (
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

def create_frame(factory, phase):
    mixing_rule_a = {"name": "MixingRule_A",
                     "cls": "NonSymmetricMixingRule",
                     "options": {"target": "ceos_a"}}
    mixing_rule_b = {"name": "MixingRule_B",
                     "cls": "LinearMixingRule",
                     "options": {"target": "ceos_b"}}
    mixing_rule_c = {"name": "MixingRule_C",
                     "cls": "LinearMixingRule",
                     "options": {"target": "ceos_c",
                                 "source": "c_i",
                                 "src_mode": "parameter"}}

    config = {
        "species": ["C3", "nC4"],
        "state": "Helmholtz",
        "contributions": [
            'H0S0ReferenceState',
            'LinearHeatCapacity',
            'StandardState',
            'IdealMix',
            'HelmholtzIdealGas',
            'CriticalParameters',
            'RedlichKwongMFactor',
            'BostonMathiasAlphaFunction',
            'RedlichKwongAFunction',
            'RedlichKwongBFunction',
            mixing_rule_a, mixing_rule_b, mixing_rule_c,
            f'RedlichKwongEOS{phase.capitalize()}']
        }
    return factory.create_frame(config)

MISS = 0.0

parameters = {
    'H0S0ReferenceState': {
        'dh_form': {'C3': -104.7e3, 'nC4': -125.6e3},
        's_0': {'C3': MISS, 'nC4': MISS},
        'T_ref': 298.15,
        'p_ref': 1e5},
    'LinearHeatCapacity': {
        'cp_a': {'C3': MISS, 'nC4': MISS},
        'cp_b': {'C3': MISS, 'nC4': MISS}},
    'CriticalParameters': {
        'T_c': {'C3': 369.8, 'nC4': 425.2},
        'p_c': {'C3': 42.5e5, 'nC4': 38.0e5},
        'omega': {'C3': 0.153, 'nC4': 0.199}},
    'BostonMathiasAlphaFunction': {
        'eta': {'C3': 0, 'nC4': 0}},
    'MixingRule_A': {'T_ref': 298.15},
    'MixingRule_C': {
        'c_i': {'C3': 0.0, 'nC4': 0.0}},
    }


def main():
    fac = create_factory()
    gas, liq = create_frame(fac, "gas"), create_frame(fac, "liquid")
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
