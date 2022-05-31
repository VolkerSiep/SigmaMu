#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Create an instance of a full RK-EOS with standard state for an actual system
and use it for something interesting, e.g. draw a simple phase diagram
"""

from yaml import load, SafeLoader

from numpy import array, ravel, log10
from numpy.linalg import solve
from casadi import SX, sum1, Function, jacobian, vertcat

from simu.thermo import ThermoFactory
from simu.utilities import flatten_dictionary, unflatten_dictionary


class ThermoNode(dict):
    def __init__(self, frame):
        self.frame = frame
        state = SX.sym('x', 2 + len(self.frame.species))
        items = zip(self.frame.property_names, self.frame(state))
        self.update({name: value for name, value in items})

    def initialise(self, T: float, p: float, n: list[float]) -> list[float]:
        return self.frame.initial_state(T, p, n)

    def relax_state(self, props, delta):
        # TODO: stupid that I need to flatten the dict again!
        props = flatten_dictionary(props)
        return self.frame.relax(props.values(), delta)


class MyThermoFactory(ThermoFactory):
    def __init__(self):
        ThermoFactory.__init__(self)
        from simu.thermo import HelmholtzState, all_contributions
        self.register(*all_contributions)
        self.register_state_definition(HelmholtzState)

        with open("frame_definitions.yml") as file:
            config = load(file, SafeLoader)
        with open("parameters.yml") as file:
            parameters = load(file, SafeLoader)

        # register frames for all defined configurations
        self.frames = {name: self.create_frame(cfg)
                       for name, cfg in config.items()}
        for frame in self.frames.values():
            frame.parameters = parameters

    def create_node(self, config_name: str) -> ThermoNode:
        return ThermoNode(self.frames[config_name])


def main():
    # create thermodynamic nodes
    factory = MyThermoFactory()
    gas = factory.create_node("Boston-Mathias-Redlich-Kwong-Gas")
    liq = factory.create_node("Boston-Mathias-Redlich-Kwong-Liquid")

    T, p, x, y = 298.15+2.0, 10e5, 0.99, 0.98

    symbols = {"thermo-nodes": {
                   "gas": gas,
                   "liq": liq
                   },
               "residuals": {
                    "T-gas": gas["T"] - T,
                    "T-liq": liq["T"] - T,
                    "p-gas": gas["p"] - p,
                    "p-liq": liq["p"] - p,
                    "equil": gas["mu"] - liq["mu"],  # 2 equations
                    "N-liq": sum1(liq["n"]) - 1,
                    "N-gas": sum1(gas["n"]) - 1,
                    }
               }
    # add jacobian to list of required symbols
    state = vertcat(liq["state"], gas["state"])
    residuals = vertcat(*symbols["residuals"].values())
    symbols["residuals"] = residuals  # overwrite residuals
    symbols["dr_dx"] = jacobian(residuals, state)

    # create a function with all the stuff
    symbols_flat = flatten_dictionary(symbols)
    func = Function("model", [state], symbols_flat.values())

    # now do the numerics
    state = array(liq.initialise(T, p, [x, 1 - x]) +
                  gas.initialise(T, p, [y, 1 - y]))
    gamma = 0.9

    for iter in range(30):
        items = zip(symbols_flat, func(state))
        res = unflatten_dictionary({name: value for name, value in items})

        residuals = ravel(res["residuals"])
        delta_x = -solve(res["dr_dx"], residuals)
        alphas = [gas.relax_state(res["thermo-nodes"]["gas"], delta_x[4:]),
                  liq.relax_state(res["thermo-nodes"]["liq"], delta_x[:4]),
                  1 / gamma]
        alpha = gamma * min(alphas)
        state += alpha * delta_x
        print(f"{iter}\t{alpha:.2g}\t{log10(max(residuals)+1e-10):.1f}")


# TODO:
#    - also validate correctness of this binary phase diagram at 10 bar
#      maybe trace the range between the pure boiling points
#       propane: T_boil = 300.09 K ( propane = x)
#       n-butane: T_boil = 352.62 (butane = 1-x)

if __name__ == "__main__":
    main()
