#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Create an instance of a full RK-EOS with standard state for an actual system
and use it for something interesting, e.g. draw a simple phase diagram
"""

from yaml import load, SafeLoader

from numpy import array, dot, ravel, log10
from numpy.linalg import solve
from casadi import SX, sum1, Function, jacobian, vertcat

from simu.thermo import ThermoFactory, ThermoFrame
from simu.utilities import flatten_dictionary, unflatten_dictionary


class ThermoNode(dict):
    """Class sharing a frame with other nodes of same type, and keeping
    the symbolic results by being a dictionary"""

    def __init__(self, frame: ThermoFrame):
        self.frame = frame
        state = SX.sym('x', 2 + len(self.frame.species))
        items = zip(self.frame.property_names, self.frame(state))
        self.update({name: value for name, value in items})


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
    nodes = {"liq": factory.create_node("Boston-Mathias-Redlich-Kwong-Liquid"),
             "gas": factory.create_node("Boston-Mathias-Redlich-Kwong-Gas")}

    T, p, x, y = 298.15+2.0, 10e5, 0.99, 0.98

    gas, liq = nodes["gas"], nodes["liq"]  # to define residuals easier
    symbols = {"thermo-nodes": nodes,
               "residuals": vertcat(
                    (gas["T"] - T) / 1e-7,
                    (liq["T"] - T) / 1e-7,
                    (gas["p"] - p) / 1,
                    (liq["p"] - p) / 1,
                    (gas["mu"] - liq["mu"]) / 1e-7,  # 2 equations
                    (sum1(liq["n"]) - 1) / 1e-7,
                    (sum1(gas["n"]) - 1) / 1e-7)
               }
    # add jacobian to list of required symbols
    state = vertcat(*[n["state"] for n in nodes.values()])
    symbols["dr_dx"] = jacobian(symbols["residuals"], state)

    # create a function with all the stuff
    symbols_flat = flatten_dictionary(symbols)
    func = Function("model", [state], symbols_flat.values())

    # now do the numerics
    state = array(liq.frame.initial_state(T, p, [x, 1 - x]) +
                  gas.frame.initial_state(T, p, [y, 1 - y]))

    gamma = 0.9  # relaxation distance

    for iter in range(30):
        # evaluate casadi function and unflatten dictionary
        items = zip(symbols_flat, func(state))
        res = unflatten_dictionary({name: value for name, value in items})

        # did it converge?
        residuals = ravel(res["residuals"])
        err = dot(residuals, residuals)
        if err < 1:
            print(f"{iter:2d}  end {log10(err):>5.1f}")
            break

        # do Newton step
        delta_x = -solve(res["dr_dx"], residuals)

        # find relaxation factor
        alpha = 1 / gamma
        idx = 0
        for name, node in nodes.items():
            length = node["state"].rows()
            new_alpha = node.frame.relax(res["thermo-nodes"][name].values(),
                                         delta_x[idx:idx + length])
            alpha = min(alpha, new_alpha)
            idx += length
        alpha = float(gamma * alpha)

        # apply scaled step
        state += alpha * delta_x
        print(f"{iter:2d} {alpha:4.2g} {log10(err):>5.1f}")
    else:
        raise ValueError("Not converged, sorry!")

# TODO:
#    - validate correctness of this binary phase diagram at 10 bar
#      maybe trace the range between the pure boiling points
#       propane: T_boil = 300.09 K ( propane = x)
#       n-butane: T_boil = 352.62 (butane = 1-x)


if __name__ == "__main__":
    main()
