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
from simu.thermo import HelmholtzState, all_contributions


class ThermoNode(dict):
    """Class sharing a frame with other nodes of same type, and keeping
    the symbolic results by being a dictionary"""

    def __init__(self, frame: ThermoFrame):
        self.frame = frame
        state = SX.sym('x', 2 + len(self.frame.species))
        dict.__init__(self, zip(self.frame.property_names, self.frame(state)))


class MyThermoFactory(ThermoFactory):
    """This factory holds a ``ThermoFrame`` object for each defined model,
    ready to be deployed in a node."""

    def __init__(self):
        """Default and only constructor"""
        ThermoFactory.__init__(self)
        self.register(*all_contributions)
        self.register_state_definition(HelmholtzState)

        with open("frame_definitions.yml", encoding='UTF-8') as file:
            config = load(file, SafeLoader)
        with open("parameters.yml", encoding='UTF-8') as file:
            parameters = load(file, SafeLoader)

        # register frames for all defined configurations
        self.frames = {name: self.create_frame(cfg)
                       for name, cfg in config.items()}
        for frame in self.frames.values():
            frame.parameters = parameters

    def create_node(self, config_name: str) -> ThermoNode:
        """Based on the name of the thermodynamic model (ThermoFrame),
        create a wrapper object that hosts the properties for the instance."""
        return ThermoNode(self.frames[config_name])


def relax(nodes: dict, result: dict, delta_x: list) -> float:
    """find relaxation factor"""
    gamma = 0.9  # relaxation distance
    alpha = 1 / gamma
    idx = 0
    for name, node in nodes.items():
        length = node["state"].rows()
        new_alpha = node.frame.relax(result[name].values(),
                                    delta_x[idx:idx + length])
        alpha = min(alpha, new_alpha)
        idx += length
    return float(gamma * alpha)


def define_symbols(nodes: dict[str, ThermoNode], params: dict) -> dict:
    """Define the symbol dictionary to be calculated by casadi"""
    gas, liq = nodes["gas"], nodes["liq"]  # to define residuals easier
    return {"thermo-nodes": nodes,
            "residuals": vertcat(
                (gas["T"] - params["T"]) / 1e-7,
                (liq["T"] - params["T"]) / 1e-7,
                (gas["p"] - params["p"]) / 1,
                (liq["p"] - params["p"]) / 1,
                (gas["mu"] - liq["mu"]) / 1e-7,  # 2 equations
                (sum1(liq["n"]) - params["N"]) / 1e-7,
                (sum1(gas["n"]) - params["N"]) / 1e-7)
            }


def main():
    """Main entry of the script"""
    # create thermodynamic nodes
    factory = MyThermoFactory()
    nodes = {"liq": factory.create_node("Boston-Mathias-Redlich-Kwong-Liquid"),
             "gas": factory.create_node("Boston-Mathias-Redlich-Kwong-Gas")}
    params = {"T": 298.15 + 2, "p": 10e5, "x": 0.99, "y": 0.98, "N": 1}
    symbols = define_symbols(nodes, params)

    # add jacobian to list of required symbols
    state = vertcat(*[n["state"] for n in nodes.values()])
    symbols["dr_dx"] = jacobian(symbols["residuals"], state)

    # create a function with all the stuff
    symbols_flat = flatten_dictionary(symbols)
    func = Function("model", [state], symbols_flat.values())

    # now do the numerics
    state = array(nodes["liq"].frame.initial_state(params["T"], params["p"],
                  [params["x"], 1 - params["x"]]) +
                  nodes["gas"].frame.initial_state(params["T"], params["p"],
                  [params["y"], 1 - params["y"]]))

    for itr in range(30):
        # evaluate casadi function and unflatten dictionary
        res = unflatten_dictionary(dict(zip(symbols_flat, func(state))))

        # did it converge?
        residuals = ravel(res["residuals"])
        err = dot(residuals, residuals)
        if err < 1:
            print(f"{itr:2d}  end {log10(err):>5.1f}")
            break

        # do Newton step
        delta_x = -solve(res["dr_dx"], residuals)
        alpha = relax(nodes, res["thermo-nodes"], delta_x)
        print(f"{itr:2d} {alpha:4.2g} {log10(err):>5.1f}")

        # apply scaled step
        state += alpha * delta_x
    else:
        raise ValueError("Not converged, sorry!")


if __name__ == "__main__":
    main()
