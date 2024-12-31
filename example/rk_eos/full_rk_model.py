"""Create an instance of a full RK-EOS with standard state for an actual system
and use it for something interesting, e.g. draw a simple phase diagram
"""

# external modules
from numpy import dot, ravel, log10, hstack
from numpy.linalg import solve
from casadi import sum1, Function, jacobian, vertcat

# internal modules
from simu.core.utilities import flatten_dictionary, unflatten_dictionary

from example.rk_eos.rkt import ThermoNode, MyThermoFactory, relax


def define_symbols(nodes: dict[str, ThermoNode], params: dict) -> dict:
    """Define the symbol dictionary to be calculated by casadi"""
    gas, liq = nodes["gas"], nodes["liq"]  # to define residuals easier
    return {
        "thermo-nodes":
        nodes,
        "residuals":
        vertcat(
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
    nodes = {
        "liq": factory.create_node("LPG-RK-Liquid"),
        "gas": factory.create_node("LPG-RK-Gas")
    }
    params = {"T": 298.15 + 2, "p": 10e5, "N": 1}
    symbols = define_symbols(nodes, params)

    # add jacobian to list of required symbols
    state = vertcat(*[n["state"] for n in nodes.values()])
    symbols["dr_dx"] = jacobian(symbols["residuals"], state)

    # create a casadi function with all the stuff
    symbols_flat = flatten_dictionary(symbols)
    func = Function("model", [state], symbols_flat.values())

    # now do the numerics
    frames = {key: node.frame for key, node in nodes.items()}
    state = hstack([f.initial_state(*f.default) for f in frames.values()])

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
        alpha = relax(frames, res["thermo-nodes"], delta_x)
        print(f"{itr:2d} {alpha:4.2g} {log10(err):>5.1f}")

        # apply scaled step
        state += alpha * delta_x
    else:
        raise ValueError("Not converged, sorry!")


if __name__ == "__main__":
    main()
