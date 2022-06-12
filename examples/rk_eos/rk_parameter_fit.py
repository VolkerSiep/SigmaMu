"""Create again the pair of liquid and gas RK-EOS. Then query parameter
sensitivity."""

from numpy import dot, hstack, linspace, ravel, log10
from casadi import Function, jacobian, SX, vertcat

from simu.utilities import flatten_dictionary, unflatten_dictionary
from rkt import MyThermoFactory


def p_sat(temperature):
    """Return water saturation pressure [Pa] for given temperature [K]"""
    p_a, p_b, p_c = 3.55959, 643.748, -198.043
    return 1e5 * 10 ** (p_a - (p_b / (temperature + p_c)))


def define_symbols(nodes, parameters):
    """Define the symbol dictionary to be calculated by casadi"""
    gas, liq = nodes["gas"], nodes["liq"]  # to define residuals easier
    temperature_spec = SX.sym("T")
    pressure = p_sat(temperature_spec)
    state = vertcat(*[n["state"] for n in nodes.values()])

    residuals = vertcat(
                (gas["T"] - temperature_spec) / 1e-7,
                (liq["T"] - temperature_spec) / 1e-7,
                (gas["p"] - pressure) / 1,
                (liq["p"] - pressure) / 1,
                (liq["n"] - 1) / 1e-7,
                (gas["n"] - 1) / 1e-7)

    prop =  gas["mu"][0] - liq["mu"][0]

    return {
        "thermo-nodes": nodes,
        "residuals": residuals,
        "property": prop,
        "parameters": {
            "temperature": temperature_spec
        },
        "derivatives": {
            "dr_dx": jacobian(residuals, state),
            "dr_dp": jacobian(residuals, parameters),
            "dy_dx": jacobian(prop, state),
            "dy_dp": jacobian(prop, parameters)
        }
    }

# TODO: shall I make a general dictionary-like object that can flatten
#   and unflatten itself? The ThermoParameter class is maybe that?
#   but that dictionary cannot use non-scalar numerical objects
#   (so long at least)


def main():
    """Main entry point of the script"""

    factory = MyThermoFactory()
    nodes = {"liq": factory.create_node("Water-RK-Liquid"),
             "gas": factory.create_node("Water-RK-Gas")}

    # use common thermodynamic parameters in both models
    nodes["liq"].frame.parameters.symbols = nodes["gas"].frame.parameters.symbols

    # find parameter to optimise on
    syms = nodes["liq"].frame.parameters.flat_symbols
    parameter = syms['BostonMathiasAlphaFunction.eta.H2O']

    symbols = define_symbols(nodes, parameter)

    # create casadi-function from symbols
    state = vertcat(*[n["state"] for n in nodes.values()])
    param = vertcat(*flatten_dictionary(symbols["parameters"]).values())
    symbols_flat = flatten_dictionary(symbols)
    func = Function("model", [state, param, parameter], symbols_flat.values())

    temperatures = linspace(373, 640, num=11)

    # now do the numerics
    state = hstack([node.frame.initial_state(*node.frame.default)
                   for node in nodes.values()])
    eta = 0
    temperature = temperatures[0]


    for itr in range(1):
        # evaluate casadi function and unflatten dictionary
        res = func(state, temperature, eta)
        res = unflatten_dictionary(dict(zip(symbols_flat, res)))

        # did it converge?
        residuals = ravel(res["residuals"])
        err = dot(residuals, residuals)
        if err < 1:
            print(f"{itr:2d}  end {log10(err):>5.1f}")
            break




if __name__ == "__main__":
    main()
