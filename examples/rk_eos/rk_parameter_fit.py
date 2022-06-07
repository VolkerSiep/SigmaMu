"""Create again the pair of liquid and gas RK-EOS. Then query parameter
sensitivity."""

from numpy import hstack, linspace, ravel
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

    return {"thermo-nodes": nodes,
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


def main():
    """Main entry point of the script"""

    factory = MyThermoFactory()
    nodes = {"liq": factory.create_node("Water-RK-Liquid"),
             "gas": factory.create_node("Water-RK-Gas")}

    # TODO: can I make ThermoFrame to optionally do this stuff for me?
    #   (there is some cleaning-up to do!!! User should not have to
    #    switch all the time between flattened and unflattened dictionary)

    # Parameters deserve to be their own object that can hold names, symbols
    #  and actual values at the same time, and also support to


    parameter_names = nodes["gas"].frame.parameter_names
    parameter_sym = SX.sym("parameters", len(parameter_names))
    parameters = unflatten_dictionary(dict(zip(parameter_names, parameter_sym)))
    for node in nodes.values():
        node.parameters = parameters

    symbols = define_symbols(nodes, parameter_sym)

    # create casadi-function from symbols
    state = vertcat(*[n["state"] for n in nodes.values()])
    param = vertcat(*flatten_dictionary(symbols["parameters"]).values())
    symbols_flat = flatten_dictionary(symbols)
    func = Function("model", [state, param, parameters], symbols_flat.values())

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
        print(itr, residuals)


if __name__ == "__main__":
    main()
