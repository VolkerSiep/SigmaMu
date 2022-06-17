"""Create again the pair of liquid and gas RK-EOS. Then query parameter
sensitivity."""

from numpy import dot, hstack, linspace, ravel, log10
from numpy.linalg import solve
from casadi import Function, jacobian, SX, vertcat

from simu.utilities import (
    FlexiDictionary, flatten_dictionary, unflatten_dictionary)
from examples.rk_eos.rkt import MyThermoFactory, relax


def p_sat(temperature):
    """Return water saturation pressure [Pa] for given temperature [K]"""
    p_a, p_b, p_c = 3.55959, 643.748, -198.043
    return 1e5 * 10 ** (p_a - (p_b / (temperature + p_c)))


def define_symbols(nodes, thermo_parameters):
    """Define the symbol dictionary to be calculated by casadi"""
    gas, liq = nodes["gas"], nodes["liq"]  # to define residuals easier
    temperature_spec = SX.sym("T")
    pressure = p_sat(temperature_spec)
    state = vertcat(*[n["state"] for n in nodes.values()])

    residuals = FlexiDictionary({
        "gas temp": (gas["T"] - temperature_spec) / 1e-7,
        "liq temp": (liq["T"] - temperature_spec) / 1e-7,
        "gas pres": (gas["p"] - pressure) / 1,
        "liq pres": (liq["p"] - pressure) / 1,
        "liq size": (liq["n"] - 1) / 1e-7,
        "gas size": (gas["n"] - 1) / 1e-7},
        values_are_symbols=True, symbol=True)

    properties = FlexiDictionary({
        "dmu": gas["mu"][0] - liq["mu"][0]},
        values_are_symbols=True)

    parameters = FlexiDictionary({
        "T_spec": temperature_spec},
        values_are_symbols=True, symbol=True)

    symbols = {
        "thermo-nodes": nodes,
        "residuals": residuals,
        "property": properties,
        "parameters": parameters,
        "jacobians": {
            "dr_dx": jacobian(residuals.symbols, state),
            "dr_dp": jacobian(residuals.symbols, thermo_parameters),
            "dy_dx": jacobian(properties.symbols, state),
            "dy_dp": jacobian(properties.symbols, thermo_parameters)
        }
    }
    return state, symbols


def main():
    """Main entry point of the script"""

    factory = MyThermoFactory()
    nodes = {"liq": factory.create_node("Water-RK-Liquid"),
             "gas": factory.create_node("Water-RK-Gas")}

    # use common thermodynamic parameters in both models
    nodes["liq"].frame.parameters.symbols = \
        nodes["gas"].frame.parameters.symbols

    # find parameter to optimise on
    syms = nodes["liq"].frame.parameters.view(flat=True, symbol=True)
    th_param = syms['BostonMathiasAlphaFunction.eta.H2O']

    state, symbols = define_symbols(nodes, th_param)
    print(symbols["residuals"])

    # create casadi-function from symbols
    param = symbols["parameters"].symbols  # process parameters (T_spec)
    symbols_flat = flatten_dictionary(symbols)
    func = Function("model", [state, param, th_param], symbols_flat.values())

    temperatures = linspace(373, 640, num=11)

    # now do the numerics
    state = hstack([node.frame.initial_state(*node.frame.default)
                   for node in nodes.values()])
    eta = 0
    temperature = temperatures[0]  # TODO: for temperature in temperatures:


    for itr in range(10):
        # evaluate casadi function and unflatten dictionary
        res = func(state, temperature, eta)
        res = unflatten_dictionary(dict(zip(symbols_flat, res)))

        # did it converge?
        residuals = ravel(vertcat(*res["residuals"].values()))
        # residuals = ravel(res["residuals"].values())
        err = dot(residuals, residuals)
        if err < 1:
            print(f"{itr:2d}  end {log10(err):>5.1f}")
            break

        # determine Newton step
        delta_x = -solve(res["jacobians"]["dr_dx"], residuals)
        alpha = relax(nodes, res["thermo-nodes"], delta_x)
        print(f"{itr:2d} {alpha:4.2g} {log10(err):>5.1f}")

        # apply scaled step
        state += alpha * delta_x
    else:
        raise ValueError("Not converged, sorry!")

    # TODO: now extract total derivative (derive equations first)



if __name__ == "__main__":
    main()
