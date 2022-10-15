"""Create again the pair of liquid and gas RK-EOS. Then query parameter
sensitivity."""

from numpy import array, dot, hstack, linspace, ravel, exp
from numpy.linalg import solve, LinAlgError
from casadi import Function, jacobian, SX, vertcat
import pylab

from simu.utilities import (FlexiDict, flatten_dictionary,
                            unflatten_dictionary)
from examples.rk_eos.rkt import MyThermoFactory, relax


def p_sat(temperature):
    """Return water saturation pressure [Pa] for given temperature [K]"""
    # p_a, p_b, p_c = 3.55959, 643.748, -198.043
    # return 1e5 * 10 ** (p_a - (p_b / (temperature + p_c)))

    # p_a, p_b, p_c = 0.61078, 17.27, 237.3
    # t_cel = temperature - 273.15
    # pressure_kpa = p_a * exp(p_b * t_cel / (t_cel + p_c))
    # return pressure_kpa * 1000
    t_c, p_c = 647.096, 22.06e6
    tau = 1 - temperature / t_c
    term = (-7.85951783 * tau + 1.84408259 * tau**1.5 - 11.7866497 * tau**3 +
            22.6807411 * tau**3.5 - 15.9618719 * tau**4 +
            1.80122502 * tau**7.5)
    return p_c * exp(t_c / temperature * term)


def define_symbols(frames):
    """Define the symbol dictionary to be calculated by casadi and by that
    define the entire model. This function returns all symbols defined,
    including the independent ones (states, parameters and thermodynamic
    parameters)"""

    thermo_parameters = frames["gas"].parameters
    states = {n: f.create_sym_state() for n, f in frames.items()}
    state = vertcat(*states.values())

    def func(frame, state):
        raw = frame.func(state, thermo_parameters.symbols)
        return dict(zip(frame.property_names, raw))

    nodes = {n: func(f, states[n]) for n, f in frames.items()}

    gas, liq = nodes["gas"], nodes["liq"]  # to define residuals easier
    temperature_spec = SX.sym("T")
    pressure = p_sat(temperature_spec)

    residuals = FlexiDict(
        {
            "gas temp": (gas["T"] - temperature_spec) / 1e-7,
            "liq temp": (liq["T"] - temperature_spec) / 1e-7,
            "gas pres": (gas["p"] - pressure) / 1,
            "liq pres": (liq["p"] - pressure) / 1,
            "liq size": (liq["n"] - 1) / 1e-7,
            "gas size": (gas["n"] - 1) / 1e-7
        },
        symbol=True)

    properties = FlexiDict({"dmu": gas["mu"][0] - liq["mu"][0]}, symbol=True)

    parameters = FlexiDict({"T_spec": temperature_spec}, symbol=True)

    symbols = {
        "independent": {
            "state": state,
            "parameters": parameters,
            "thermo_parameters": thermo_parameters.view(flat=False,
                                                        symbol=True),
        },
        "thermo-nodes": nodes,
        "residuals": residuals,
        "properties": properties,
        "jacobians": {
            "dr_dx": jacobian(residuals.symbols, state),
            "dr_dp": jacobian(residuals.symbols, thermo_parameters.symbols),
            "dy_dx": jacobian(properties.symbols, state),
            "dy_dp": jacobian(properties.symbols, thermo_parameters.symbols)
        }
    }
    return symbols


def solve_one_point(frames, state, function, parameters, th_param):
    """Solve one data point and calculate delta mu.

    Actually, this is already independent of the problem posed."""
    for _ in range(30):
        # evaluate casadi function and unflatten dictionary
        res = function(state, parameters, th_param)

        # did it converge?
        residuals = ravel(vertcat(*res["residuals"].values()))
        err = dot(residuals, residuals)
        if err < 1:
            # print(f"{itr:2d}  end {log10(err):>5.1f}")
            break

            # determine Newton steplog10
        delta_x = -solve(res["jacobians"]["dr_dx"], residuals)
        alpha = relax(frames, res["thermo-nodes"], delta_x)
        # print(f"{itr:2d} {alpha:4.2g} {log10(err):>5.1f}")

        # apply scaled step
        state += alpha * delta_x
    else:
        raise ValueError("Not converged, sorry!")
    return res


def main():
    """Main entry point of the script"""

    factory = MyThermoFactory()
    frames = {
        "gas": factory.frames["Water-RK-Gas"],
        "liq": factory.frames["Water-RK-Liquid"]
    }
    # th_param = syms['BostonMathiasAlphaFunction.eta.H2O']

    symbols = define_symbols(frames)

    # create casadi-function from symbols
    indep = symbols["independent"]
    state = indep["state"]
    param = indep["parameters"].symbols  # proc. param.: T_spec
    th_param = indep["thermo_parameters"].symbols
    symbols_flat = flatten_dictionary(symbols)
    func = Function("model", [state, param, th_param], symbols_flat.values())

    def fun(state, param, th_param):
        raw = func(state, param, th_param)
        return unflatten_dictionary(dict(zip(symbols_flat, raw)))

    temperatures = linspace(273, 600, num=21)

    # now do the numerics
    th_param = indep["thermo_parameters"].flat_values

    errs = []
    temperatures_p = []
    gas_const = 8.314426

    for _ in range(3):
        state = hstack(
            [frame.initial_state(*frame.default) for frame in frames.values()])
        obj = 0  # objective function (just to monitor)
        grad_obj = 0  # gradient of objective function
        hess_obj = 0  # hessian of objective function
        errs.append([])
        temperatures_p.append([])
        for temperature in temperatures:
            state = hstack([
                frame.initial_state(temperature, p_sat(temperature), [1])
                for frame in frames.values()
            ])
            try:
                res = solve_one_point(frames, state, fun, temperature,
                                      th_param)
            except (ValueError, LinAlgError):
                print(f"T = {temperature} failed to converge")
                continue
            dmu = res['properties']['dmu']
            err = exp(float(dmu) / (gas_const * temperature)) - 1
            errs[-1].append(err * 100)
            temperatures_p[-1].append(temperature)
            jac_yx, jac_yp, jac_rp, jac_rx = \
                [res["jacobians"][i] for i in "dy_dx dy_dp dr_dp dr_dx".split()]
            try:
                # calculate total differential for parameter sensitivity
                dy_dp = jac_yp - dot(jac_yx, solve(jac_rx, jac_rp))
            except LinAlgError:
                # print(f"T = {temperature} failed to show observability")
                continue
            obj += 0.5 * dmu * dmu
            grad_obj += dot(dy_dp.T, dmu)
            hess_obj += dot(dy_dp.T, dy_dp)

        d_eta = -grad_obj[0, 0] / hess_obj[0, 0]
        print("eta = ", th_param[0], f"Objective: {obj}", "delta_eta = ",
              d_eta)
        th_param[0] += d_eta

    for k, err in enumerate(errs):
        tsk = array(temperatures_p[k])
        # pylab.plot(tsk, p_sat(tsk) * dmu, "-", label=str(k))
        pylab.plot(tsk, err, "-", label=str(k))
    pylab.grid()
    pylab.legend(loc="best")
    pylab.show()


if __name__ == "__main__":
    main()
