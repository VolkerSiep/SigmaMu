from numpy import linspace, log, min
from matplotlib import pyplot

parameters = {
    "T_0": 298.15,  # K
    "white": {
        "DH": 0,  # J/mol
        "S0": 51.55,  # 39,  # J/(mol.K)  # adapted to transition temperature
        "CP": 26.99  # J/(mol.K)
    },
    "grey": {
        "DH": -2090,  # J/mol
        "S0": 44.14,  # J/(mol.K)
        "CP": 25.77  # J/(mol.K)
    }
}


def mu(temperature, species):
    p = parameters[species]
    dh, s0, cp = p["DH"], p["S0"], p["CP"]
    t_0 = parameters["T_0"]
    return dh - temperature * s0 \
        + cp * (temperature - t_0 * (1 + log(temperature / t_0)))

def main():
    T = linspace(273.15, 298.15)
    pyplot.figure(figsize=(8, 4))
    pyplot.plot(T - 273.15, min([mu(T, "white"), mu(T, "grey")], axis=0) / 1000,
                "-", color="#cccccc", lw=10, label="stable phase")
    pyplot.plot(T - 273.15, mu(T, "white") / 1000, "k-", label="white tin")
    pyplot.plot(T - 273.15, mu(T, "grey") / 1000, "b-", label="grey tin")
    pyplot.plot([13.2, 13.2], [-15.4, -14], "g:", label="Transition temperature")
    pyplot.plot([8.83, 8.83], [-15.4, -14], "r:", label="Calculated transition")
    pyplot.grid()
    pyplot.xlim([0, 25])
    pyplot.ylim([-15.4, -14])
    pyplot.legend()
    pyplot.xlabel("Temperature [$^\\circ$C]")
    pyplot.ylabel("Chemical potential [kJ/(mol$\\,$K)]")
    pyplot.show()


if __name__ == '__main__':
    main()