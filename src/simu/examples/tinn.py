from numpy import linspace, log, min
from matplotlib import pyplot

parameters = {
    "T_0": 298.15,  # K
    "white": {
        "DH": 0,  # J/mol
        "S0": 51.39,  # J/(mol.K)  # adapted to transition temperature
        "CP": 26.99  # J/(mol.K)
    },
    "grey": {
        "DH": -2090,  # J/mol
        "S0": 44.14,  # J/(mol.K)
        "CP": 25.77  # J/(mol.K)
    }
}


def mu(temperature, species):
    dh = parameters[species]["DH"]
    s0 = parameters[species]["S0"]
    cp = parameters[species]["CP"]
    t_0 = parameters["T_0"]
    return dh - temperature * s0 \
        + cp * (temperature - t_0 - log(temperature / t_0))

def main():
    T = linspace(278.15, 293.15)
    pyplot.figure(figsize=(8, 4))
    pyplot.plot(T - 273.15, min([mu(T, "white"), mu(T, "grey")], axis=0) / 1000,
                "-", color="#dddddd", lw=10, label="stable phase")
    pyplot.plot(T - 273.15, mu(T, "white") / 1000, "k-", label="white tinn")
    pyplot.plot(T - 273.15, mu(T, "grey") / 1000, "b-", label="grey tinn")
    pyplot.plot([13.2, 13.2], [-15.2, -14.8], "r:", label="Transition temperature")
    pyplot.grid()
    pyplot.xlim([5, 20])
    pyplot.ylim([-15.2, -14.8])
    pyplot.legend()
    pyplot.xlabel("Temperature [$^\\circ$C]")
    pyplot.ylabel("Chemical potential [kJ/(mol$\\,$K)]")
    pyplot.show()


if __name__ == '__main__':
    main()