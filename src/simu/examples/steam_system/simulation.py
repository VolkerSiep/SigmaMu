from pprint import pprint

from simu import NumericHandler, SimulationSolver
from process import SteamGeneration

def main():
    numeric = NumericHandler(SteamGeneration.top(), port_properties=False)

    solver = SimulationSolver(numeric)
    result = solver.solve(gamma=0.5)

    pprint(result.properties["thermo_props"])


if __name__ == '__main__':
    main()
