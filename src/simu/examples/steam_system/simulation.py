from pprint import pprint

from simu import NumericHandler, SimulationSolver
from process import SteamGeneration

def main():
    numeric = NumericHandler(SteamGeneration.top(), port_properties=False)
    # print(numeric.arguments)
    # pprint(numeric.export_state())

    solver = SimulationSolver(numeric)
    result = solver.solve()
    # pprint(result.properties["thermo_props"])
    result = solver.solve()


if __name__ == '__main__':
    main()
