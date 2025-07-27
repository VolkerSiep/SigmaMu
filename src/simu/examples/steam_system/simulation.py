from pprint import pprint

from simu import NumericHandler, SimulationSolver
from simu.examples.steam_system.process import SteamGeneration

def main():
    numeric = NumericHandler(SteamGeneration.top(), port_properties=False)

    solver = SimulationSolver(numeric)
    result = solver.solve()

    props = result.properties
    pprint(props["thermo_props"])
    pprint(props["model_props"])


if __name__ == '__main__':
    main()
