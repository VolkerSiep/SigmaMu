from pprint import pprint
from pathlib import Path
from sys import path

from simu import NumericHandler, SimulationSolver
from process import SteamGeneration

def main():
    path.insert(0, str(Path(__file__).parent))

    numeric = NumericHandler(SteamGeneration.top(), port_properties=False)

    solver = SimulationSolver(numeric)
    result = solver.solve()

    props = result.properties
    pprint(props["thermo_props"])
    pprint(props["model_props"])


if __name__ == '__main__':
    main()
