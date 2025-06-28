from simu import NumericHandler, SimulationSolver

from process import SteamGeneration

def main():
    numeric = NumericHandler(SteamGeneration.top())
    # print(numeric.arguments)

    solver = SimulationSolver(numeric)
    result = solver.solve()

if __name__ == '__main__':
    main()