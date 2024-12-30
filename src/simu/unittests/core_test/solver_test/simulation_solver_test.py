from simu import NumericHandler, SimulationSolver
from simu.examples.material_model import Source
from simu.core.utilities import assert_reproduction




def test_instantiate():
    numeric = NumericHandler(Source.top())
    _ = SimulationSolver(numeric)


def test_solve():
    numeric = NumericHandler(Source.top())
    solver = SimulationSolver(numeric)
    solver.solve()