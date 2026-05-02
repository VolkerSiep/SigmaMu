from numpy.testing import assert_allclose
from simu.core.solver.linear import NumpySolver, ScaledLinearSparseSolver


def test_numpy_solver(linear_system):
    matrix, rhs, expected_x = linear_system
    solver = NumpySolver()
    x = solver.solve(matrix, rhs)
    assert_allclose(x, expected_x)


def test_scaled_linear_sparse_solver(linear_system):
    matrix, rhs, expected_x = linear_system
    solver = ScaledLinearSparseSolver()
    x = solver.solve(matrix, rhs)
    assert_allclose(x, expected_x)
