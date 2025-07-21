# stdlib
from __future__ import annotations
from typing import TYPE_CHECKING
from collections.abc import Sequence

# external
from numpy import array, squeeze
from numpy.linalg import solve
from casadi import SX, vertcat, jacobian, Function, jtimes

# internal
from simu.core.utilities.structures import flatten_dictionary

if TYPE_CHECKING:
    from simu.core.thermo.frame import ThermoFrame
    from simu.core.thermo.state import InitialState
    from simu.core.utilities.quantity import Quantity
    from simu.core.utilities.types import NestedMap

GAMMA = 0.9
MAX_ITER = 30
REL_TOL = 1e-9


def __define_functions(frame: ThermoFrame, parameters: NestedMap[Quantity],
                      initial_state: InitialState) -> (Function, Function):
    x = SX.sym("x", len(frame.species) + 2)
    result = frame(x, parameters, squeeze_results=False)
    props = result["props"]

    t_res = (props["T"] / initial_state.temperature).to("").m - 1
    p_res = (props["p"] / initial_state.pressure).to("").m - 1
    n_res = (props["n"] / initial_state.mol_vector).to("").m - 1

    def unpack(name):
        flat = flatten_dictionary(result.get("name", {}))
        return [v.m for v in flat.values()]

    bounds = vertcat(*unpack("bounds"))
    residuals = vertcat(*[t_res, p_res, n_res])

    dr_dx = jacobian(residuals, x)

    # function for relaxation factor
    dx = SX.sym("dx", len(frame.species) + 2)
    alpha = -bounds / jtimes(bounds, x, dx)

    return (Function("r", [x], [residuals, dr_dx], ["x"], ["r", "dr_dx"]),
            Function("a", [x, dx], [alpha], ["x", "dx"], ["alpha"]))


def refine_initial_state(
        frame: ThermoFrame,
        parameters: NestedMap[Quantity],
        initial_state: InitialState,
        estimate: Sequence[float]) -> Sequence[float]:
    r"""
    The strategy for initialising thermodynamic models relies on the demand that
    all models can provide an estimate of their initial state based on
    temperature, pressure, and molar quantities, such that a Newton-Raphson
    solver, respecting the model's domain boundaries, will find a refinement to
    exactly meet the specified :math:`T, p, \vec n` specifications.

    This function implements the specialised version of the solver, refining the
    initial estimate into a converged result.

    :param frame: The thermodynamic model
    :param parameters: The current set of thermodynamic parameters
    :param initial_state: The desired initial state in terms of
      :math:`T, p, \vec n`
    :param estimate: The estimate of the internal state vector.
    :return: The refined internal state vector

    This function is not to be called for Gibbs models, as the internal state
    is trivially defined in :math:`T, p, \vec n` . The most common application
    is to refine Helmholtz models on :math:`T, V, \vec n`. Though the
    non-trivial part is in that case reduced to finding
    :math:`p(V) - p_\mathrm{init} = 0`, the implementation is kept general in
    order to allow for arbitrary internal state definitions.
    """
    f_r, f_a = __define_functions(frame, parameters, initial_state)
    for iteration in range(30):
        r, dr_dx = f_r(estimate)
        if r.T @ r < REL_TOL ** 2:
            break
        dx = -squeeze(solve(dr_dx, r))
        a = squeeze(array(f_a(estimate, dx)))
        a = a[0 < a]
        if len(a):
            alpha = min(min(a) * GAMMA, 1.0)
            dx = alpha * dx
        estimate = estimate + dx
    else:
        msg = "State estimate not sufficiently close for convergence"
        raise ValueError(msg)
    return estimate
