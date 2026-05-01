from __future__ import annotations
from typing import TYPE_CHECKING, Any
from pydantic import BaseModel, Field
from numpy import sqrt, zeros_like
from numpy.linalg import norm, solve
from scipy.sparse import diags
from scipy.sparse.linalg import spsolve

if TYPE_CHECKING:
    from numpy.typing import NDArray
    from scipy.sparse import csr_array

class ScaledLinearSparseSolverConfig(BaseModel):
    """Configuration options for the :class:`ScaledLinearSparseSolver`.
    See the solver description for interpretation of the available options.
    """
    num_scale: int = Field(default=5, ge=0, le=10)
    """The maximum number of iterative scaling cycles"""
    dense_limit: int = Field(default=20_000, ge=0)
    """The maximum system size to even attempt ``numpy.linalg.solve``"""
    res_tol: float = Field(default=1e-6, gt=0.0, lt=1.0)
    """The tolerance of reproduced right hand-side accepted for the solution"""
    norm_tol: float = Field(default=1e-4, gt=0.0, lt=1.0)
    """The tolerance criteria to terminate scaling cycles"""


class ScaledLinearSparseSolver:
    r"""In the core, this solver utilizes the ``scipy.sparse.linalg.spsolve``
    function. However, as equation systems that represent thermodynamic
    systems are often badly scaled, the system is pre-scaled by linear
    transformations on columns and rows. This helps the pivoting heuristics
    of the solver to remain efficient and robust.

    By experience - a lot more without the scaling - pivoting heuristics
    in sparse solvers can easily inflate numerical imprecision for these
    matrices. The core solver will then not generate any error, but return an
    array that is far from the actual solution.

    The scaling applies at most
    :attr:`numscale <ScaledLinearSparseSolverConfig.num_scale>`
    iterative scaling of both rows and columns, such that

    .. math::

        A\,x = r \quad\Rightarrow\quad
        S_r\,A\,S_x\,S_x^{-1}\,x = S_r\,r \quad\Rightarrow\quad
        A' = S_r\,A\,S_x,\quad x' = S_x^{-1}\,x, \quad r' = S_r\,r

    with (:math:`\varepsilon_n =`
    :attr:`norm_tol <ScaledLinearSparseSolverConfig.norm_tol>`)

    .. math::

        |1 - \sum_i {a'}^2_{ij}| < \varepsilon_n\ \forall_j \quad\text{and}\quad
        |1 - \sum_j {a'}^2_{ij}| < \varepsilon_n\ \forall_i

    After solving, the residual condition
    :math:`||r' - A'\,x') / ||r'|| < \varepsilon_t =`
    :attr:`res_tol <ScaledLinearSparseSolverConfig.res_tol>` is evaluated.

    If the condition is not met, the dense
    ``numpy.linalg.solve`` function is tried as a fall-back up to a system
    size of :attr:`dense_limit <ScaledLinearSparseSolverConfig.dense_limit>`.

    This dense solver is normally robust unless the matrix truly is bad
    conditioned.
    """
    def __init__(self,
                 config: ScaledLinearSparseSolverConfig | None = None,
                 **options: Any):
        config = config or ScaledLinearSparseSolverConfig()
        self.config = config.model_copy(update=options)

    def reset(self):
        pass

    def solve(self, matrix: csr_array, rhs: NDArray) -> NDArray:
        s_x = 1
        # scale system if demanded
        if self.config.num_scale:
            matrix, s_r, s_x = self._scale(matrix)
            rhs = rhs / s_r

        if norm(rhs) == 0.0:
            return zeros_like(rhs)

        # solve with sparse solver
        x = spsolve(matrix, rhs)
        error = _residual(matrix, rhs, x)
        if error < self.config.res_tol:
            return x / s_x

        if rhs.shape[0] > self.config.dense_limit:
            msg = (f"Linear solver error, remaining error norm: {error:.3g}; "
                   f"System too large (N = {rhs.shape[0]})for dense solver")
            raise RuntimeError(msg)

        # retry with dense solver
        x = solve(matrix, rhs)
        error = _residual(matrix, rhs, x)
        if error < self.config.res_tol:
            return x / s_x

        msg = f"Linear solver error, remaining error norm: {error:.3g}"
        raise RuntimeError(msg)

    def _scale(self, matrix: csr_array):
        total_col_norms = 1
        total_row_norms = 1
        size = matrix.shape[0]
        for i in range(self.config.num_scale):
            col_norms = sqrt(matrix.power(2).sum(axis=0))
            col_norms[col_norms == 0] = 1.0
            total_col_norms *= col_norms
            matrix = matrix @ diags(1.0 / col_norms)

            row_norms = sqrt(matrix.power(2).sum(axis=1))
            row_norms[row_norms == 0] = 1.0
            total_row_norms *= row_norms
            matrix = (matrix.T @ diags(1.0 / row_norms)).T
            residual = 0.5 * (_scale_error(row_norms) + _scale_error(col_norms))
            if residual < self.config.norm_tol * size:
                break
        return matrix, total_row_norms, total_col_norms


class NumpySolver:
    r"""This solver is just a wrap around ``numpy.linalg.solve`` and as such
    suitable for systems smaller than :math:`\approx 10^3`.

    For such small systems, this can be faster that using sparse solvers.
    Due to rigorous pivoting, the solution can also be expected more numerically
    stable.

    Typical applications are small models used for evaluating samples of
    thermodynamic data sets.
    """
    def __init__(self, res_tol=1e-6):
        self.res_tol = res_tol

    def reset(self):
        pass

    def solve(self, matrix: csr_array, rhs: NDArray) -> NDArray:
        x = solve(matrix.todense(), rhs)
        error = _residual(matrix, rhs, x)
        if error > self.res_tol:
            msg = f"Linear solver error, remaining error norm: {error:.3g}"
            raise RuntimeError(msg)
        return x


def _residual(matrix: csr_array, rhs: NDArray, x: NDArray) -> float:
    return float(norm(rhs - matrix @ x) / norm(rhs))

def _scale_error(norms: NDArray) -> float:
    e = norms - 1
    return float(e @ e)