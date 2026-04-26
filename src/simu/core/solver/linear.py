from __future__ import annotations
from typing import TYPE_CHECKING, Any
from pydantic import BaseModel, Field
from numpy import sqrt
from numpy.linalg import norm, solve
from scipy.sparse import diags
from scipy.sparse.linalg import spsolve

if TYPE_CHECKING:
    from numpy.typing import NDArray
    from scipy.sparse import csr_array

class ScaledLinearSparseSolverConfig(BaseModel):
    num_scale: int = Field(default=5, ge=0, le=10)
    dense_limit: int = Field(default=20_000, ge=0)
    res_tol: float = Field(default=1e-6, gt=0.0, lt=1.0)
    norm_tol: float = Field(default=1e-4, gt=0.0, lt=1.0)


class ScaledLinearSparseSolver:
    def __init__(self, config: ScaledLinearSparseSolverConfig, **options: Any):
        config = config or ScaledLinearSparseSolverConfig()
        self.config = config.model_copy(update=options)

    def solve(self, matrix: csr_array, rhs: NDArray) -> NDArray:
        s_x = None
        # scale system if demanded
        if self.config.num_scale:
            matrix, s_r, s_x = self._scale(matrix)
            rhs = rhs @ s_r

        # solve with sparse solver
        x = spsolve(matrix, rhs)
        error = norm(rhs - matrix @ x) / norm(rhs)
        if error < self.config.res_tol:
            return x if s_x is None else x @ s_x

        if rhs.shape[0] > self.config.dense_limit:
            msg = (f"Linear solver error, remaining error norm: {error:.3g}; "
                   f"System too large (N = {rhs.shape[0]})for dense solver")
            raise RuntimeError(msg)

        # retry with dense solver
        x = solve(matrix, rhs)
        error = norm(rhs - matrix @ x) / norm(rhs)
        if error < self.config.res_tol:
            return x if s_x is None else x @ s_x

        msg = f"Linear solver error, remaining error norm: {error:.3g}"
        raise RuntimeError(msg)

    def _scale(self, matrix: csr_array):
        total_col_norms = 1
        total_row_norms = 1
        for i in range(self.config.num_scale):
            col_norms = sqrt(matrix.multiply(matrix).sum(axis=0))
            col_norms[col_norms == 0] = 1.0
            total_col_norms = total_col_norms * col_norms
            matrix = matrix @ diags(1.0 / col_norms)

            row_norms = sqrt(matrix.multiply(matrix).sum(axis=1))
            row_norms[row_norms == 0] = 1.0
            total_row_norms = total_row_norms * row_norms
            matrix = (matrix.T @ diags(1.0 / row_norms)).T
            qlt = (sum((row_norms - 1) ** 2) +
                   sum((col_norms - 1) ** 2)) / matrix.shape[0]
            if qlt < self.config.norm_tol:
                break
        return matrix, diags(1.0 / total_row_norms), diags(1.0 / total_col_norms)