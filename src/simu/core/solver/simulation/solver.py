# stdlib
from symtable import Function
from typing import Any
from copy import deepcopy
from time import time
from collections.abc import Sequence, Iterator
from dataclasses import dataclass

# external
from casadi import SX, jacobian, jtimes, Function
from numpy import array, argmin, argmax, abs, squeeze, isfinite
from numpy.typing import NDArray
from scipy.sparse import csr_array

# internal
from simu.core.model.numeric import NumericHandler, NHKeys
from simu.core.utilities.quantity import Quantity, QFunction
from simu.core.utilities.output import ProgressTableOutput
from simu.core.utilities.types import NestedMutMap
from simu.core.utilities.errors import (
    IterativeProcessInterrupted, NonSquareSystem)

from .report import (
    SimulationSolverReport, SimulationSolverIterationReport, PropertyFunction)
from .config import SimulationSolverConfig

@dataclass
class _FunctionCollection:
    f_r: Function
    f_b: Function
    f_y: PropertyFunction


_OUTPUT_TABLE_DEFINITION = {
    "iteration": ("Iter", "{: 5d}"),
    "lmet": ("LMET", "{:5.1f}"),
    "relax_factor": ("Alpha", "{:7.2g}"),
    "duration": ("Time", "{:6.2f}"),
    "min_alpha_name": ("Limit on bound", "{:>50s}"),
    "max_res_name": ("Max residual", "{:>50s}")
}

class SimulationSolver:
    r"""
    The simulation solver assumes both thermodynamic and model parameters to
    be constant, aiming to find the state variable values such that all
    residuals evaluate to zero within their tolerance.
    """

    def __init__(self, model: NumericHandler,
                 config: SimulationSolverConfig | None = None, **options: Any):
        r"""On construction, the solver object requires a
        :class:`~simu.NumericHandler` object. The solver object can then be
        reused for multiple solver runs, for instance with variable parameter
        values (sensitivity study).

        :param model: The numeric handler of a model, in most cases obtained by
          the expression ``NumericHandler(ModelClass.top())``.
        :param config: Options for the solver as defined in
          :class:`~simu.core.solver.simulation.config.SimulationSolverConfig`.
        :param options: overwriting individual configurations directly
        """
        self._model = model
        self._config = config or SimulationSolverConfig()
        self.set_options(**options)

        args = model.arguments
        # store size of state
        self.__state_size = args[NHKeys.VECTORS][NHKeys.STATES].m.size()[0]
        res_size = len(model.vector_res_names(NHKeys.RESIDUALS))

        if self.__state_size != res_size:
            raise NonSquareSystem(self.__state_size, res_size)

        # user shall not think that putting a state here has any effect
        del args[NHKeys.VECTORS][NHKeys.STATES]
        self._model_parameters: NestedMutMap[Quantity] = args

        self._funcs : _FunctionCollection | None = None
        self._x : Quantity | None = None

    def set_options(self, config: SimulationSolverConfig | None = None,
                    **options: Any):
        """Overwrite configuration for subsequent solver runs

       :param config: Options for the solver as defined in
          :class:`~simu.core.solver.simulation.config.SimulationSolverConfig`.
       :param options: overwriting individual configurations directly
       """
        config = config or self._config
        self._config = config.model_copy(update=options)

    def solve(self, **options: Any) -> SimulationSolverReport:
        """
        This method triggers iterative the solving process. This takes less than
        0.1 seconds for small models, and increases with model complexity.
        For large models, the time per iteration is due to the solving of linear
        systems cubic in system size, though the model structure might render
        this a conservative estimate.

        With `pypardiso`_ installed, the solving of the linear systems is
        performed on all available CPU cores. However, their solver sometimes
        chokes and returns a wrong solution. Therefore, the norm of the
        solution is checked, and ``scipy.sparse.linalg.spsolve`` is used in
        those instances.

        :param options: overwriting individual configurations for this solver
          run.

        :return: The report including the iteration sequence
        """
        config = self._config.model_copy(update=options)
        table = ProgressTableOutput(
            _OUTPUT_TABLE_DEFINITION,
            output=config.output
        )
        reports = []
        for iter_report in self.solve_iter(config):
            reports.append(iter_report)
            table.row(iter_report)

            # callback
            if config.call_back_iter is not None:
                cb_result = config.call_back_iter(
                    iter_report, self._x.magnitude, self._funcs.f_y
                )
                if not cb_result:
                    msg = "Solver iterations interrupted by callback"
                    raise IterativeProcessInterrupted(msg)

        # retain state if desired
        if self._config.retain_solution:
            thermo_param = self.model_parameters[NHKeys.THERMO_PARAMS]
            self._model.retain_state(self._x.nonzeros(), thermo_param)

        return SimulationSolverReport(
            iterations=reports,
            final_state=self._x,
            prop_func=self._funcs.f_y
        )

    def solve_iter(self, config: SimulationSolverConfig = None
                   ) -> Iterator[SimulationSolverIterationReport]:
        """Run individual iterations and return control flow back to the client
        code after each iteration. This allows for finer control in a
        multithreaded environment, for instance to update a GUI with trends
        about the convergence progress and intermediate values, or to allow
        interactive pausing and cancelling of the simulation run.

        :param config: An optional opportunity to overwrite solver options.
          Note that

          - ``output`` will be ignored, as no output is written
          - ``call_back_iter`` will be ignored, as no callback will be called

        :return: An iterator over all generated iteration reports.
        """
        model = self._model
        config = config or SimulationSolverConfig()
        start_time = time()
        residual_names = model.vector_res_names(NHKeys.RESIDUALS)
        bound_names = model.vector_res_names(NHKeys.BOUNDS)

        self._funcs = funcs = self._prepare_functions()
        self._x = x = self.initial_state

        for iteration in range(config.max_iter):
            # evaluate system (matrix and rhs)
            r, dr_dx = funcs.f_r(x)
            r = squeeze(array(r))
            dr_dx = csr_array(dr_dx)

            if not_final := _not_final(r, residual_names):
                nf = ", ".join(not_final)
                msg = f"Non-finite values in the following residuals: {nf}"
                raise ValueError(msg)

            max_err, max_res_name = _assess_residuals(r, residual_names)
            if max_err < 1:
                break

            # calculate full update
            dx = -self._config.linear_solver.solve(dr_dx, r)

            # find relaxation factor
            b, a = funcs.f_b(x, dx)
            alpha, min_alpha_name = self._relax(b, a, bound_names)
            if alpha < config.wall:
                msg = f"Relaxation factor is below {config.wall}, " \
                      "no solution found"
                raise ValueError(msg)

            # apply update
            self._x = x = x + alpha * dx

            # reporting
            duration = time() - start_time
            report = SimulationSolverIterationReport(
                iteration=iteration,
                max_err=float(max_err),
                max_res_name=max_res_name,
                relax_factor=float(alpha),
                min_alpha_name=min_alpha_name,
                duration=duration
            )
            yield report
        else:
            msg = f"Model did not converge after {config.max_iter} iterations"
            raise ValueError(msg)

        # reporting
        duration = time() - start_time
        yield SimulationSolverIterationReport(
            iteration=iteration,
            max_err=float(max_err),
            max_res_name=max_res_name,
            relax_factor=1,
            min_alpha_name="",
            duration=duration
        )

    def _prepare_functions(self) -> _FunctionCollection:
        # prepare
        #  - a casadi function x -> (r, dr/dx)
        #  - a casadi function: (x, dx) -> (a_i = b_i / (db_i/dx_j) * dx_j)
        # prepare a QFunction x -> (y_m, y_t)
        param = deepcopy(self._model_parameters)
        x = SX.sym("x", self.__state_size)
        param[NHKeys.VECTORS][NHKeys.STATES] = Quantity(x)
        res = self._model.function(param, squeeze_results=False)  # EXPENSIVE!!
        vectors = res[NHKeys.VECTORS]
        r, b = vectors[NHKeys.RESIDUALS].m, vectors[NHKeys.BOUNDS].m
        dx = SX.sym("dx", self.__state_size)
        f_y = QFunction({"x": Quantity(x)}, res)  # EXPENSIVE!!

        return _FunctionCollection(
            f_r=Function("f_r", [x], [r, jacobian(r, x)]),
            f_b=Function("f_b", [x, dx], [b, -b / jtimes(b, x, dx)]),
            f_y=lambda z: f_y({"x": Quantity(z)})
        )

    def _relax(self, b: NDArray, a: NDArray,
               bound_names: Sequence[str]) -> tuple[float, str]:
        config = self._config
        a, b = [squeeze(array(x)) for x in (a, b)]
        # are there bounds violated?
        invalid = [n for n, m_i in zip(bound_names, b <= 0) if m_i]

        if invalid:
            msg = f"Bound violation of: {', '.join(invalid)}"
            raise ValueError(msg)

        mask = (a > 0)
        a = a[mask]
        alpha, min_alpha_name = 1.0, ""
        if not len(a):
            return alpha, min_alpha_name

        min_a_idx = int(argmin(a))
        if a[min_a_idx] * config.gamma < 1:
            alpha = a[min_a_idx] * config.gamma
            bn = [b for b, m in zip(bound_names, mask) if m]
            min_alpha_name = bn[min_a_idx]
        return alpha, min_alpha_name

    @property
    def initial_state(self):
        """Freshly extract the initial values from the model. These might have
        been changed after the solver class was instantiated"""
        args = self._model.arguments
        return args[NHKeys.VECTORS][NHKeys.STATES]

    @property
    def model_parameters(self) -> NestedMutMap[Quantity]:
        """A convenience property to access the parameters of the model as a
        mutable object. The state variables are removed in this instance, as
        these are rather provided by the solver during the iterative solving
        process."""
        return self._model_parameters


def _not_final(vector: NDArray, names: Sequence[str]) -> Sequence[str]:
    finite = isfinite(vector)
    if False in isfinite(finite):
        return [names[i] for i, f in enumerate(finite) if not f]
    return []

def _assess_residuals(vector: NDArray,
                      names: Sequence[str]) -> tuple[float, str]:
    if len(vector):
        # assess error
        max_err_idx = int(argmax(abs(vector)))
        max_name = names[max_err_idx]
        max_err = float(abs(vector[max_err_idx]))
        return max_err, max_name
    else:  # trivial model, nothing to solve
        return 0, ""
