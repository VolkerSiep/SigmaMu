# stdlib
from typing import Any
from time import time
from collections.abc import Generator, Sequence, Callable
from dataclasses import dataclass

# external
from casadi import jacobian, jtimes
from numpy import array, squeeze, zeros_like
from numpy.typing import NDArray
from scipy.sparse import csr_array

from simu import SymbolQuantity
# internal
from simu.core.model.numeric import NumericHandler, NHKeys
from simu.core.utilities.quantity import Quantity
from simu.core.utilities.qfunction import QFunction
from simu.core.utilities.output import ProgressTableOutput
from simu.core.utilities.types import Map, NestedMap, NestedMutMap
from simu.core.utilities.errors import IterativeProcessInterrupted

from ..common import relax, not_finite, assess_residuals, check_model_square
from .config import SimulationSolverConfig
from .report import (SimulationSolverReport, SimulationSolverIterationReport)
from ...utilities.structures import FlatDictionary, path_to_flat_key

type SolverFunction = Callable[[FlatDictionary], NestedMutMap[Quantity]]

class _SolverFunctionWrapper:
    def __init__(self, func: QFunction, keys: Sequence[str]):
        self._func = func
        self._keys = keys

    def __call__(self, args: FlatDictionary) -> Sequence[NDArray]:
        result = self._func.evaluate_prepared(args)
        return [result[k].m for k in self._keys]


@dataclass
class _FunctionCollection:
    f_r: _SolverFunctionWrapper  # args -> r, r_x
    f_b: _SolverFunctionWrapper  # args -> a, b
    f_y: SolverFunction  # args -> y


_OUTPUT_TABLE_DEFINITION = {
    "iteration": ("Iter", "{: 5d}"),
    "lmet": ("LMET", "{:5.1f}"),
    "relax_factor": ("Alpha", "{:7.2g}"),
    "duration": ("Time", "{:6.2f}"),
    "min_alpha_name": ("Limit on bound", "{:>50s}"),
    "max_res_name": ("Max residual", "{:>50s}")
}


class SimulationSolver:
    """
    The simulation solver assumes both thermodynamic and model parameters to
    be constant, aiming to find the state variable values such that all
    residuals evaluate to zero within their tolerance.
    """

    def __init__(self, model: NumericHandler,
                 config: SimulationSolverConfig | None = None, **options: Any):
        """On construction, the solver object requires a
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
        self._config = (config or SimulationSolverConfig()).update(**options)

        check_model_square(model)

        # define functions and prepare arguments for efficiency
        self._funcs, self._model_parameters = self._prepare_functions()
        self._x : NDArray | None = None

    def set_options(self, config: SimulationSolverConfig | None = None,
                    **options: Any):
        """Overwrite configuration for subsequent solver runs

       :param config: Options for the solver as defined in
          :class:`~simu.core.solver.simulation.config.SimulationSolverConfig`.
       :param options: overwriting individual configurations directly
       """
        self._config = (config or self._config).update(**options)

    def solve(self, **options: Any) -> SimulationSolverReport:
        """
        This method triggers iterative the solving process. This takes less than
        0.1 seconds for small models, and increases with model complexity.
        For large models, the time per iteration is due to the solving of linear
        systems cubic in system size, though the model structure might render
        this a conservative estimate.

        :param options: overwriting individual configurations for this solver
          run.

        :return: The report including the iteration sequence
        """
        config = self._config.update(**options)
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
                assert self._x is not None
                assert self._funcs is not None
                cb_result = config.call_back_iter(
                    iter_report, array(self._x),
                    lambda: self._funcs.f_y(self._model_parameters)
                )
                if not cb_result:
                    msg = "Solver iterations interrupted by callback"
                    raise IterativeProcessInterrupted(msg)

        # retain state if desired
        if config.retain_solution:
            thermo_param = self.model_parameters[NHKeys.THERMO_PARAMS]
            self._model.retain_state(array(self._x), thermo_param)

        assert self._funcs is not None
        return SimulationSolverReport(
            iterations=reports,
            final_state=array(self._x),
            prop_func=lambda: self._funcs.f_y(self._model_parameters)
        )

    def solve_iter(
            self, config: SimulationSolverConfig | None = None,
            **options: Any
    ) -> Generator[SimulationSolverIterationReport, Map[Any] | None, None]:
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
        config = (config or self._config).update(**options)
        start_time = time()
        residual_names = model.vector_res_names(NHKeys.RESIDUALS)
        bound_names = model.vector_res_names(NHKeys.BOUNDS)
        funcs, args = self._funcs, self._model_parameters
        x_key = path_to_flat_key([NHKeys.VECTORS, NHKeys.STATES])
        dx_key = path_to_flat_key([NHKeys.VECTORS, NHKeys.DELTA_STATES])

        self._x = x = self.initial_state

        config.linear_solver.reset()

        for iteration in range(config.max_iter):
            # evaluate system (matrix and rhs)
            args[x_key] = x
            r, dr_dx = funcs.f_r(args)
            r = squeeze(array(r))
            dr_dx = csr_array(dr_dx)

            max_err, max_res_name = assess_residuals(r, residual_names)
            if max_err < 1:
                break

            # calculate full update
            dx = -config.linear_solver.solve(dr_dx, r)
            args[dx_key] = Quantity(dx)
            # find relaxation factor
            b, a = funcs.f_b(args)
            alpha, min_alpha_name = relax(array(b), array(a),
                                          bound_names, config.gamma)
            if alpha < config.wall:
                msg = f"Relaxation factor is below {config.wall}, " \
                      "no solution found"
                raise ValueError(msg)

            # apply update
            self._x = x = x + alpha * dx

            # reporting
            duration = time() - start_time
            options = yield SimulationSolverIterationReport(
                iteration=iteration,
                max_err=float(max_err),
                max_res_name=max_res_name,
                relax_factor=float(alpha),
                min_alpha_name=min_alpha_name,
                duration=duration
            )
            if options is not None:
                config = config.update(**options)
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

    def _prepare_functions(self) -> tuple[_FunctionCollection, FlatDictionary]:
        # TODO:
        #   Under construction:
        #   - extend args with DX as an argument vector
        #   - create QFunction objects that take all arguments, and return
        #     what the three functions returned before.
        #   - prepare the arguments using self.model_parameters
        #     I do not need to remove the state anymore, as the object is not
        #     exposed the the user anymore.
        #   When updating the parameters:
        #   - allow only single model parameter updates at a time.
        #   - This looks up the unit in the QFunction, converts, and inserts
        #     the magnitude into the prepared structure.
        #   During iterations:
        #   - just patch state in prepared parameters.
        #   When reporting:
        #   - return lambda function that calls f_y with final parameters, no
        #     meaning to require even the state as an argument, hence f_y has
        #     no parameters

        args = construct_extended_arguments(self._model)
        res = self._model.sym_res
        vectors = res[NHKeys.VECTORS]
        x = args[NHKeys.VECTORS][NHKeys.STATES]
        dx = args[NHKeys.VECTORS][NHKeys.DELTA_STATES]
        r, b = vectors[NHKeys.RESIDUALS], vectors[NHKeys.BOUNDS]
        jac = Quantity(jacobian(r.m, x.m))
        a = Quantity(-b.m / jtimes(b.m, x.m, dx.m))
        f_r = QFunction(args, {"r": r, "dr_dx": jac})
        f_b = QFunction(args, {"b": b, "a": a})
        f_y = QFunction(args, res)

        args = self._model.arguments
        dx = zeros_like(args[NHKeys.VECTORS][NHKeys.STATES].m)
        args[NHKeys.VECTORS][NHKeys.DELTA_STATES] = Quantity(dx)

        functions = _FunctionCollection(
            f_r=_SolverFunctionWrapper(f_r, ["r", "dr_dx"]),
            f_b=_SolverFunctionWrapper(f_b, ["b", "a"]),
            f_y=f_y.evaluate_prepared
        )
        model_parameters = f_r.prepare_arguments(args)
        return functions, model_parameters

    @property
    def initial_state(self) -> NDArray:
        """Freshly extract the initial values from the model. These might have
        been changed after the solver class was instantiated"""
        args = self._model.arguments
        return args[NHKeys.VECTORS][NHKeys.STATES].magnitude

    @property
    def model_parameters(self) -> NestedMutMap[Quantity]:
        """A convenience property to access the parameters of the model as a
        mutable object. The state variables are removed in this instance, as
        these are rather provided by the solver during the iterative solving
        process."""
        return self._model_parameters


def assess_max_res(r: NDArray, names: Sequence[str]):
    if not_finite_res := not_finite(r, names):
        nf = ", ".join(not_finite_res)
        msg = f"Non-finite values in the following residuals: {nf}"
        raise ValueError(msg)
    return assess_residuals(r, names)

def construct_extended_arguments(model: NumericHandler) -> NestedMap[Quantity]:
    args = model.sym_args
    size = len(model.vector_arg_names(NHKeys.STATES))
    vectors = args[NHKeys.VECTORS]
    dx_sym = SymbolQuantity(
        "dx", "dimensionless",
        ["dx_{i:05d}" for i in range(size)]
    )
    return {
        NHKeys.MODEL_PARAMS: args[NHKeys.MODEL_PARAMS],
        NHKeys.THERMO_PARAMS: args[NHKeys.THERMO_PARAMS],
        NHKeys.VECTORS: {
           NHKeys.STATES: vectors[NHKeys.STATES],
            NHKeys.DELTA_STATES: dx_sym,
        }
    }