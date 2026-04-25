# stdlib
import sys
from symtable import Function
from typing import Callable, Any
from copy import deepcopy
from time import time
from io import TextIOBase

# external
from casadi import SX, jacobian, jtimes, Function
from numpy import array, argmin, argmax, abs, squeeze, isfinite, sqrt
from numpy.linalg import norm, solve
from scipy.sparse import csr_array, diags
from scipy.sparse.linalg import spsolve as scipy_spsolve

try:  # use pypardiso if installed
    from pypardiso import spsolve
except ImportError:  # use scipy if not
    spsolve = scipy_spsolve

# internal
from simu.core.model.numeric import NumericHandler, NHKeys
from simu.core.utilities.quantity import Quantity, QFunction
from simu.core.utilities.output import ProgressTableOutput
from simu.core.utilities.types import Map, NestedMutMap
from simu.core.utilities.configurable import Configurable
from simu.core.utilities.errors import (
    IterativeProcessInterrupted, NonSquareSystem)

from .report import SimulationSolverReport, SimulationSolverIterationReport
from .config import SimulationSolverCallback


class SimulationSolver(Configurable):
    r"""
    The simulation solver assumes both thermodynamic and model parameters to
    be constant, aiming to find the state variable values such that all
    residuals evaluate to zero within their tolerance.
    """

    # noinspection PyUnusedLocal
    # Options are parsed via inspection
    def __init__(self, model: NumericHandler, *,
                 max_iter: int = 30,
                 gamma: float = 0.9,
                 wall: float = 1e-20,
                 output: TextIOBase | str = "stdout",
                 call_back_iter: SimulationSolverCallback = None,
                 retain_solution: bool = True):
        r"""On construction, the solver object requires a
        :class:`~simu.NumericHandler` object. The solver object can then be
        reused for multiple solver runs, for instance with variable parameter
        values (sensitivity study).

        :param model: The numeric handler of a model, in most cases obtained by
          the expression ``NumericHandler(ModelClass.top())``.
        :param kwargs: Options for the solver as defined in
          :class:`~simu.core.solver.simulation.config.SimulationSolverConfig`.
        """
        super().__init__(exclude=["model"])
        self._model = model

        # store arguments (parameters) so the user can change them
        args = deepcopy(model.arguments)
        # store size of state
        self.__state_size = args[NHKeys.VECTORS][NHKeys.STATES].m.size()[0]
        res_size = len(model.vector_res_names(NHKeys.RESIDUALS))

        if self.__state_size != res_size:
            raise NonSquareSystem(self.__state_size, res_size)

        # user shall not think that putting a state here has any effect
        del args[NHKeys.VECTORS][NHKeys.STATES]
        self.__model_parameters: NestedMutMap[Quantity] = args

    def solve(self, **kwargs: Any) -> SimulationSolverReport:
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

        For each given keyword argument, the
        :meth:`~simu.core.utilities.configurable.Configurable.set_option`
        method is invoked, giving the same effect as if the option was provided
        with the constructor.

        :return: The report including the iteration sequence
        """
        for name, value in kwargs.items():
            self.set_option(name, value)
        opt = self.options
        model = self._model
        start_time = time()
        residual_names = model.vector_res_names(NHKeys.RESIDUALS)
        bound_names = model.vector_res_names(NHKeys.BOUNDS)
        reports = []

        output = self.__find_output()

        table = ProgressTableOutput({
            "lmet": ("LMET", "{:5.1f}"),
            "relax_factor": ("Alpha", "{:7.2g}"),
            "duration": ("Time", "{:6.2f}"),
            "min_alpha_name": ("Limit on bound", "{:>50s}"),
            "max_res_name": ("Max residual", "{:>50s}")
        }, row_dig=5, row_head="Iter", stream=output)

        funcs = self._prepare_functions()
        x = self.initial_state

        for iteration in range(opt["max_iter"]):
            # evaluate system (matrix and rhs)
            r, dr_dx = funcs["f_r"](x)
            r = squeeze(array(r))
            r_finite = isfinite(r)
            if False in isfinite(r):
                names = [residual_names[i]
                         for i, f in enumerate(r_finite) if not f]
                nf = ", ".join(names)
                msg = f"Non-finite values in the following residuals: {nf}"
                raise ValueError(msg)

            dr_dx = csr_array(dr_dx)

            if len(r):
                # assess error
                max_err_idx = int(argmax(abs(r)))
                max_res_name = residual_names[max_err_idx]
                max_err = abs(r[max_err_idx])
                if max_err < 1:
                    break
            else:  # trivial model, nothing to solve
                max_err = 0
                max_res_name = ""
                break

            # calculate full update
            dx = self._solve_linear(dr_dx, r)

            # find relaxation factor
            b, a = map(lambda z: squeeze(array(z)), funcs["f_b"](x, dx))
            # are there bounds violated?
            invalid = [n for n, m_i in zip(bound_names, b <= 0) if m_i]
            if invalid:
                msg = f"Bound violation of: {', '.join(invalid)}"
                raise ValueError(msg)

            mask = (a > 0)
            a = a[mask]
            alpha, min_alpha_name = 1, ""
            if len(a):
                min_a_idx = int(argmin(a))
                if a[min_a_idx] * opt["gamma"] < 1:
                    alpha = a[min_a_idx] * opt["gamma"]
                    bn = [b for b, m in zip(bound_names, mask) if m]
                    min_alpha_name = bn[min_a_idx]
                if alpha < opt["wall"]:
                    msg = f"Relaxation factor is below {opt['wall']}, " \
                          "no solution found"
                    raise ValueError(msg)
            # apply update
            x = x + alpha * dx

            # reporting
            duration = time() - start_time
            reports.append(SimulationSolverIterationReport(
                max_err=float(max_err),
                max_res_name=max_res_name,
                relax_factor=float(alpha),
                min_alpha_name=min_alpha_name,
                duration=duration
            ))
            if opt["call_back_iter"] is not None:
                cb_result = opt["call_back_iter"](
                    iteration, reports[-1], x.magnitude,
                    lambda x_arg: funcs["f_y"]({"x": Quantity(x_arg)})
                )
                if not cb_result:
                    msg = "Solver iterations interrupted by callback"
                    raise IterativeProcessInterrupted(msg)
            table.row(reports[-1], iteration)
        else:
            msg = f"Model did not converge after {opt['max_iter']} iterations"
            raise ValueError(msg)

        # reporting
        duration = time() - start_time
        reports.append(SimulationSolverIterationReport(
            max_err=float(max_err),
            max_res_name=max_res_name,
            relax_factor=1,
            min_alpha_name="",
            duration=duration
        ))
        table.row(reports[-1], iteration)

        # retain state if desired
        if opt["retain_solution"]:
            thermo_param = self.model_parameters["thermo_params"]
            model.retain_state(x.nonzeros(), thermo_param)

        return SimulationSolverReport(
            iterations=reports,
            final_state=x,
            prop_func=lambda z: funcs["f_y"]({"x": Quantity(z)})
        )

    def __find_output(self) -> TextIOBase | None:
        output = self.options["output"]
        if isinstance(output, str):
            opts = {"stdout": sys.stdout, "none": None}
            try:
                return opts[output.lower()]
            except KeyError:
                raise ValueError(f"Invalid stream name '{output}'")
        return output

    def _prepare_functions(self) -> Map[Callable]:
        # prepare
        #  - a casadi function x -> (r, dr/dx)
        #  - a casadi function: (x, dx) -> (a_i = b_i / (db_i/dx_j) * dx_j)
        # prepare a QFunction x -> (y_m, y_t)
        param = deepcopy(self.__model_parameters)
        x = SX.sym("x", self.__state_size)
        param[NHKeys.VECTORS][NHKeys.STATES] = Quantity(x)
        res = self._model.function(param, squeeze_results=False)  # EXPENSIVE!!
        vectors = res[NHKeys.VECTORS]
        r, b = vectors[NHKeys.RESIDUALS].m, vectors[NHKeys.BOUNDS].m
        dx = SX.sym("dx", self.__state_size)
        f_y = QFunction({"x": Quantity(x)}, res)  # EXPENSIVE!!
        return {
            "f_r": Function("f_r", [x], [r, jacobian(r, x)]),
            "f_b": Function("f_b", [x, dx], [b, -b / jtimes(b, x, dx)]),
            "f_y": f_y
        }

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
        return self.__model_parameters

    @property
    def _arg_validations(self):
        between = Configurable._validate_between
        return {
            "max_iter": between(1, 10000),
            "gamma": between(0, 0.999),
            "call_back_iter": {
                "f": lambda x: x is None or callable(x),
                "msg": "must be callable",
            },
            "output": {
                "f": lambda x: isinstance(x, str) or isinstance(x, TextIOBase),
                "msg": "must be a stream or a qualified string"
            }
        }

    @staticmethod
    def _scale(matrix: csr_array, num=10):
        total_col_norms = 1
        total_row_norms = 1
        for i in range(num):
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
            if qlt < 0.001:
                break
        return matrix, diags(1.0 / total_row_norms), diags(1.0 / total_col_norms)

    @staticmethod
    def _solve_linear(dr_dx: csr_array, r):
        dr_dx, s_r, s_x = SimulationSolver._scale(dr_dx, num=5)
        n = r.shape[0]
        r = r @ s_r
        dx = -spsolve(dr_dx, r)
        dr = r + dr_dx @ dx
        if (nr := norm(dr)) > 0.1 * dr.shape[0]:
            # print(f"Remaining norm: {nr:.2f} - using scipy fallback")
            dx = -scipy_spsolve(dr_dx, r)
            dr = r + dr_dx @ dx
        if (nr := norm(dr)) > 0.1 * dr.shape[0]:
            if n < 1000:
                # print(f"Remaining norm: {nr:.2f} - using numpy fallback")
                return -solve(dr_dx.toarray(), r) @ s_x
            msg = f"Linear solver error, remaining residual: {nr:.2f}"
            raise ValueError(msg)
        return dx @ s_x
