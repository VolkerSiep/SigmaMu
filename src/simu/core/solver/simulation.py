from symtable import Function
from typing import Callable, Sequence, Tuple
from copy import deepcopy
from dataclasses import dataclass, field
from time import time
from sys import stdout
from io import TextIOBase

from casadi import MX, jacobian, jtimes, Function
from numpy import array, argmin, argmax, abs, squeeze, log10
from scipy.sparse import csc_array
from pypardiso import spsolve

from ..model.numeric import NumericHandler
from ..utilities import Quantity, QFunction
from ..utilities.types import Map, NestedMutMap, NestedMap
from ..utilities.configurable import Configurable

_VEC, _STATE = NumericHandler.VECTORS, NumericHandler.STATE_VEC
_RES, _BOUND = NumericHandler.RES_VEC, NumericHandler.BOUND_VEC

# todo:
#  - move to utilities.table
#  - can I allow logging to the right logger additionally?
#    (just give it as a bloody optional argument, you dork!)
#  - document the whole solver thingy
class Table:
    def __init__(self, columns: Map[Tuple[str, str]], row_dig = 3,
                 row_head = "Row", stream=stdout):
        self.__cols = columns
        self._first = True
        self._row_fmt = None if row_dig is None else f"{{:{row_dig}d}}"
        self._row_head = row_head
        self._write = (lambda t: None) if (stream is None) else stream.write
        self._row = 1

    def row(self, data: object, row: int = None):
        write = self._write
        row = self._row if row is None else row

        elem = [c[1].format(getattr(data, k)) for k, c in self.__cols.items()]
        if self._row_fmt is not None:
            row = self._row_fmt.format(row)

        if self._first:
            headings = [f"{c[0]:>{len(e)}s}"
                        for e, c in zip(elem, self.__cols.values())]
            if self._row_fmt is not None:
                row_head = f"{self._row_head:{len(row)}s}"
                headings = [row_head] + headings

            write(" ".join(headings) + "\n")
            write(" ".join("-" * len(h) for h in headings) + "\n")
            self._first = False

        if self._row_fmt is not None:
            elem = [row] + elem
        write(" ".join(elem) + "\n")




@dataclass
class SimulationSolverIterationReport:
    max_err: float  # maximum error to tolerance in residuals
    max_res_name: str  # name of residual with maximum error
    relax_factor: float  # applied relaxation factor
    min_alpha_name: str  # name of most step-constraining bound variable
    duration: float  # accumulative seconds past during solving process
    lmet: float = field(init=False)  # logarithmic max error to tolerance

    def __post_init__(self):
        self.lmet = log10(self.max_err + 1e-8)

@dataclass
class SimulationSolverReport:
    iterations: Sequence[SimulationSolverIterationReport]
    final_state: Sequence[float]
    result: NestedMap[Quantity]


CALL_BACK_TYPE = Callable[
    [int, SimulationSolverIterationReport, Quantity, QFunction], bool]


class SimulationSolver(Configurable):
    r"""
    The simulation solver assumes both thermodynamic and model parameters to
    be constant
    """
    # noinspection PyUnusedLocal
    # Options are parsed via inspection
    def __init__(self, model: NumericHandler, *,
                 max_iter: int = 30,
                 gamma: float = 0.9,
                 wall: float = 1e-20,
                 output: TextIOBase|None = stdout,
                 call_back_iter: CALL_BACK_TYPE = None):
        """

        """
        super().__init__(exclude=["model"])
        self._model = model

        # store arguments (parameters) so the user can change them
        args = deepcopy(model.arguments)
        # store size of state
        self.__state_size = args[_VEC][_STATE].magnitude.size()[0]
        # user shall not think that putting a state here has any effect
        del args[_VEC][_STATE]
        self.__model_parameters : NestedMutMap[Quantity] = args

    def solve(self) -> SimulationSolverReport:
        start_time = time()
        opt = self.options
        model = self._model
        residual_names = model.vector_res_names(_RES)
        bound_names = model.vector_res_names(_BOUND)
        reports = []

        table = Table({
            "lmet": ("LMET", "{:5.1f}"),
            "relax_factor": ("Alpha", "{:7.2g}"),
            "duration": ("Time", "{:6.1g}"),
            "min_alpha_name": ("Limit on bound", "{:>40s}"),
            "max_res_name": ("Max residual", "{:>40s}")
        }, row_dig=5, row_head="Iter", stream=opt["output"])

        funcs = self._prepare_functions()
        x = self.initial_state

        for iteration in range(opt["max_iter"]):
            # evaluate system (matrix and rhs)
            r, dr_dx = funcs["f_r"](x)
            r = squeeze(array(r))
            dr_dx = csc_array(dr_dx)

            # assess error
            max_err_idx = argmax(abs(r))
            max_res_name = residual_names[max_err_idx]
            max_err = abs(r[max_err_idx])
            if max_err < 1:
                break

            # calculate full update
            dx = -spsolve(dr_dx, r)

            # find relaxation factor
            a = squeeze(array(funcs["f_b"](x, dx)))
            a = a[0 < a]
            alpha, min_alpha_name = 1, ""
            if len(a):
                min_a_idx = int(argmin(a))
                if a[min_a_idx] * opt["gamma"] < 1:
                    alpha = a[min_a_idx] * opt["gamma"]
                    min_alpha_name = bound_names[min_a_idx]
                if alpha < opt["wall"]:
                    msg = f"Relaxation coefficient is below {opt["wall"]}, " \
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
                    iteration, reports[-1], x,
                    lambda x_arg: funcs["f_y"]({"x": Quantity(x_arg)})
                )
                if not cb_result:
                    msg = "Solver iterations interrupted by callback"
                    raise ValueError(msg)
            table.row(reports[-1], iteration)
        else:
            msg = f"Model did not converge after {opt["max_iter"]} iterations"
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
        return SimulationSolverReport(
            iterations=reports,
            final_state=x,
            result=funcs["f_y"]({"x": Quantity(x)})
        )

    def _prepare_functions(self) -> Map[Callable]:
        # prepare
        #  - a casadi MX function x -> (r, dr/dx)
        #  - a casadi MX function: (x, dx) -> (a_i = b_i / (db_i/dx_j) * dx_j)
        # prepare a QFunction x -> (y_m, y_t)
        param = deepcopy(self.__model_parameters)

        param[_VEC][_STATE] = (Quantity(x := MX.sym("x", self.__state_size)))
        res = self._model.function(param, squeeze_results=False)
        r, b = res[_VEC][_RES], res[_VEC][_BOUND]
        dx = Quantity(MX.sym("x", self.__state_size))
        return {
            "f_r": Function("f_r", [x], [r, jacobian(r, x)]),
            "f_b": Function("f_b", [x, dx], [-b / jtimes(b, x, dx)]),
            "f_y": QFunction({"x": Quantity(x)}, res)
        }

    @property
    def initial_state(self):
        """Freshly extract the initial values from the model. These might have
        been changed after the solver class was instantiated"""
        args = self._model.arguments
        return args[NumericHandler.VECTORS][NumericHandler.STATE_VEC]

    @property
    def model_parameters(self) -> NestedMutMap[Quantity]:
        return self.__model_parameters

    @property
    def _arg_validations(self):
        between = Configurable._validate_between
        return {
            "max_iter": between(1, 10000),
            "gamma": between(0.1, 0.999),
            "call_back_iter": {
                "f": lambda x: x is None or callable(x),
                "msg": "must be callable",
            },
            "output": {
                "f": lambda x: x is None or isinstance(x, TextIOBase),
                "msg": "must be a stream"
            }
        }
