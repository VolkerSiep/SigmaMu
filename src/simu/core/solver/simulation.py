from symtable import Function
from typing import Callable, Sequence, Any
from copy import deepcopy
from dataclasses import dataclass, field
from time import time
import sys
from io import TextIOBase

from casadi import MX, jacobian, jtimes, Function
from numpy import array, argmin, argmax, abs, squeeze, log10
from scipy.sparse import csc_array

try:  # use pypardiso if installed
    from pypardiso import spsolve
except ImportError:  # use scipy if not
    from scipy.sparse.linalg import spsolve

from ..model.numeric import NumericHandler
from ..utilities import Quantity, QFunction
from ..utilities.output import ProgressTableOutput
from ..utilities.types import Map, NestedMutMap, NestedMap
from ..utilities.configurable import Configurable
from ..utilities.errors import IterativeProcessInterrupted, NonSquareSystem

_VEC, _STATE = NumericHandler.VECTORS, NumericHandler.STATE_VEC
_RES, _BOUND = NumericHandler.RES_VEC, NumericHandler.BOUND_VEC


@dataclass
class SimulationSolverIterationReport:
    """This data class object is provided for each iteration during a
    :class:`SimulationSolver` run.
    """
    max_err: float
    r"""For each :class:`~simu.core.utilities.residual.Residual`, the 
    quotient of residual value :math:`r_i` and tolerance :math:`t_i` is 
    calculated. ``max_err`` is the maximum absolute value of these quotients:
    
    .. math:: \mathrm{MET} = \max_i \frac{r_i}{t_i}
    """

    max_res_name: str
    """The name of the :class:`~simu.core.utilities.residual.Residual` which
    causes the value of :attr:`max_err`"""

    relax_factor: float
    """The applied relaxation factor to stay within the domain of the process
    model and the thermodynamic models, according to the defined bounds."""

    min_alpha_name: str
    """The name of the bound that is most limiting and therefore causing the
    value of :attr:`relax_factor`"""

    duration: float
    """The accumulative duration of the solving process inclusive the given
    iteration"""

    lmet: float = field(init=False)  # logarithmic max error to tolerance
    r"""The logarithmic (base 10) value of :attr:`max_err` :math:`r_i / t_i`,
    practically defined as
    
    .. math::
        
        \mathrm{LMET} = \log_{10} \left (
            \max_i \frac{r_i}{t_i} + 10^{-8} \right )

    The offset is introduced to not cause ``NaN`` values for the lucky case in
    which all residuals are exactly zero. This can however easily happen for
    linear systems. Anyhow, :math:`\mathrm{LMET} < 1` is already a sufficient
    condition for convergence.    
    """

    def __post_init__(self):
        self.lmet = log10(self.max_err + 1e-8)


@dataclass
class SimulationSolverReport:
    """The data class object returned from a :class:`SimulationSolver` run
    """
    iterations: Sequence[SimulationSolverIterationReport]
    """A :class:'SimulationSolverIterationReport` object for each performed
    iteration"""

    final_state: Sequence[float]
    """The numerical final state of the model.
    
    .. important::
    
        This state is not suitable for handling initial values in a robust way,
        as it can be very sensitive to minor model changes that for instance
        impact the liquid volumes of equations of state.
        
        Instead, use :meth:`simu.NumericHandler.export_state`,
        :meth:`~simu.NumericHandler.import_state` and 
        :meth:`~simu.NumericHandler.retain_state`.
    """

    prop_func: Callable[[Sequence[float]], NestedMap[Quantity]]
    """The function to calculate all properties of the model as function of
    the state."""

    @property
    def properties(self) -> NestedMap[Quantity]:
        """This property returns all properties of the model, evaluated on
        the :attr:`final_state` attribute. For larger models, this causes
        noticeable computational effort. For this reason, this evaluation is
        only done on demand."""
        return self.prop_func(self.final_state)


SimulationSolverCallback = Callable[
    [int,
     SimulationSolverIterationReport,
     Sequence[float],
     Callable[
         [Sequence[float]],
         NestedMap[Quantity]]
     ],
     bool]
"""A function (type) to act as a call-back in the :class:`SimulationSolver`
solving process. The arguments are as follows:

- ``iteration``: The iteration number as integer, incrementing from zero
- ``report``: The (:class:`SimulationSolverIterationReport`) object
- ``state``: The internal state of the model at given iteration as a sequence of
  floats
- ``prop_func``: A function to calculate all model properties for the given
  state.

The last argument is provided instead of the property structure itself, as it
might be expensive to calculate the entire property structure in each iteration.
This way it can be done on demand.

The callback shall return ``True``, if the solver is to continue, or ``False``
otherwise. In the latter case, the solver will raise a 
:class:`~simu.core.utilities.errors.IterativeProcessInterrupted` exception.

Example:

.. code-block::
   :linenos:

    from pprint import pprint
    
    def my_callback(iteration, report, state, prop_func):
        # This can be a lot to print
        all_properties = prop_func(state)
        pprint(all_properties)

"""

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
                 output: TextIOBase|str = "stdout",
                 call_back_iter: SimulationSolverCallback = None,
                 retain_solution: bool = True):
        r"""On construction, the solver object requires a
        :class:`~simu.NumericHandler` object. The solver object can then be
        reused for multiple solver runs, for instance with variable parameter
        values (sensitivity study).

        :param model: The numeric handler of a model, in most cases obtained by
          the expression ``NumericHandler(ModelClass.top())``.
        :param max_iter: The maximum number of iterations (default 30).

          .. note::

            Normally, 30 iterations should be sufficient. In other words, if the
            model is not converged after 30 iterations, chances are quite low
            that it still will converge at all. The advice would be to try to
            improve the starting values and to investigate whether the model is
            properly posed.

        :param gamma: As described above, :math:`\gamma` (default 0.9) is the
          fraction of the step-length applied by the solver before hitting the
          domain boundary. Normally, changing the value is not required.
          Generally, a lower value makes the model more robust against
          non-linear domain boundaries (and thus linearisation errors causing
          the state to exit the domain. A higher value yields slightly faster
          convergence, if the solution is in comparison with the initial values
          very close to the domain boundary.
        :param wall: Either if there is no solution within the domain of the
          model (for instance: The material balance forces some of the species
          flows in a stream to be negative), or if the solver for other reasons
          is forced to try to leave the model domain, the state will move closer
          and closer to the domain boundary and not revert. At some point,
          :math:`\gamma` becomes ridiculously small, and we need to give up.
          This threshold value is defined by ``wall`` (default ``1e-20``).
        :param output: The io-stream to direct the solver output to, or a
          descriptive string (case-insensitive):

          - ``"stdout"``: The output will be written to standard out (default)
          - ``"none"``: No output will be printed.

          .. note::

            Instead of printing, one might either analyse the returned
            :class:`~simu.core.solver.simulation.SimulationSolverReport` project
            after the run, or utilise the ``call_back_iter`` callback and
            process the iteration progress from there.

        :param call_back_iter: A callback function (default ``None``),
          see :data:`~simu.core.solver.simulation.SimulationSolverCallback`,
          to intercept the solving process. The returned boolean variable
          determines whether the solver iteration is continued or not.

        :param retain_solution: Whether the solver shall, on success, retain the
          obtained state in the model, such that it can be exported via
          :meth:`~simu.NumericHandler.export_state` and be reused as the initial
          values for the next solving process.
        """
        super().__init__(exclude=["model"])
        self._model = model

        # store arguments (parameters) so the user can change them
        args = deepcopy(model.arguments)
        # store size of state
        self.__state_size = args[_VEC][_STATE].magnitude.size()[0]
        res_size = len(model.vector_res_names(NumericHandler.RES_VEC))

        if self.__state_size != res_size:
            raise NonSquareSystem(self.__state_size, res_size)

        # user shall not think that putting a state here has any effect
        del args[_VEC][_STATE]
        self.__model_parameters : NestedMutMap[Quantity] = args

    def solve(self, **kwargs: Any) -> SimulationSolverReport:
        """
        This method triggers iterative the solving process. This takes less than
        0.1 seconds for small models, and increases with model complexity.
        For large models, the time per iteration is due to the solving of linear
        systems cubic in system size, though the model structure might render
        this a conservative estimate.

        With `pypardiso`_ installed, the solving of the linear systems is
        performed on all available CPU cores.

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
        residual_names = model.vector_res_names(_RES)
        bound_names = model.vector_res_names(_BOUND)
        reports = []

        output = self.__find_output()

        table = ProgressTableOutput({
            "lmet": ("LMET", "{:5.1f}"),
            "relax_factor": ("Alpha", "{:7.2g}"),
            "duration": ("Time", "{:6.1g}"),
            "min_alpha_name": ("Limit on bound", "{:>40s}"),
            "max_res_name": ("Max residual", "{:>40s}")
        }, row_dig=5, row_head="Iter", stream=output)

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
                    msg = f"Relaxation factor is below {opt["wall"]}, " \
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

        # retain state if desired
        if opt["retain_solution"]:
            model.retain_state(x, self.model_parameters["thermo_params"])

        return SimulationSolverReport(
            iterations=reports,
            final_state=x,
            prop_func=lambda z: funcs["f_y"]({"x": Quantity(z)})
        )

    def __find_output(self) -> TextIOBase|None:
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
            "gamma": between(0.1, 0.999),
            "call_back_iter": {
                "f": lambda x: x is None or callable(x),
                "msg": "must be callable",
            },
            "output": {
                "f": lambda x: isinstance(x, str)  or isinstance(x, TextIOBase),
                "msg": "must be a stream or a qualified string"
            }
        }
