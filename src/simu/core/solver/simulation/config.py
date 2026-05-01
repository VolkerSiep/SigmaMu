import sys
from typing import Self
from collections.abc import Callable, Sequence

from pydantic import BaseModel, Field, ConfigDict

from simu import Quantity
from simu.core.utilities.types import NestedMap, OutputIOStream, LinearSolver
from .report import SimulationSolverIterationReport
from ..linear import ScaledLinearSparseSolver


type SimulationSolverCallback = Callable[
    [SimulationSolverIterationReport,
     Sequence[float],
     Callable[
         [Sequence[float]],
         NestedMap[Quantity]]
     ],
    bool]
"""A function (type) to act as a call-back in the 
:class:`~simu.SimulationSolver` solving process. The arguments are as follows:

- ``report``: The data class containing information regarding the current
  iteration
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

    def my_callback(report, state, prop_func):
        # This can be a lot to print
        all_properties = prop_func(state)
        pprint(all_properties)
        return True

"""

class SimulationSolverConfig(BaseModel):
    max_iter: int = Field(default=30, ge=1)
    """The maximum number of iterations (default 30).

    .. note::

      Normally, 30 iterations should be sufficient. In other words, if the
      model is not converged after 30 iterations, chances are quite low
      that it still will converge at all. The advice would be to try to
      improve the starting values and to investigate whether the model is
      properly posed.
    """

    gamma: float = Field(default=0.9, gt=0.0, lt=1.0)
    r""":math:`\gamma` (default 0.9) is the
    fraction of the step-length applied by the solver before hitting the
    domain boundary. Normally, changing the value is not required.
    Generally, a lower value makes the model more robust against
    non-linear domain boundaries (and thus linearisation errors causing
    the state to exit the domain). A higher value yields slightly faster
    convergence, if the solution is in comparison with the initial values
    very close to the domain boundary.
    """

    wall: float = Field(default=1e-20, ge=0.0, lt=0.01)
    r"""Either if there is no solution within the domain of the
    model (for instance: The material balance forces some of the species
    flows in a stream to be negative), or if the solver for other reasons
    is forced to try to leave the model domain, the state will move closer
    and closer to the domain boundary and not revert. At some point,
    :math:`\gamma` becomes ridiculously small, and we need to give up.
    This threshold value is defined by ``wall`` (default ``1e-20``).
    """

    output: OutputIOStream | None = Field(default_factory=lambda: sys.stdout)
    """The stream to direct the solver output to, by default ``sys.stdout``.
    ``None`` suppresses output. 

    The stream can be any object that supports a ``write`` method that consumes
    a string argument.

    .. note::

      Instead of printing, one might either analyse the returned
      :class:`~simu.core.solver.simulation.report.SimulationSolverReport`
      project after the run, or utilise the ``call_back_iter`` callback and
      process the iteration progress from there.
    """

    call_back_iter: SimulationSolverCallback | None = Field(default=None)
    """A callback function (default ``None``),
    see :data:`~simu.core.solver.simulation.config.SimulationSolverCallback`,
    to intercept the solving process. The returned boolean variable
    determines whether the solver iteration is continued or not.
    """

    retain_solution: bool = Field(default=True)
    """Whether the solver shall, on success, retain the
    obtained state in the model, such that it can be exported via
    :meth:`~simu.NumericHandler.export_state` and be reused as the initial
    values for the next solving process.
    """

    linear_solver: LinearSolver = \
        Field(default_factory=ScaledLinearSparseSolver)
    r"""An option to provide any other linear solver for solving the Newton-type
    updates. This can be useful for instance for very large models (say larger
    than 50000 variables, when multi-core and/or iterative solvers become
    superior in terms of performance and robustness.
    """

    model_config = ConfigDict(arbitrary_types_allowed=True)

    def update(self, **options) -> Self:
        data = self.model_dump() | options
        if "output" not in options:
            # reverse undesired irreversible serialization of IO stream
            data["output"] = self.output
        return SimulationSolverConfig.model_validate(data)
