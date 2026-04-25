import sys
from collections.abc import Callable, Sequence
from io import TextIOBase
from typing import Any

from pydantic import BaseModel, Field, field_validator

from simu import Quantity
from simu.core.utilities.types import NestedMap
from .report import SimulationSolverIterationReport


type SimulationSolverCallback = Callable[
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

    output: TextIOBase | str = Field(default="stdout")
    """The io-stream to direct the solver output to, or a
    descriptive string (case-insensitive):

    - ``"stdout"``: The output will be written to standard out (default)
    - ``"none"``: No output will be printed.

    .. note::

      Instead of printing, one might either analyse the returned
      :class:`~simu.core.solver.simulation.SimulationSolverReport` project
      after the run, or utilise the ``call_back_iter`` callback and
      process the iteration progress from there.
    """

    call_back_iter: SimulationSolverCallback | None = Field(default=None)
    """A callback function (default ``None``),
    see :data:`~simu.core.solver.simulation.SimulationSolverCallback`,
    to intercept the solving process. The returned boolean variable
    determines whether the solver iteration is continued or not.
    """

    retain_solution: bool = Field(default=True)
    """Whether the solver shall, on success, retain the
    obtained state in the model, such that it can be exported via
    :meth:`~simu.NumericHandler.export_state` and be reused as the initial
    values for the next solving process.
    """

    @field_validator("output", mode="before")
    @classmethod
    def _validate_output(cls, value: Any) -> TextIOBase | None:
        if isinstance(value, str):
            opts = {"stdout": sys.stdout, "none": None}
            try:
                return opts[value.lower()]
            except KeyError:
                raise ValueError(f"Invalid stream name '{value}'")
        return value
