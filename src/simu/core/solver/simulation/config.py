from collections.abc import Callable, Sequence
from io import TextIOBase

from pydantic import BaseModel, Field

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
    gamma: float = Field(default=0.9, gt=0.0, lt=1.0)
    wall: float = Field(default=1e-20, ge=0.0, lt=0.01)
    output: TextIOBase | str = Field(default="stdout")
    call_back_iter: SimulationSolverCallback | None = Field(default=None)
    retain_solution: bool= Field(default=True)
