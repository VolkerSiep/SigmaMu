from dataclasses import dataclass
from typing import Sequence

from numpy.typing import NDArray
from scipy.sparse import csr_array

from simu import Quantity
from simu.core.utilities.types import NestedMutMap

@dataclass
class ThermoFitOuterIterationReport:
    """
    Captures the state and metrics of a single outer iteration within the
    ThermoFit solver process.

    This report tracks convergence metrics, relaxation parameters, and the
    parameter state at a specific iteration step.
    """

    iteration: int
    """The number of the iteration, starting with 1"""

    q_norm: float
    r"""The norm of the penalty contributions"""

    stationarity: float
    r"""The closeness to the stationary condition"""

    relax_factor: float
    """The applied relaxation factor to stay within the domain of the process
    model and the thermodynamic models, according to the defined bounds."""

    duration: float
    """The accumulative duration of the solving process inclusive the given
    iteration"""

    tau: NestedMutMap[Quantity]
    """The current set of parameters"""

    num_failed: int
    """The number of failed data point evaluations during this iteration"""


@dataclass
class ThermoFitReport:
    """
    Aggregates the results and history of a complete ThermoFit solver execution.

    This report contains the full history of iterations, the total number of
    data points processed, and the final optimized parameter set.
    """
    iterations: Sequence[ThermoFitOuterIterationReport]
    """A sequence of reports for each outer iteration performed."""

    num_data_points: int
    """The total number of data points used in the fitting process."""

    final_parameters: NestedMutMap[Quantity]
    """The final set of parameters resulting from the optimization."""


@dataclass
class ContributionResult:
    q: Sequence[NDArray]
    dq_dt: Sequence[NDArray]
    num_failed: int


@dataclass
class DataPointResult:
    x: NDArray
    dr_dx: csr_array
