from dataclasses import dataclass
from typing import Sequence

from numpy.typing import NDArray
from scipy.sparse import csr_array

from simu import Quantity
from simu.core.utilities.types import NestedMutMap

@dataclass
class ThermoFitOuterIterationReport:
    """

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


@dataclass
class ThermoFitReport:
    """

    """
    iterations: Sequence[ThermoFitOuterIterationReport]
    num_data_points: int
    final_parameters: NestedMutMap[Quantity]


@dataclass
class ContributionResult:
    q: Sequence[NDArray]
    dq_dt: Sequence[NDArray]
    num_failed: int


@dataclass
class DataPointResult:
    x: NDArray
    dr_dx: csr_array
