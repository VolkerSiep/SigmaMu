from collections.abc import Sequence
from dataclasses import dataclass
from simu import Quantity
from simu.core.utilities.types import NestedMutMap

@dataclass
class ThermoFitOuterIterationReport:
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


@dataclass
class ThermoFitReport:
    iterations: Sequence[ThermoFitOuterIterationReport]
    parameters: NestedMutMap[Quantity]


