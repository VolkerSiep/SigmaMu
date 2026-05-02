from collections.abc import Sequence, Callable
from dataclasses import dataclass, field
from math import log10
from simu import Quantity
from simu.core.utilities.types import NestedMutMap


type PropertyFunction = Callable[[Sequence[float]], NestedMutMap[Quantity]]
"""A function to calculate model properties from a flat state vector"""

@dataclass
class SimulationSolverIterationReport:
    """This data class object is provided for each iteration during a
    :class:`~simu.SimulationSolver` run.
    """
    iteration: int
    """The number of the iteration, starting with 1"""

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
    """The data class object returned from a :class:`~simu.SimulationSolver` run
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

    prop_func: PropertyFunction
    """The function to calculate all properties of the model as function of
    the state."""

    @property
    def properties(self) -> NestedMutMap[Quantity]:
        """This property returns all properties of the model, evaluated on
        the :attr:`final_state` attribute. For larger models, this causes
        noticeable computational effort. For this reason, this evaluation is
        only done on demand."""
        return self.prop_func(self.final_state)

