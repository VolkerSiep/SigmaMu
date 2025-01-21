from ..utilities import Quantity
from ..utilities.types import Map


class BoundHandler(Map[Quantity]):
    """The boundary handler keeps process model properties that must be strictly
    positive in order to remain within the model's mathematical domain.

    Thermodynamic (material) properties, such as temperature, pressure and
    molar quantities do not need to be added, as their bounds are to be defined
    in the :class:`simu.ThermoContribution` object(s) that constraint(s) their
    domain.

    .. note::

        These bounds do not represent inequality constraints in an optimization
        problem, but true domain boundaries that need to be respected by any
        solver in order to keep a feasible vector of independent variables.

        It is also generally not a good idea to impose bounds that are based on
        the expected range for a solved model, as a solver might need to step
        through values that are not physically reasonable but yet within the
        model domain.

    Typical examples of bounds defined here are:

      - Temperature difference over a heat exchanger when using the logarithmic
        mean temperature. A temperature cross-over would cause negative
        arguments to the ``log`` function.
      - Process model parameters that might be subject to change, such as
        geometric dimensions of equipment.

    A **bad** example is using these bounds to keep a temperature value between
    a lower and an upper limit.
    """
    def __init__(self):
        self.__bounds = {}

    def add(self, name: str, bound: Quantity):
        """Add a quantity to the bound handler to signal the solver that its
        value must remain strictly positive."""
        if name in self.__bounds:
            raise KeyError(f"Bound {name} already defined")
        self.__bounds[name] = bound

    def __getitem__(self, key: str):
        return self.__bounds[key]

    def __len__(self):
        return len(self.__bounds)

    def __iter__(self):
        return iter(self.__bounds)


