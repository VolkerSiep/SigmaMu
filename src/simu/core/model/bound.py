from ..utilities import Quantity
from ..utilities.types import Map

BoundProxy = Map["Quantity"]


class BoundHandler(BoundProxy):
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
        value must remain strictly positive.

        If the unit includes an offset, such as ``degC`` or ``barg``, this
        offset is eliminated:

        >>> from simu import SymbolQuantity
        >>> # in a real case, this better be a dependent property
        >>> b = SymbolQuantity("T", "degC")
        >>> handler = BoundHandler()
        >>> handler.add("T", b)
        >>> print(f"{handler["T"]}")
        T delta_degree_Celsius

        If the intention is really to define a bound at a given absolute value,
        this must explicitly be done as follows:

        >>> from simu import Quantity
        >>> handler.add("T2", b - Quantity(31.4159, "degC"))
        >>> print(f"{handler["T2"]}")
        (T-31.4159) delta_degree_Celsius
        """
        if name in self.__bounds:
            raise KeyError(f"Bound {name} already defined")

        # eliminate impact of offset in units like degC and barg
        self.__bounds[name] = bound - Quantity(0.0, bound.units)

    def __getitem__(self, key: str):
        return self.__bounds[key]

    def __len__(self):
        return len(self.__bounds)

    def __iter__(self):
        return iter(self.__bounds)

    def create_proxy(self) -> BoundProxy:
        """Create a proxy object"""
        return self
