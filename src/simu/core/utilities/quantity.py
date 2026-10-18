"""This module defines classes and functionality around ``pint`` quantities.
These quantities can be symbolic (hosting ``casadi.SX`` as magnitudes), or
numeric."""

# stdlib modules
from collections.abc import Mapping

# external modules
import casadi as cas  # import as module to distinguish name-equal functions
from pint import UnitRegistry, set_application_registry, Quantity as QtyType
from pint.util import UnitsContainer

# internal modules
from simu.core.data import DATA_DIR
from .types import NestedMap
from .errors import DimensionalityError


def _create_registry() -> UnitRegistry:
    """Instantiate the application's unit registry"""
    reg = UnitRegistry(autoconvert_offset_to_baseunit=True)
    set_application_registry(reg)
    unit_file = DATA_DIR / "uom_definitions.txt"
    reg.load_definitions(unit_file)
    return reg


unit_registry = _create_registry()
_Q = unit_registry.Quantity

del _create_registry

#  variables to simplify units of measurements
__DERIVED_UNITS = \
    "joule watt newton pascal coulomb volt farad ohm " \
    "siemens weber tesla henry lumen lux".split()
__simplify_quantity_cache = {}


class Quantity(_Q):
    """Proper quantity base-class for sub-classing.

    Being a subclass of ``pint.Quantity``, this class only really adds the
    ``__json__`` method to return its JSON representation.

    The constructor is used as for ``pint.Quantity``.
    """

    def __new__(cls, *args, **kwargs):
        obj = super().__new__(_Q, *args, **kwargs)
        obj.__class__ = cls
        return obj

    def __json__(self) -> str:
        """Custom method to export to JSON for testing and serialization"""
        try:
            return f"{self:~.16g}"
        except TypeError:  # magnitude cannot be formatted, probably symbol
            return f"{self:~}"

    # operators needed for proper type checking in PyCharm (so long at least)
    def __add__(self, other) -> Quantity:
        return super().__add__(other)

    def __radd__(self, other) -> Quantity:
        return super().__radd__(other)

    def __sub__(self, other) -> Quantity:
        return super().__sub__(other)

    def __rsub__(self, other) -> Quantity:
        return super().__rsub__(other)

    def __mul__(self, other) -> Quantity:
        return super().__mul__(other)

    def __rmul__(self, other) -> Quantity:
        return super().__rmul__(other)

    def __matmul__(self, other):
        try:
            factor = other.magnitude
            unit = self.units * other.units
        except AttributeError:
            unit = self.units
            factor = other
        return Quantity(self.m @ factor, unit)

    def __rmatmul__(self, other):
        try:
            factor = other.magnitude
            unit = self.units * other.units
        except AttributeError:
            factor = other
            unit = self.units
        return Quantity(factor @ self.m, unit)

    def __truediv__(self, other) -> Quantity:
        return super().__truediv__(other)

    def __rtruediv__(self, other) -> Quantity:
        return super().__rtruediv__(other)

    def __pow__(self, other) -> Quantity:
        return super().__pow__(other)

    def __rpow__(self, other) -> Quantity:
        return super().__rpow__(other)


# write back class, so it's used as result of quantity operations
unit_registry.Quantity = Quantity

class SymbolQuantity(Quantity):
    """A quantity class specialized to host casadi symbols (SX) of in
    particular, but not necessarily, independent variables."""

    def __new__(cls, *args, **kwargs):
        """Really generate an object of type ``Quantity``. This is just a
        hacky way to specialize the constructor, due to the way the base-class
        is implemented. The arguments are:

        - ``name`` (str): The name of the ``casadi.SX`` symbol
        - ``units`` (str): The unit of measurement conform to ``pint`` units
        - ``sub_keys`` (Iterator[str]): If not ``None``, defines a vectorial
          quantity with given sub-keys.

            >>> s = SymbolQuantity("speed", "m/s")
            >>> print(f"{s:~}")
            speed m / s

            >>> v = SymbolQuantity("vel", "m/s", "xyz")
            >>> print(f"{v:~}")
            [vel.x, vel.y, vel.z] m / s
        """

        def attributes(name: str, units: str,
                       sub_keys: list[str] | int | None = None):
            """Turn the constructor arguments into arguments for the
            baseclass"""
            if sub_keys is None:
                magnitude = cas.SX.sym(name)
            elif isinstance(sub_keys, int):
                magnitude = cas.SX.sym(name, sub_keys)
            else:
                magnitude = [cas.SX.sym(f"{name}.{s}") for s in sub_keys]
                magnitude = cas.vertcat(*magnitude)
            return magnitude, units

        return super().__new__(Quantity, *attributes(*args, **kwargs))

    def __json__(self) -> str:
        """Custom method to export to JSON for testing"""
        return f"{str(self.magnitude)}{self.units:~}"


# redefine Jacobian to propagate pint units
def jacobian(dependent: Quantity, independent: Quantity) -> Quantity:
    """Calculate the casadi Jacobian and reattach the units of measurements.

    >>> x = SymbolQuantity("x", "m")
    >>> y = (x * x) / 2
    >>> dy_dx = jacobian(y, x)
    >>> print(f"{dy_dx:~}")
    x m
    """
    magnitude = cas.jacobian(dependent.magnitude, independent.magnitude)
    return Quantity(cas.simplify(magnitude),
                    dependent.units / independent.units)


def qsum(quantity: Quantity) -> Quantity:
    """Sum a symbol vector quantity same was a ``casadi.sum1``, considering
    units of measurements. This function only applies to quantity objects with
    ``casadi.SX`` objects as magnitudes.

    .. note::

        This function sums the elements of the single, but vectorial argument,
        and is in that different from the builtin ``sum`` function that sums
        over an iterator of objects.
    """
    return Quantity(cas.sum1(quantity.magnitude), quantity.units)


def qsqrt(quantity: Quantity | float) -> Quantity:
    """Determine square root of symbolic quantity, considering units of
    measurement"""
    if not isinstance(quantity, Quantity):
        quantity = Quantity(quantity)
    return Quantity(cas.sqrt(quantity.magnitude), quantity.units**0.5)


def qpow(base: Quantity | float, exponent: Quantity | float) -> Quantity:
    """Determine power of ``base`` to ``exponent``, considering units of
    measurements. Both arguments must be dimensionless. For the special case
    that ``exponent`` is a constant scalar (not a symbol), use the ``**``
    operator
    """
    if not isinstance(base, Quantity):
        base = Quantity(base)
    if not isinstance(exponent, Quantity):
        exponent = Quantity(exponent)

    base = base.to_base_units()  # else e.g. cm / m is dimensionless
    exponent = exponent.to_base_units()
    if not base.dimensionless:
        raise DimensionalityError(base.units, "dimensionless")
    if not exponent.dimensionless:
        raise DimensionalityError(exponent.units, "dimensionless")
    return Quantity(base.magnitude**exponent.magnitude)


def conditional(condition: cas.SX, negative: Quantity,
                positive: Quantity) -> Quantity:
    """Element-wise branching into the positive and negative branch
    depending on the condition and considering the units of measurements.

        >>> x = SymbolQuantity("x", "m")
        >>> y = conditional(x > 0, -x, x)  # abs function
        >>> print(y)
        @1=0, @2=((@1<x)==@1), ((@2?(-x):0)+((!@2)?x:0)) meter

    As non-recommended as it is to use this function excessively, the resulting
    ``casadi`` expression is indeed overcomplicated and could simplify to ::

        (x>0)?(x):(-x) meter

    .. note::

        You cannot just code ``x if x > 0 else -x``, as this would not allow
        ``casadi`` branching dynamically dependent on the value of ``x``
        later-on.
    """
    units = positive.units
    negative = negative.to(units)  # convert to same unit
    tup_branches = [condition, negative.magnitude, positive.magnitude]
    tup = zip(*map(cas.vertsplit, tup_branches))

    magnitude = cas.vertcat(*[cas.conditional(c, [n], p) for c, n, p in tup])
    return Quantity(magnitude, units)


def qvertcat(*quantities: Quantity) -> Quantity:
    """Concatenate a bunch of scalar symbolic quantities with compatible units
    to a vector quantity"""
    units = quantities[0].units
    magnitudes = [q.to(units).magnitude for q in quantities]
    return Quantity(cas.vertcat(*magnitudes), units)


def base_unit(unit: str) -> str:
    """Create the base unit of given unit.

    >>> print(base_unit("light_year"))
    m
    >>> print(base_unit("week"))
    s
    """
    if unit == "":
        unit = "dimensionless"
    base = Quantity(unit).to_base_units().units
    return f"{base:~}"


def simplify_quantity(quantity: Quantity) -> Quantity:
    """Try to convert the unit into a more compact notation, allowing
    derived units, such as Watt, Joule and Pascal to be used.

    Examples:

    >>> q = Quantity(1.0, "kg * m**2 / mol / s**2")
    >>> print(f"{simplify_quantity(q):~}")
    1.0 J / mol

    >>> q = Quantity(1.0, "kg / K / s**3")
    >>> print(f"{simplify_quantity(q):~}")
    1.0 W / K / m ** 2

    >>> q = Quantity(1.0, "kg / m**2 / s**2")
    >>> print(f"{simplify_quantity(q):~}")
    1.0 Pa / m

    >>> q = Quantity(1.0, "kg * m**2 / s**3 / A**2")
    >>> print(f"{simplify_quantity(q):~}")
    1.0 Ω

    >>> q = Quantity(1.0, "m**3 / s")
    >>> print(f"{simplify_quantity(q):~}")
    1.0 m ** 3 / s
    """
    def weight_factors(dimensionality: UnitsContainer) -> float:
        """Count the factors, but weight an exponent unequal unity less hard
        than a completely new factor
         """
        return sum(1 + abs(e) for e in dimensionality.values())

    quantity = quantity.to_base_units()
    original_dim = tuple(quantity.dimensionality.items())
    try:
        new_unit = __simplify_quantity_cache[original_dim]
        return quantity.to(new_unit)
    except KeyError:
        pass
    count = weight_factors(quantity.dimensionality)
    for unit in __DERIVED_UNITS:
        trial = Quantity(1, unit)
        candidate = (quantity / trial).to_base_units()
        new_count = weight_factors(candidate.dimensionality)
        if new_count < count:
                quantity = candidate * trial
                count = new_count
    __simplify_quantity_cache[original_dim] = quantity.units
    return quantity


def base_magnitude[T](quantity: QtyType[T]) -> T:
    """Return the magnitude of the quantity in base units. This works for
    symbolic and numeric quantities, and for scalars and vectors

    >>> base_magnitude(Quantity("1 km"))
    1000.0
    >>> base_magnitude(Quantity("20 degC"))
    293.15
    >>> int(base_magnitude(Quantity("1 barg")))
    201325
    >>> base_magnitude(Quantity("speed_of_light"))
    299792458.0

    The base units are likely the SI unit system, but code shall not rely on
    this fact - only that it is a consistent (and offset-free) unit system.

    .. note::

        "**Use SI or I will BTU you with my feet!**"

        I saw this sentence once on the T-shirt of a nerd.
        It turns out it was a mirror.
    """
    return quantity.to_base_units().magnitude


def extract_units_dictionary(structure: NestedMap[Quantity] | Quantity) \
        -> NestedMap[str] | str:
    """Based on a nested dictionary of Quantities, create a new nested
    dictionaries with only the (base) units of measurement

    >>> d = {"a": {"b": Quantity("1 m"), "c": Quantity("1 min")}}
    >>> extract_units_dictionary(d)
    {'a': {'b': 'm', 'c': 's'}}
    """
    if isinstance(structure, Mapping):
        return {k: extract_units_dictionary(v) for k, v in structure.items()}
    assert isinstance(structure, Quantity)
    return f"{simplify_quantity(structure).units:~}"
