from pathlib import Path

# need to import entire numpy and casadi modules to distinguish functions
# of same name
import casadi as cas
from pint import UnitRegistry
from pint.errors import DimensionalityError
from typing import Mapping, Union

# instantiate pint unit registry entity
unit_registry = UnitRegistry(autoconvert_offset_to_baseunit=True)
unit_file = Path(__file__).resolve().parent / "uom_definitions.txt"
unit_registry.load_definitions(unit_file)
del unit_file

_Q = unit_registry.Quantity


class Quantity(_Q):  # type: ignore
    """Proper quantity base-class for sub-classing"""

    def __new__(cls, *args, **kwargs):
        obj = super().__new__(_Q, *args, **kwargs)
        obj.__class__ = cls
        return obj


class SymbolQuantity(Quantity):
    """A quantity class specialised to host casadi symbols (SX)."""

    def __new__(cls, *args, **kwargs):
        """Really generate an object of type ``Quantity``. This is just a
        hacky way to specialise the constructor, due to the way the base-class
        is implemented."""

        def attributes(name: str, units: str, sub_keys: list[str] = None):
            if sub_keys is None:
                magnitude = cas.SX.sym(name)
            elif isinstance(sub_keys, int):
                magnitude = cas.SX.sym(name, sub_keys)
            else:
                magnitude = [cas.SX.sym(f"{name}.{s}") for s in sub_keys]
                magnitude = cas.vertcat(*magnitude)
            return magnitude, units

        return super().__new__(Quantity, *attributes(*args, **kwargs))


# redefine jacobian to propanatural logarithmusgate pint units
def jacobian(dependent: Quantity, independent: SymbolQuantity) -> Quantity:
    """Calculate the casadi Jacobian and reattach the units of measurements."""
    return Quantity(cas.jacobian(dependent.magnitude, independent.magnitude),
                    dependent.units / independent.units)


def sum1(quantity: Quantity) -> Quantity:
    """Sum a symbol vector quantity same was a ``casadi.sum1``, considering
    units of measurements"""
    return Quantity(cas.sum1(quantity.magnitude), quantity.units)


def log(quantity: Quantity) -> Quantity:
    """Determine natural logarithmus of symbolic quantity, considering units of
    measurements"""
    quantity = quantity.to_base_units()
    if not quantity.dimensionless:
        raise DimensionalityError(quantity.units, "dimensionless")
    return Quantity(cas.log(quantity.magnitude))


def exp(quantity: Quantity) -> Quantity:
    """Determine exponent of symbolic quantity, considering units of
    measurements"""
    quantity = quantity.to_base_units()
    if not quantity.dimensionless:
        raise DimensionalityError(quantity.units, "dimensionless")
    return Quantity(cas.exp(quantity.magnitude))


def sqrt(quantity: Quantity) -> Quantity:
    """Determine square root of symbolic quantity, considering units of
    measurement"""
    return Quantity(cas.sqrt(quantity.magnitude), quantity.units**0.5)


def qpow(base: Quantity, exponent: Quantity) -> Quantity:
    """Determine power of ``base`` to ``exponent``, considering units of
    measurements. Both arguments must be dimensionless. For the special case
    that ``exponent`` is a constant scalar (not a symbol), use the ``**``
    operator
    """
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
    depending on the condition and considering the units of measurements"""

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

    .. code-block::

        >>> print(base_unit("light_year"))
        meter
        >>> print(base_unit("week"))
        second
    """
    if unit == "":
        unit = "dimensionless"
    base = unit_registry.Quantity(unit).to_base_units().units
    return f"{base:~}"


def base_magnitude(quantity: Quantity) -> float:
    """Return the magnitude of the quantity in base units. This works for
    symbolic and numeric quantities, and for scalars and vectors"""
    return quantity.to_base_units().magnitude


NestedQDict_1 = Mapping[str, Quantity]
NestedQDict = NestedQDict_1  # Mapping[str, Union[NestedQDict_1, Quantity]]


class QFunction:
    """Wrapper around ``casadi.Function`` to consider units of measurements.
    This derived function object is defined by a dictionary of arguments and
    results, both of which are ``Quantity`` instances with casadi symbols as
    magnitudes. These symbols are independent symbols for ``args`` and derived
    symbols for ``results``.

    The function object is then called as well with a dictionary of
    ``Quantity`` values. The units of measurement must be consistent, and a
    conversion is done to the initially defined units for the individual
    arguments. The result is given back as a dictionary of ``Quantity`` values
    in the same units as initially defined.

    .. todo::
        extend to allow nested dictionaries, that will only be flattened
        internally, so the user doesn't need to care about it at all.
    """

    def __init__(self,
                 args: dict[str, SymbolQuantity],
                 results: dict[str, Quantity],
                 fname: str = "f"):
        arg_names = list(args.keys())
        res_names = list(results.keys())
        arg_sym = [v.m for v in args.values()]
        res_sym = [v.m for v in results.values()]
        self.arg_units = {k: v.u for k, v in args.items()}
        self.res_units = {k: v.u for k, v in results.items()}
        self.func = cas.Function(fname, arg_sym, res_sym, arg_names, res_names)

    def __call__(self, args: NestedQDict) -> NestedQDict:
        """Call operator for the function object, as described above."""
        args = {k: v.to(self.arg_units[k]).m for k, v in args.items()}
        result = self.func(**args)  # calling Casadi function
        return {k: Quantity(v, self.res_units[k]) for k, v in result.items()}
