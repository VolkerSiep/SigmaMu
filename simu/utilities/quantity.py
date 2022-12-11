# stdlib modules
from pathlib import Path
from re import split, escape

# external modules
# need to import entire casadi module to distinguish functions of same name
import casadi as cas
from numpy import squeeze
from pint import UnitRegistry, set_application_registry
from pint.errors import DimensionalityError

# instantiate pint unit registry entity
unit_registry = UnitRegistry(autoconvert_offset_to_baseunit=True)
set_application_registry(unit_registry)
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

    def __json__(self):
        """Custom method to export to json for testing"""
        return f"{self:~}"


# write back class, so it's used as result of quantity operations
unit_registry.Quantity = Quantity


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

    def __json__(self):
        """Custom method to export to json for testing"""
        return f"{str(self.magnitude)}{self.units:~}"


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
        m
        >>> print(base_unit("week"))
        s
    """
    if unit == "":
        unit = "dimensionless"
    base = unit_registry.Quantity(unit).to_base_units().units
    return f"{base:~}"


def base_magnitude(quantity: Quantity) -> float:
    """Return the magnitude of the quantity in base units. This works for
    symbolic and numeric quantities, and for scalars and vectors"""
    return quantity.to_base_units().magnitude


# Typing of recursive structures is basically impossible, more so because
# mypy doesn't accept Quantity as a class (with its methods).

QDict = dict[str, Quantity]

_SEPARATOR = "/"  # separator when (un-)flattening dictionaries


def flatten_dictionary(structure, prefix: str = "") -> QDict:
    r"""Convert the given structure into a flat list of key value pairs,
    where the keys are ``SEPARATOR``-separated concatonations of the paths,
    and values are the values of the leafs. Non-string keys are converted
    to strings. Occurances of ``SEPARATOR`` are escaped by ``\``.
    """
    try:
        items = structure.items()  # is this dictionary enough for us?
    except AttributeError:  # doesn't seem so, this is just a value
        return {prefix: structure}  # type: ignore

    result: QDict = {}
    # must sort to create the same sequence every time
    # (dictionary might have content permutated)
    for key, value in sorted(items):
        key = str(key).replace(_SEPARATOR, rf"\{_SEPARATOR}")  # esc. separator
        key = f"{prefix}{_SEPARATOR}{key}" if prefix else key
        result.update(flatten_dictionary(value, key))
    return result


def unflatten_dictionary(flat_structure: QDict):
    """This is the reverse of :func:`flatten_dictionary`, inflating the
    given one-depth dictionary into a nested structure."""
    result = {}  # type: ignore

    def insert(struct, keys, value):
        first = keys.pop(0)
        if keys:
            if first not in struct:
                struct[first] = {}
            insert(struct[first], keys, value)
        else:
            struct[first] = value

    # split by non-escaped separators and unescape escaped separators
    regex = rf'(?<!\\){escape(_SEPARATOR)}'
    for key, value in flat_structure.items():
        keys = [
            k.replace(rf"\{_SEPARATOR}", _SEPARATOR)
            for k in split(regex, key)
        ]
        insert(result, keys, value)
    return result


def extract_units_dictionary(structure):
    """Based on a nested dictionary of Quantities, create a new nested
    dictionaries with only the units of measurement"""
    try:
        items = structure.items()  # is this dictionary enough for us?
    except AttributeError:  # doesn't seem so, this is just a quantity
        return f"{structure.units:~}"

    return {key: extract_units_dictionary(value) for key, value in items}


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
    """

    def __init__(self, args, results, fname: str = "f"):
        args_flat = flatten_dictionary(args)
        results_flat = flatten_dictionary(results)
        arg_names = list(args_flat.keys())
        res_names = list(results_flat.keys())
        arg_sym = [v.magnitude for v in args_flat.values()]
        res_sym = [v.magnitude for v in results_flat.values()]
        self.arg_units = {k: v.units for k, v in args_flat.items()}
        self.res_units = {k: v.units for k, v in results_flat.items()}
        self.func = cas.Function(fname, arg_sym, res_sym, arg_names, res_names)

    def __call__(self, args, squeeze_results=True):
        """Call operator for the function object, as described above."""
        args_flat = {
            key: value.to(self.arg_units[key]).magnitude
            for key, value in flatten_dictionary(args).items()
        }
        result = self.func(**args_flat)  # calling Casadi function
        if squeeze_results:
            result = {k: squeeze(v) for k, v in result.items()}
        result = {k: Quantity(v, self.res_units[k]) for k, v in result.items()}
        return unflatten_dictionary(result)

    @property
    def result_structure(self):
        """Return the result structure as a nested dictionary, only including
        the units of measurements as values of end nodes"""
        units = {k: f"{v:~}" for k, v in self.res_units.items()}
        return unflatten_dictionary(units)
