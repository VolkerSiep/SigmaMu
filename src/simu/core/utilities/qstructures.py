"""This module contains data structures that build on the quantity datatype"""

# stdlibs
from typing import Union, TypeVar, Self
from collections.abc import Callable, Iterable, Mapping

# external libs
import casadi as cas
from pint.errors import DimensionalityError

# internal libs
from . import base_unit, Quantity, SymbolQuantity, qvertcat, qpow, qsqrt
from .types import NestedMap


class ParameterDictionary(dict):
    """This class is a nested dictionary of SymbolQuantities to represent
    parameters with functionality to be populated using the ``register_*``
    methods.
    """

    class SparseMatrix(dict):
        """This helper class represents a nested dictionary that contains
        two levels of keys and values representing a quantity."""

        def pair_items(self):
            """Return an iterator yielding a scalar quantity with the key pair
            for each element in the sub-structure. The elements have the
            shape ``(key_1, key_2, quantity)``."""
            for key_1, second in self.items():
                for key_2, quantity in second.items():
                    yield key_1, key_2, quantity

    def register_scalar(self, key: str, unit: str):
        """Create a scalar quantity and add the structure to the dictionary.
        The given unit is converted to base units before being applied. Calling
        the method returns the created quantity

            >>> pdict = ParameterDictionary()
            >>> print(pdict.register_scalar("speed", "cm/h"))
            speed meter / second

        In this output, ``speed`` is the name of the ``casadi.SX`` node
        representing the magnitude of returned Quantity. The dictionary then
        contains the following entry:

            >>> print(pdict)
            {'speed': <Quantity(speed, 'meter / second')>}
        """
        unit = base_unit(unit)
        quantity = SymbolQuantity(key, unit)
        self[key] = quantity
        return quantity

    def register_vector(self, key: str, sub_keys: Iterable[str],
                        unit: str) -> Quantity:
        """Create a quantity vector with symbols and add the structure to
        the dictionary. The given unit is converted to base units before being
        applied. Calling the method returns the created quantity

            >>> pdict = ParameterDictionary()
            >>> print(pdict.register_vector("velocity", "xyz", "knot"))
            [velocity.x, velocity.y, velocity.z] meter / second

        The dictionary then contains the following entries:

            >>> from pprint import pprint
            >>> pprint(pdict)
            {'velocity': {'x': <Quantity(velocity.x, 'meter / second')>,
                          'y': <Quantity(velocity.y, 'meter / second')>,
                          'z': <Quantity(velocity.z, 'meter / second')>}}
        """
        unit = base_unit(unit)
        self[key] = {s: SymbolQuantity(f"{key}.{s}", unit) for s in sub_keys}
        return qvertcat(*self[key].values())

    def register_sparse_matrix(self, key: str, pairs: Iterable[tuple[str, str]],
                               unit: str) -> NestedMap[Quantity]:
        """Create a sparse matrix quantity and add the structure to the
        dictionary. The given unit is converted to base units before being
        applied.

            >>> pdict = ParameterDictionary()
            >>> binaries = [("H2O", "CO2"), ("H2O", "CH4")]
            >>> from pprint import pprint
            >>> pprint(pdict.register_sparse_matrix("K_ij", binaries, "K"))
            {'H2O': {'CH4': <Quantity(K_ij.H2O.CH4, 'kelvin')>,
                     'CO2': <Quantity(K_ij.H2O.CO2, 'kelvin')>}}

        After above call, the dictionary contains the following entries:

            >>> from pprint import pprint
            >>> pprint(pdict)
            {'K_ij': {'H2O': {'CH4': <Quantity(K_ij.H2O.CH4, 'kelvin')>,
                              'CO2': <Quantity(K_ij.H2O.CO2, 'kelvin')>}}}

        """
        unit = base_unit(unit)
        res = ParameterDictionary.SparseMatrix({f: {} for f, _ in pairs})
        for first, second in pairs:
            quantity = SymbolQuantity(f"{key}.{first}.{second}", unit)
            res[first][second] = quantity
        self[key] = res
        return res

    def get_quantity(self, *keys):
        """Extract a quantity from the given sequence of key. Being a nested
        dictionary, each key from the argument list is used to navigate into
        the structure. The value of the most inner addressed key is returned.
        For normal usage, this should be of type ``Quantity``."""
        entry = self
        for key in keys:
            entry = entry[key]
        return entry

    def get_vector_quantity(self, *keys):
        """Extract a vector quantity from the given sequence of keys. The
        method extracts the values of the structure below the sequence of
        argument keys, and concatenates them as a single vector property.
        """
        entry = self
        for key in keys:
            entry = entry[key]
        return qvertcat(*entry.values())


_OType = Union[float, Quantity, Mapping[str, Quantity]]
_SType = Union[float, cas.SX]


class QuantityDict(dict[str, Quantity]):
    """Many properties on process modelling level are vectorial. This includes
    any species-specific properties, such as for instance mole fractions,
    chemical potentials or partial enthalpy. By keeping such data in instances
    of this class, they can always be accessed as a dictionary, using the
    bracket-operator (``__get_item__``).

    Additionally, this class supports most arithmetic operations, such as
    ``+, -, *, /, **`` - all of them interpreted element-wise. As two instances
    can have deviating keys (mostly species), there are some rules:

    - missing elements are assumed as zero, and structural zeros are
      omitted, i.e.

      >>> a = QuantityDict({
      ...         "A": Quantity("1 m"),
      ...         "B": Quantity("50 cm")})
      >>> b = QuantityDict({
      ...         "B": Quantity("1 m"),
      ...         "C": Quantity("50 cm")})
      >>> y = a + b
      >>> for key, value in y.items(): print(f"{key}: {value:~}")
      A: 1 m
      B: 150.0 cm
      C: 50 cm

      >>> y = a * b
      >>> for key, value in y.items(): print(f"{key}: {value:~}")
      B: 50 cm * m

    - A missing denominator element in division directly raises
      ``ZeroDivisionError``

      >>> y = a / b
      Traceback (most recent call last):
      ...
      ZeroDivisionError: Missing denominator element in QuantityDict division

    - Operations can be mixed with scalar Quantities

      >>> y = a["A"] * b
      >>> for key, value in y.items(): print(f"{key}: {value:~}")
      B: 1 m ** 2
      C: 50 cm * m

    - floats as second operands in binary operators act as dimensionless
      quantities

      >>> y = 3 * a
      >>> for key, value in y.items(): print(f"{key}: {value:~}")
      A: 3 m
      B: 150 cm

      >>> y = 3 + a
      Traceback (most recent call last):
      ...
      pint.errors.DimensionalityError: ...
    """
    @classmethod
    def from_vector_quantity(cls, quantity: Quantity,
                             keys: list[str]) -> Self:
        """As the magnitude of a ``Quantity`` can be a container itself,
        this convenience factory method combines such a vector quantity with
        a set of given keys into a ``QuantityDict`` object

        >>> raw = Quantity([1, 2, 3], "m")
        >>> dic = QuantityDict.from_vector_quantity(raw, ["A", "B", "C"])
        >>> for name, value in dic.items():
        ...     print(f"{name}: {value:~}")
        A: 1 m
        B: 2 m
        C: 3 m
        """
        try:
            l_magnitude = len(quantity.magnitude)
        except TypeError:
            l_magnitude = quantity.magnitude.size1()
        if l_magnitude != len(keys):
            raise ValueError("Dimension mismatch in resolving vector quantity: "
                             f"{l_magnitude} != {len(keys)}")
        return cls({key: quantity[k] for k, key in enumerate(keys)})

    def __add__(self, other: _OType) -> Self:
        try:
            items = other.items()
        except AttributeError:
            return QuantityDict({k: v + other for k, v in self.items()})

        result = self.copy()
        for key, value in items:
            result[key] = (result[key] + value) if key in self else value
        return QuantityDict(result)

    def __radd__(self, other: _OType) -> Self:
        return self + other

    def __mul__(self, other: _OType) -> Self:
        try:
            other.items()
        except AttributeError:
            return QuantityDict({k: v * other for k, v in self.items()})

        result = {k: v * other[k] for k, v in self.items() if k in other}
        return QuantityDict(result)

    def __rmul__(self, other: _OType) -> Self:
        return self * other

    def __pos__(self) -> Self:
        return self

    def __neg__(self) -> Self:
        return QuantityDict({k: -v for k, v in self.items()})

    def __sub__(self, other: _OType) -> Self:
        return self + (-other)

    def __rsub__(self, other: _OType) -> Self:
        return (-self) + other

    def __truediv__(self, other: _OType) -> Self:
        try:
            other.items()
        except AttributeError:
            return QuantityDict({k: v / other for k, v in self.items()})

        try:
            result = {k: v / other[k] for k, v in self.items()}
        except KeyError:
            msg = "Missing denominator element in QuantityDict division"
            raise ZeroDivisionError(msg) from None
        return QuantityDict(result)

    def __rtruediv__(self, other: _OType) -> Self:
        try:
            items = other.items()
        except AttributeError:
            return QuantityDict({k: other / v for k, v in self.items()})

        try:
            result = {k: v / self[k] for k, v in items}
        except KeyError:
            msg = "Missing denominator element in QuantityDict division"
            raise ZeroDivisionError(msg) from None
        return QuantityDict(result)

    def __matmul__(self, other: Self) -> Quantity:
        """Can assume other to be QuantityDict, otherwise use
        *-operator instead."""
        return (self * other).sum()

    def __rmatmul__(self, other) -> Quantity:
        """For 2 vectors, scalar product is commutative"""
        return other @ self

    def __pow__(self, other: _OType) -> Self:
        try:
            items = other.items()
        except AttributeError:
            result = {k: qpow(v, other) for k, v in self.items()}
        else:
            result = {k: qpow(v, other.get(k, 0)) for k, v in self.items()}
            for k, v in items:
                if k not in self:
                    result[k] = qpow(0.0, v)
        return QuantityDict(result)

    def __rpow__(self, other: _OType) -> Self:
        try:
            items = other.items()
        except AttributeError:
            result = {k: qpow(other, v) for k, v in self.items()}
        else:
            result = {k: qpow(other.get(k, 0) , v) for k, v in self.items()}
            for k, v in items:
                if k not in self:
                    result[k] = qpow(v, 0.0)
        return QuantityDict(result)

    def sum(self) -> Quantity:
        """Sum all elements of the QuantityDict object. Naturally, all
        elements must have equal physical dimensions

        >>> a = QuantityDict({
        ...         "B": Quantity("1 m"),
        ...         "C": Quantity("50 cm")})
        >>> print(f"{a.sum():~}")
        1.5 m
        """
        return sum(self.values())  # type: ignore[type-error]


_V = TypeVar("_V", float, Quantity, Mapping[str, Quantity])
UNFUNC_TYPE = Callable[[_V], _V]


def unary_func(quantity: _V, func: UNFUNC_TYPE) -> _V:
    """Call a unary function that requires a dimensionless argument on the
    argument, accepting that argument to be a float, a scalar quantity,
    a symbolic quantity, or a QuantityDict object (well, any dictionary
    of the above, actually)."""
    def scalar_func(value):
        """Call the unary function on the value's magnitude"""
        try:
            value = value.to_base_units()
            magnitude = value.magnitude
        except AttributeError:
            magnitude = value
        else:
            if not value.dimensionless:
                raise DimensionalityError(value.units, "dimensionless")
        return Quantity(func(magnitude))

    try:
        items = quantity.items()
    except AttributeError:
        return scalar_func(quantity)
    return QuantityDict({k: unary_func(v, func) for k, v in items})


def sqrt(quantity: _V) -> _V:
    """The square root function is a special case of a unary function in that
    the argument is not required to be dimensionless.

    >>> a = QuantityDict({
    ...         "B": Quantity("1 m**2"),
    ...         "C": Quantity("2500 cm**2")})
    >>> print(sqrt(a))
    {'B': <Quantity(1.0, 'meter')>, 'C': <Quantity(50.0, 'centimeter')>}
    """
    try:
        items = quantity.items()
    except AttributeError:
        return qsqrt(quantity)
    return QuantityDict({k: qsqrt(v) for k, v in items})


def log(quantity: _V) -> _V:
    """Determine natural logarithms, considering units of
    measurements. The main intent is to use this version for symbolic
    quantities and QuantityDict objects, but it also works on floats.

    >>> x = Quantity(10.0, "cm/m")
    >>> log(x)
    <Quantity(-2.30258509, 'dimensionless')>

    >>> a = {"A": SymbolQuantity("A", "dimless"),
    ...      "B": SymbolQuantity("B", "dimless")}
    >>> y = log(a)
    >>> for key, value in y.items(): print(f"{key}: {value:~}")
    A: log(A)
    B: log(B)

    >>> log(10)
    <Quantity(2.30258509, 'dimensionless')>

    The other unary functions are defined in the same manner.
    """
    return unary_func(quantity, cas.log)


log10: UNFUNC_TYPE = lambda x: unary_func(x, cas.log10)
exp: UNFUNC_TYPE = lambda x: unary_func(x, cas.exp)
sin: UNFUNC_TYPE = lambda x: unary_func(x, cas.sin)
cos: UNFUNC_TYPE = lambda x: unary_func(x, cas.cos)
tan: UNFUNC_TYPE = lambda x: unary_func(x, cas.tan)
arcsin: UNFUNC_TYPE = lambda x: unary_func(x, cas.arcsin)
arccos: UNFUNC_TYPE = lambda x: unary_func(x, cas.arccos)
arctan: UNFUNC_TYPE = lambda x: unary_func(x, cas.arctan)
sinh: UNFUNC_TYPE = lambda x: unary_func(x, cas.sinh)
cosh: UNFUNC_TYPE = lambda x: unary_func(x, cas.cosh)
tanh: UNFUNC_TYPE = lambda x: unary_func(x, cas.tanh)
arcsinh: UNFUNC_TYPE = lambda x: unary_func(x, cas.arcsinh)
arccosh: UNFUNC_TYPE = lambda x: unary_func(x, cas.arccosh)
arctanh: UNFUNC_TYPE = lambda x: unary_func(x, cas.arctanh)


def parse_quantities_in_struct(struct: Union[NestedMap[str], str])\
        -> Union[Quantity, NestedMap[Quantity]]:
    """Return a new struct that contains parsed quantities at the leaf
    values of the given input structure.

    The structure can be a nested dictionary, given the keys as strings and the
    leaf values as strings that can be parsed as ``pint`` quantities.
    For example:

    >>> from pprint import pprint
    >>> y = parse_quantities_in_struct({
    ...    'speed': {
    ...        'car': '100 km/hr',
    ...        'snail': '1 cm/min',
    ...        'fingernail': '1.2 mm/day'},
    ...    'weight': {
    ...        'car': '1.5 t',
    ...        'snail': '10 g',
    ...        'fingernail': '300 mg'}
    ... })
    >>> pprint(y)
    {'speed': {'car': <Quantity(100.0, 'kilometer / hour')>,
               'fingernail': <Quantity(1.2, 'millimeter / day')>,
               'snail': <Quantity(1.0, 'centimeter / minute')>},
     'weight': {'car': <Quantity(1.5, 'metric_ton')>,
                'fingernail': <Quantity(300, 'milligram')>,
                'snail': <Quantity(10, 'gram')>}}
    """
    try:
        items = struct.items()
    except AttributeError:
        return Quantity(struct)
    return {key: parse_quantities_in_struct(value) for key, value in items}


def quantity_dict_to_strings(struct: Quantity | NestedMap[Quantity],
                             significant_digits: int = 17) \
        -> str | NestedMap[str]:
    """Return a new structure with the quantity instances replaced by a string
    representation that is parsable by the :class:`simu.Quantity` constructor.

    Example:

    >>> from pprint import pprint
    >>> from simu import Quantity
    >>> struct = {'speed': {'car': Quantity(400 / 3, 'kilometer / hour'),
    ...                     'fingernail': Quantity(1.2, 'millimeter / day'),
    ...                     'snail': Quantity(1.0, 'centimeter / minute')},
    ...           'weight': {'car': Quantity(1.5, 'metric_ton'),
    ...                      'fingernail': Quantity(300, 'milligram'),
    ...                      'snail': Quantity(10, 'gram')}}
    >>> pprint(quantity_dict_to_strings(struct))
    {'speed': {'car': '133.33333333333334 km / h',
               'fingernail': '1.2 mm / d',
               'snail': '1 cm / min'},
     'weight': {'car': '1.5 t', 'fingernail': '300 mg', 'snail': '10 g'}}

    """
    try:
        items = struct.items()
    except AttributeError:
        return f"{struct:.{significant_digits}g~}"
    return {key: quantity_dict_to_strings(value, significant_digits)
            for key, value in items}




def extract_sub_structure(source: NestedMap[Quantity],
                          structure: NestedMap[str]) -> NestedMap[Quantity]:
    """Given a nested structure map ``structure`` that defines the units of
    measurement of leaf value quantities, extract those quantities from a source
    structure ``source``. A ``KeyError`` is raised if the source structure does
    not contain the requested data, and a ``DimensionalityError`` is raised if
    the queried quantity is of incompatible dimensions.

    Example:

    >>> from simu import Quantity
    >>> src = {"a": {"b": Quantity(1, "km")},
    ...        "c": Quantity(3, "degC"),
    ...        "d": {"e": Quantity(2, "s"), "f": Quantity(3, "kJ")}}
    >>> struct = {"a": {"b": "m"}, "d": {"e": "s"}}
    >>> print(extract_sub_structure(src, struct))
    {'a': {'b': <Quantity(1, 'kilometer')>}, 'd': {'e': <Quantity(2, 'second')>}}
    """
    def prepare(name: str, key: str, query: NestedMap[str],
                src: NestedMap[Quantity]) -> NestedMap[Quantity]:
        name = f"{name}.{key}" if name else key
        try:
            items = query.items()
        except AttributeError:
            try:
                src[key].to(query)
            except DimensionalityError as err:
                err.extra_msg = f" - Error fetching thermo parameter '{name}'."
                raise err from None
            return src[key]
        return {k: prepare(name, k, q, src[key]) for k, q in items}
    return {k: prepare("", k, s, source) for k, s in structure.items()}

