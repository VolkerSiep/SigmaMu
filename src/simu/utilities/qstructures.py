"""This module contains data structures that build on the quantity datatype"""

# stdlibs
from typing import Iterable, Tuple, Union
from collections.abc import Callable

# external libs
import casadi as cas
from pint.errors import DimensionalityError

# internal libs
from . import base_unit, Quantity, SymbolQuantity, qvertcat


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

    def register_vector(self, key: str, sub_keys: Iterable[str], unit: str):
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

    def register_sparse_matrix(self, key: str,
                               pairs: Iterable[Tuple[str, str]], unit: str):
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


_OType = Union[float, Quantity, "QuantityDict"]
_RType = Union[Quantity,"QuantityDict"]
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
                             keys: list[str]) -> "QuantityDict":
        """As the magnitude of a ``Quantity`` can be a container itself,
        this convenience factory method combines such a vector quantity with
        a set of given keys into a ``QuantityDict`` object

        >>> a = Quantity([1,2,3], "m")
        >>> y = QuantityDict.from_vector_quantity(a, ["A", "B", "C"])
        >>> for key, value in y.items(): print(f"{key}: {value:~}")
        A: 1 m
        B: 2 m
        C: 3 m
        """
        return cls({key: quantity[k] for k, key in enumerate(keys)})

    def __add__(self, other: _OType) -> "QuantityDict":
        try:
            items = other.items()
        except AttributeError:
            return QuantityDict({k: v + other for k, v in self.items()})

        result = self.copy()
        for key, value in items:
            result[key] = (result[key] + value) if key in self else value
        return QuantityDict(result)

    def __radd__(self, other: _OType) -> "QuantityDict":
        return self + other

    def __mul__(self, other: _OType) -> "QuantityDict":
        try:
            other.items()
        except AttributeError:
            return QuantityDict({k: v * other for k, v in self.items()})

        result = {k: v * other[k] for k, v in self.items() if k in other}
        return QuantityDict(result)

    def __rmul__(self, other: _OType) -> "QuantityDict":
        return self * other

    def __pos__(self) -> "QuantityDict":
        return self

    def __neg__(self) -> "QuantityDict":
        return QuantityDict({k: -v for k, v in self.items()})

    def __sub__(self, other: _OType) -> "QuantityDict":
        return self + (-other)

    def __rsub__(self, other: _OType) -> "QuantityDict":
        return (-self) + other

    def __truediv__(self, other: _OType) -> "QuantityDict":
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

    def __rtruediv__(self, other: _OType) -> "QuantityDict":
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

# TODO: exponential operator, dot product @, unary functioons, ...




def unary_func(quantity: _OType, func: Callable[[_SType], _SType]) -> _RType:
    """Call a unary function that requires a dimensionless argument on the
    argument, accepting that argument to be a float, a scalar quantity,
    a symbolic quantity, or a QuantityDict object (well, any dictionary
    of the above, actually."""
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


def log(quantity: _OType) -> _RType:
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


log10 = lambda x: unary_func(x, cas.log10)
exp = lambda x: unary_func(x, cas.exp)
sin = lambda x: unary_func(x, cas.sin)
cos = lambda x: unary_func(x, cas.cos)
tan = lambda x: unary_func(x, cas.tan)
arcsin = lambda x: unary_func(x, cas.arcsin)
arccos = lambda x: unary_func(x, cas.arccos)
arctan = lambda x: unary_func(x, cas.arctan)
sinh = lambda x: unary_func(x, cas.sinh)
cosh = lambda x: unary_func(x, cas.cosh)
tanh = lambda x: unary_func(x, cas.tanh)
arcsinh = lambda x: unary_func(x, cas.arcsinh)
arccosh = lambda x: unary_func(x, cas.arccosh)
arctanh = lambda x: unary_func(x, cas.arctanh)
