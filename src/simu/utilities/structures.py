"""
This module contains general helper functions that are useful on several
levels, while relying to maximal degree on standard python structures.
"""

# stdlib modules
from typing import Tuple
from collections import Counter
from collections.abc import Iterable

# internal modules
from .quantity import base_unit, qvertcat, SymbolQuantity


class MCounter(Counter):
    """This is a slight extention of the ``Collections.Counter`` class
    to also allow multiplication with integers:

        >>> a = MCounter({"a": 1})
        >>> b = MCounter({"b": 1})
        >>> a + 2 * b
        MCounter({'b': 2, 'a': 1})
    """

    def __mul__(self, other):
        if not isinstance(other, int):
            raise TypeError("Non-int factor")
        return MCounter({k: other * v for k, v in self.items()})

    def __rmul__(self, other):
        return self * other  # call __mul__

    def __add__(self, other):
        return MCounter(super().__add__(other))

    def __pos__(self):
        return self

    # I don't think I need this!?
    # @classmethod
    # def fromkeys(cls, iterable, v=None):
    #     raise NotImplementedError()


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
