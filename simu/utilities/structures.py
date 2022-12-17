"""
This module contains general helper functions that are useful on several
levels, while relying to maximal degree on standard python structures.
"""

# stdlib modules
from collections import Counter
from collections.abc import Mapping
from typing import Collection

# external modules
from casadi import SX, vertcat

# internal modules
from .quantity import (base_unit, flatten_dictionary, qvertcat,
                       unflatten_dictionary, SymbolQuantity)


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

    @classmethod
    def fromkeys(cls, iterable, v=None):
        raise NotImplementedError()


class ParameterDictionary(dict):
    """This class is a nested dictionary of SymbolQuantities to represent
    parameters with functionality to be populated using the ``register_*``
    methods.
    """

    class SparseMatrix(dict):
        """This helper class represents a nested dictionary that contains
        two levels of keys and values representing a quantity."""

        def pair_items(self):
            """Return a scalar quantity with the key pair for each element
            in the sub-structure"""
            for key_1, second in self.items():
                for key_2, quantity in second.items():
                    yield (key_1, key_2, quantity)

    def register_scalar(self, key: str, unit: str):
        """Create a scalar quantity and add the structure to the dictionary"""
        unit = base_unit(unit)
        quantity = SymbolQuantity(key, unit)
        self[key] = quantity
        return quantity

    def register_vector(self, key: str, species: list[str], unit: str):
        """Create a quantity vector with symbols and add the structure to
        the dictionary"""
        unit = base_unit(unit)
        self[key] = {s: SymbolQuantity(f"{key}.{s}", unit) for s in species}
        return qvertcat(*self[key].values())

    def register_sparse_matrix(self, key: str, pairs: dict, unit: str):
        """Create a sparse matrix quantity and add the structure to the
        dictionary. The given unit is converted to base units before being
        applied, i.e. """
        unit = base_unit(unit)
        res = ParameterDictionary.SparseMatrix({f: {} for f, _ in pairs})
        for first, second in pairs:
            quantity = SymbolQuantity(f"{key}.{first}.{second}", unit)
            res[first][second] = quantity
        self[key] = res
        return res

    def get_quantity(self, *keys):
        """Extract a quantity from the given sequence of keys"""
        entry = self
        for key in keys:
            entry = entry[key]
        return entry

    def get_vector_quantity(self, *keys):
        """Extract a vector quantity from the given sequence of keys"""
        entry = self
        for key in keys:
            entry = entry[key]
        return qvertcat(*entry.values())


class FlexiDict(Mapping):  # maybe not needed anymore
    """Using Casadi, we often deal with nested dictionaries that hold both
    casadi symbols of type ``SX`` and floating point values. These structures
    need to be used both in flat and in nested configuration.

    The purpose of this class is to describe such structures, being able to
    switch between the representations (flat or not, and symbolic or not).
    """

    def __init__(self, dictionary: dict, flat=True, symbol=False):
        """Create a flexible dictionary from a source.

        :param dictionary: The data source, either a flat or a nested
            dictionary, and either with float or with casadi symbol values.
        :param flat: Flag being ``True`` if the provided dictionary is flat, and
            ``False`` if the dictionary is nested.
        :param symbol: Flag being ``True`` if the provided dictionary has values
            of type ``casadi.SX``, and ``False``, if floats are provided.

        The internal representation is always a flattened dictionary with
        a list of floats and a vector-shaped ``casadi.SX`` object.
        """
        self._is_flat = flat
        self._is_symbol = symbol
        if not flat:
            dictionary = flatten_dictionary(dictionary)

        if symbol:
            data = {
                "keys": dictionary.keys(),
                "symbols": vertcat(*dictionary.values()),
                "values": [0.0] * len(dictionary),
            }
        else:
            data = {
                "keys": dictionary.keys(),
                "symbols": SX.sym("param", len(dictionary)),
                "values": list(dictionary.values()),
            }
        self._data = data

    def view(self, flat: bool = None, symbol: bool = None):
        """Changes the view of the dictionary.

        :param flat: ``True`` if the dictionary should appear flat.
        :param symbol: ``True`` if the provided dictionary should appear
            with casadi symbols as values.

        If any parameter is left at ``None``, the original view attribute is
        kept. Switching views is not changing the internal data representation.
        """
        flat = self._is_flat if flat is None else flat
        symbol = self._is_symbol if symbol is None else symbol
        result = FlexiDict({}, flat=flat, symbol=symbol)
        result._data = self._data  # pylint: disable=protected-access
        return result

    def __getitem__(self, key):
        return self._dict[key]

    def __len__(self):
        return len(self._dict)

    def __iter__(self):
        return iter(self._dict)

    @property
    def symbols(self) -> SX:
        """The flat list of symbols"""
        return self._data["symbols"]

    @symbols.setter
    def symbols(self, symbols: SX):
        if symbols.rows() != len(self._data["keys"]):
            raise ValueError("Incompatible length of values")
        if symbols.columns() != 1:
            raise ValueError("Only one-column symbol vectors allowed")
        self._data["symbols"] = symbols

    @property
    def flat_values(self) -> list[float]:
        """The flat list of numeric values"""
        return self._data["values"]

    @flat_values.setter
    def flat_values(self, values: Collection[float]):
        if len(values) != len(self._data["keys"]):
            raise ValueError("Incompatible length of values")
        self._data["values"] = values

    def set_struct_values(self, struct: dict):
        """Set float values, given in a possibly nested dictionary.

        :param struct: A dictionary mapping the keys of the dictionary
          to the values to be set. The given dictionary must be complete and not
          contain keys that are not part of this instance.
        """
        flat = flatten_dictionary(struct)
        keys = flat.keys()
        if len(keys) != len(self._data["keys"]):
            raise ValueError("Unequal length of parameter structures")
        for n_1, n_2 in zip(keys, self._data["keys"]):
            if n_1 != n_2:
                raise ValueError(f"Unequal parameter name: '{n_1}' vs '{n_2}'")
        self._data["values"] = list(flat.values())

    @property
    def _dict(self):
        if self._is_symbol:
            val = self._data["symbols"].elements()
        else:
            val = self._data["values"]
        result = dict(zip(self._data["keys"], val))
        if not self._is_flat:
            result = unflatten_dictionary(result)
        return result


class SpeciesDict(dict):  # better as a global function?
    """A specialised dictionary that can be reduced to only host a subset of
    keys. This is primarily used to deal with data sets containing unary,
    binary and higher order data on a set of chemical species."""

    def __init__(self, *args, **kwargs):
        dict.__init__(self, *args, **kwargs)

    def reduce(self, species: Collection[str]):
        """Perform a shallow reduction on this instance.

        :param species: The list of species to reduce the structure to.
        :return: A copy of this dictionary, but only with the keys listed in
          ``species``"""
        return SpeciesDict(
            {key: value
             for key, value in self.items() if key in species})

    @staticmethod
    def deep_reduce(dictionary: dict, species: Collection[str]):
        """Reduce the species of the given dictionary recursively.

        If the dictionary itself is a ``SpeciesDict``, apply :py:meth:`reduce`.
        Then (independent on above condition), apply reduction to all values of
        ``dictionary`` recursively.
        """
        # is this a species dict?
        try:
            dictionary = dictionary.reduce(species)  # type: ignore
        except AttributeError:  # no, just keep it
            dictionary = dict(dictionary)  # make a shallow copy

        # is this still iterable?
        try:
            items = dictionary.items()  # yes, continue below
        except AttributeError:  # no, return object as is
            return dictionary

        return {
            key: SpeciesDict.deep_reduce(value, species)
            for key, value in items
        }


def create_vector_struct(keys, unit):
    """Create a parameter structure with ``None`` values for given keys and
    unit of measurement. For instance

    .. code-block::

        >>> create_vector_struct(["A", "B", "C"], "J/mol")
        [{'A': None, 'B': None, 'C': None}, 'J/mol']

    """
    return [{k: None for k in keys}, unit]


def create_scalar_struct(unit):
    """Create a scalar parameter structure with ``None`` values for a given
    unit of measurement. For instance

    .. code-block::

        >>> create_scalar_struct("J/mol")
        [None, 'J/mol']
    """
    return [None, unit]


def create_vector(dictionary: dict, keys: Collection[str]) -> SX:
    """During :meth:`define`, sub-directories often require to be converted
    into a single ``casadi.SX`` vector. This method provides such
    functionality.

    .. code-block::

        >>> struct = {"A": SX.sym("x"),
        ...           "C": SX.sym("y"),
        ...           "B": SX.sym("z")}
        >>> create_vector(struct, ["A", "B", "C"])
        SX([x, z, y])

    :param dictionary: The parameter structure, pointing to the node where
        to extract the vector from
    :param keys: The keys of the sub-dictionary to collect the ``casadi``
        symbols from. If left as ``None``, the species names of the
        contribution are assumed.
    :return: The ``casadi.SX`` object containing the collected symbols
    """
    return vertcat(*[dictionary[key] for key in keys])
