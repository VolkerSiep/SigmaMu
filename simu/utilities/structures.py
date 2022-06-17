#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This module contains general helper functions that are useful on several
levels, while relying to maximal degree on standard python structures.
"""

# stdlib modules
from typing import Collection
from collections.abc import Mapping

# external modules
from casadi import DM, SX, vertcat

SEPARATOR = "."


def iter_binary_parameters(species, parameters, name):
    """For a sparse binary parameter structure in nested dictionaries, this is
    a generator to yield tuples of indices and values."""
    try:
        current = parameters[name].items()
    except KeyError:
        return
    for first, rest in current:
        idx_1 = species.index(first)
        for second, value in rest.items():
            idx_2 = species.index(second)
            yield (idx_1, idx_2, value)


def flatten_dictionary(structure, prefix=''):
    r"""Convert the given structure into a flat list of key value pairs,
    where the keys are ``SEPARATOR``-separated concatonations of the paths,
    and values are the values of the leafs. Non-string keys are converted
    to strings. Occurances of ``SEPARATOR`` are escaped by ``\``.
    """
    try:
        items = structure.items()  # is this dictionary enough for us?
    except AttributeError:  # doesn't seem so
        return {prefix: structure}  # this is just a value

    result = {}
    # must sort to create the same sequence every time
    # (dictionary might have content permutated)
    for key, value in sorted(items):
        key = str(key).replace(".", r"\.")
        key = f"{prefix}{SEPARATOR}{key}" if prefix else key
        result.update(flatten_dictionary(value, key))
    return result


def unflatten_dictionary(flat_structure):
    """This is the reverse of :func:`flatten_dictionary`, inflating the
    given one-depth dictionary into a nested structure."""
    result = {}

    def insert(struct, keys, value):
        first = keys.pop(0)
        while first.endswith("\\"):
            first = f"{first[:-1]}.{keys.pop(0)}"
        if keys:
            if first not in struct:
                struct[first] = {}
            insert(struct[first], keys, value)
        else:
            struct[first] = value

    for key, value in flat_structure.items():
        insert(result, key.split(SEPARATOR), value)
    return result


class FlexiDictionary(Mapping):
    def __init__(self, dictionary: dict,
                 values_are_symbols: bool =False,
                 flat: bool = False,
                 symbol: bool = False):
        self._is_flat = flat
        self._is_symbol = symbol
        flat_dict = flatten_dictionary(dictionary)
        if values_are_symbols:
            data = {"keys": flat_dict.keys(),
                    "symbols": vertcat(*flat_dict.values()),
                    "values":  [0.0] * len(flat_dict)}
        else:
            data = {"keys": flat_dict.keys(),
                    "symbols": SX.sym("param", len(flat_dict)),
                    "values":  list(flat_dict.values())}
        self._data = data

    def view(self, flat: bool = True, symbol: bool = False):
        result = FlexiDictionary({}, flat=flat, symbol = symbol)
        result._data = self._data  # copy reference to data only, ...
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
        """The flat list of numeric values """
        return self._data["values"]

    @flat_values.setter
    def flat_values(self, values: Collection[float]):
        if len(values) != len(self._data["keys"]):
            raise ValueError("Incompatible length of values")
        self._data["values"] = values

    def set_struct_values(self, struct: dict):
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


class SpeciesDict(dict):
    def __init__(self, *args, **kwargs):
        dict.__init__(self, *args, **kwargs)

    def filter(self, species):
        return SpeciesDict({key: value for key, value in self.items()
                            if key in species})

    @staticmethod
    def filter_species(dictionary, species):
        # is this a species dict?
        try:
            dictionary = dictionary.filter(species)
        except AttributeError:  # no, just keep it
            dictionary = dict(dictionary)  # make a shallow copy

        # is this still iterable?
        try:
            items = dictionary.items()  # yes, continue below
        except AttributeError:  # no, return object as is
            return dictionary

        return {key: SpeciesDict.filter_species(value, species)
                for key, value in items}



# TODO: make below part of unit tests and test all of above!!


def main():
    t = {"A":
         {"B.C": 1,
          "C": {"D": 2,
                "E": 3
              }
          },
         "F": 4,
         "10": 5
         }

    flat = flatten_dictionary(t)
    print(unflatten_dictionary(flat))


if __name__ == "__main__":
    main()
