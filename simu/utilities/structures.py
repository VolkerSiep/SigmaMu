#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This module contains general helper functions that are useful on several
levels, while relying to maximal degree on standard python structures.
"""

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
    and values are the values of the leafs. Non-string keys will be converted
    to strings.

    :raises ValueError: If any of the keys contains the ``SEPARATOR`` character
    """
    try:
        items = structure.items()  # is this dictionary enough for us?
    except AttributeError:  # doesn't seem so
        return {prefix: structure}  # this is just a value

    result = {}
    for key, value in items:
        if SEPARATOR in str(key):
            raise ValueError(f"Separator '{SEPARATOR}' in key '{key}'")
        key = f"{prefix}{SEPARATOR}{key}" if prefix else str(key)
        result.update(flatten_dictionary(value, key))
    return result


def unflatten_dictionary(flat_structure):
    """This is the reverse of :func:`flatten_dictionary`, inflating the
    given one-depth dictionary into a nested structure."""
    result = {}

    def insert(struct, keys, value):
        first = keys.pop(0)
        if keys:
            if first not in struct:
                struct[first] = {}
            insert(struct[first], keys, value)
        else:
            struct[first] = value

    for key, value in flat_structure.items():
        insert(result, key.split(SEPARATOR), value)
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



# TODO: make below part of unit tests


def main():
    t = {"A":
         {"B": 1,
          "C": {"D": 2,
                "E": 3
              }
          },
         "F": 4,
         10: 5
         }

    flat = flatten_dictionary(t)
    print(unflatten_dictionrary(flat))


if __name__ == "__main__":
    main()
