"""
This module contains general helper functions that are useful on several
levels, while relying to maximal degree on standard python structures.
"""
from re import escape, split
from typing import TypeVar, Callable
from collections import Counter

# internal modules
from .types import NestedMap, MutMap, Map, NestedMutMap

_V = TypeVar("_V")
_R = TypeVar("_R")
FLATTEN_SEPARATOR = "/"  # separator when (un-)flattening dictionaries


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


def flatten_dictionary(structure: NestedMap[_V], prefix: str = "") -> Map[_V]:
    r"""Convert the given structure into a flat list of key value pairs,
    where the keys are ``SEPARATOR``-separated concatonations of the paths,
    and values are the values of the leafs. Non-string keys are converted
    to strings. Occurances of ``SEPARATOR`` are escaped by ``\``.

    >>> d: NestedMap[int] = {"a": {"b": 1, "c": 2}, "d": {"e/f": 3}}
    >>> flatten_dictionary(d)
    {'a/b': 1, 'a/c': 2, 'd/e\\/f': 3}
    """
    try:
        items = structure.items()  # is this dictionary enough for us?
    except AttributeError:  # doesn't seem so, this is just a value
        return {prefix: structure}  # type: ignore

    result: MutMap[_V] = {}
    # must sort to create the same sequence every time
    # (dictionary might have content permuted)
    for key, value in sorted(items):
        key = str(key).replace(FLATTEN_SEPARATOR, rf"\{FLATTEN_SEPARATOR}")  # esc. separator
        key = f"{prefix}{FLATTEN_SEPARATOR}{key}" if prefix else key
        result.update(flatten_dictionary(value, key))
    return result


def unflatten_dictionary(flat_structure: Map[_V]) -> NestedMap[_V]:
    r"""This is the reverse of :func:`flatten_dictionary`, inflating the
    given one-depth dictionary into a nested structure.

    >>> d = {"a/b": 1, "a/c": 2, r"d/e\/f": 3}
    >>> unflatten_dictionary(d)
    {'a': {'b': 1, 'c': 2}, 'd': {'e/f': 3}}
    """
    result: NestedMutMap[_V] = {}

    def insert(struct, sub_keys, sub_value):
        """insert one element into the nested structure"""
        first = sub_keys.pop(0)
        if sub_keys:
            if first not in struct:
                struct[first] = {}
            insert(struct[first], sub_keys, sub_value)
        else:
            struct[first] = sub_value

    # split by non-escaped separators and unescape escaped separators
    regex = rf'(?<!\\){escape(FLATTEN_SEPARATOR)}'
    for key, value in flat_structure.items():
        keys = [
            k.replace(rf"\{FLATTEN_SEPARATOR}", FLATTEN_SEPARATOR)
            for k in split(regex, key)
        ]
        insert(result, keys, value)
    return result


def nested_map(structure: NestedMap[_V],
               function: Callable[[_V], _R]) -> NestedMap[_R]:
    """Apply a unary function to each leaf values of the given nested
    dictionary, and return the same structure with the function's values as
    leafs.
    >>> a = {'a': {'b': 1, 'c': 2}, 'd': {'e/f': 3}}
    >>> nested_map(a, lambda x: 2 * x)
    {'a': {'b': 2, 'c': 4}, 'd': {'e/f': 6}}
    """
    try:
        items = structure.items()
    except AttributeError:
        return function(structure)
    return {k: nested_map(v, function) for k, v in items}