"""
This module contains general helper functions that are useful on several
levels, while relying to maximal degree on standard python structures.
"""
# stdlib
from re import escape, split
from collections import Counter
from collections.abc import Mapping, Callable, Sequence, MutableMapping
from typing import Iterator

# external modules
from numpy import hstack

# internal modules
from .types import NestedMap, MutMap, Map, NestedMutMap

FLATTEN_SEPARATOR = "/"  # separator when (un-)flattening dictionaries


class MCounter(Counter):
    """This is a slight extension of the ``Collections.Counter`` class
    to also allow multiplication with scalar numbers:

    >>> a = MCounter({"a": 1})
    >>> b = MCounter({"b": 1})
    >>> a + 2.5 * b
    MCounter({'b': 2.5, 'a': 1})

    Note that we stretch the use of MCounter to allow floats. Therefore the
    ``elements`` method is removed as a compromise in design.
    """

    def __mul__(self, other):
        return MCounter({k: other * v for k, v in self.items()})

    def __rmul__(self, other):
        return self * other  # call __mul__

    def __add__(self, other):
        return MCounter(super().__add__(other))

    def __pos__(self):
        return self

    def elements(self):
        raise NotImplemented


class FlatDictionary[V](MutableMapping[str, V]):
    """Effective wrapping of a nested dictionary to access flat value array
    without overhead, but still allowing to set new values:

    >>> source = {"a": [1, 2, 3], "b": [4, 5], "c": 6}
    >>> d = FlatDictionary(source)
    >>> d["a"]
    array([1, 2, 3])
    >>> d.values
    array([1, 2, 3, 4, 5, 6])

    >>> d["b"] = [7, 8]
    >>> d.values
    array([1, 2, 3, 7, 8, 6])
    """
    def __init__(self, source: Map[V]):
        flat = flatten_dictionary(source)
        start = 0
        self._index = {key: (start, start := start + general_length(value))
                       for key, value in flat.items()}
        self.values = hstack(list(flat.values()))

    def __setitem__(self, key: str, value: V):
        s, e = self._index[key]
        assert general_length(value) == e - s
        self.values[s: e] = value

    def __delitem__(self, key: str):
        raise RuntimeError("Deletion from FlatDictionary is not permitted")

    def __getitem__(self, key: str) -> V:
        s, e = self._index[key]
        return self.values[s: e]

    def __len__(self) -> int:
        return len(self._index)

    def __iter__(self) -> Iterator[str]:
        return iter(self._index)


def general_length(arr) -> int:
    try:
        return len(arr)
    except TypeError:
        return 1


def flatten_dictionary[V](structure: NestedMap[V], prefix: str = "") -> Map[V]:
    r"""Convert the given structure into a flat list of key value pairs,
    where the keys are ``SEPARATOR``-separated concatenations of the paths,
    and values are the values of the leafs. Non-string keys are converted
    to strings. Occurrences of ``SEPARATOR`` are escaped by ``\``.

    >>> d: NestedMap[int] = {"a": {"b": 1, "c": 2}, "d": {"e/f": 3}}
    >>> flatten_dictionary(d)
    {'a/b': 1, 'a/c': 2, 'd/e\\/f': 3}
    """
    if not isinstance(structure, Mapping):
        return {prefix: structure}

    result: MutMap[V] = {}
    # We must sort to create the same sequence each time
    #   (dictionary might have content permuted)
    for key, value in sorted(structure.items()):
        # esc. separator
        key = str(key).replace(FLATTEN_SEPARATOR, rf"\{FLATTEN_SEPARATOR}")
        key = f"{prefix}{FLATTEN_SEPARATOR}{key}" if prefix else key
        result.update(flatten_dictionary(value, key))
    return result


def unflatten_dictionary[V](flat_structure: Map[V]) -> NestedMutMap[V]:
    r"""This is the reverse of :func:`flatten_dictionary`, inflating the
    given one-depth dictionary into a nested structure.

    >>> d = {"a/b": 1, "a/c": 2, r"d/e\/f": 3}
    >>> unflatten_dictionary(d)
    {'a': {'b': 1, 'c': 2}, 'd': {'e/f': 3}}
    """
    result: NestedMutMap[V] = {}

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


def path_to_flat_key(path: Sequence[str]) -> str:
    """Convert a path, representing the key sequence of the non-flattened
    dictionary, into the equivalent key that can be used to address the same
    element in the flattened dictionary
    """
    path = [p.replace(FLATTEN_SEPARATOR, rf"\{FLATTEN_SEPARATOR}")
            for p in path]
    return FLATTEN_SEPARATOR.join(path)


def nested_map[V, R](structure: NestedMap[V],
               function: Callable[[V], R]) -> NestedMap[R]:
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