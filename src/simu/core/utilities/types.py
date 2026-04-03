"""This module defines types / base classes of complex data structures"""

# stdlib
# from typing import TypeVar
from collections.abc import Mapping, MutableMapping

# T = TypeVar("T")
# """A generic type variable"""

type Map[T] = Mapping[str, T]
"""A mapping of strings to another type"""

type NestedMap[T] = Map[T | "NestedMap[T]"]
"""A nested mapping of strings to another type"""

type MutMap[T] = MutableMapping[str, T]
"""A mutable mapping of strings to another type"""

type NestedMutMap[T] = MutMap[T | "NestedMutMap[T]"]
"""A nested mutable mapping of strings to another type"""
