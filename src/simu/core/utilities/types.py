"""This module defines types / base classes of complex data structures"""

# stdlib
from typing import TypeVar, Generic
from collections.abc import Mapping, MutableMapping

T = TypeVar("T")
"""A generic type variable"""

Map = Mapping[str, T]
"""A mapping of strings to another type"""

NestedMap = Map[T | "NestedMap[T]"]
"""A nested mapping of strings to another type"""

MutMap = MutableMapping[str, T]
"""A mutable mapping of strings to another type"""

NestedMutMap = MutMap[T | "NestedMutMap[T]"]
"""A nested mutable mapping of strings to another type"""
