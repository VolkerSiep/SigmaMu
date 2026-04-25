"""This module defines types of complex data structures"""

from collections.abc import Mapping, MutableMapping
from typing import Protocol, runtime_checkable

type Map[T] = Mapping[str, T]
"""A mapping of strings to another type"""

type NestedMap[T] = Map[T | "NestedMap[T]"]
"""A nested mapping of strings to another type"""

type MutMap[T] = MutableMapping[str, T]
"""A mutable mapping of strings to another type"""

type NestedMutMap[T] = MutMap[T | "NestedMutMap[T]"]
"""A nested mutable mapping of strings to another type"""


@runtime_checkable
class OutputIOStream(Protocol):
    """A Protocol class to describe anything that can write a string"""
    def write(self, line: str):  ...

