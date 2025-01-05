"""This module defines types of complex data structures"""

from typing import Union, TypeVar, Self
from collections.abc import Mapping, MutableMapping

__V = TypeVar("__V")
"""An arbitrary type"""

Map = Mapping[str, __V]
"""A mapping of strings to another type"""

MutMap = MutableMapping[str, __V]
"""A mutable mapping of strings to another type"""

NestedMap = Map[__V | Self]
"""A nested mapping of strings to another type"""

NestedMutMap = MutMap[__V | Self]
"""A nested mutable mapping of strings to another type"""
