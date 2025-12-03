"""This module defines types of complex data structures"""

# stdlib
from typing import TypeVar, Self, TypeAliasType
from collections.abc import Mapping, MutableMapping

VT = TypeVar("VT")
"""An arbitrary value type"""

Map = Mapping[str, VT]
"""A mapping of strings to another type"""

MutMap = MutableMapping[str, VT]
"""A mutable mapping of strings to another type"""

NestedMap = TypeAliasType("NestedMap", Map[VT | "NestedMap[VT]"],
                          type_params=(VT, ))
"""A nested mapping of strings to another type"""

NestedMutMap = TypeAliasType("NestedMutMap", MutMap[VT | "NestedMutMap[VT]"],
                             type_params=(VT, ))
"""A nested mutable mapping of strings to another type"""
