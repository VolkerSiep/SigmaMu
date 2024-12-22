"""This module defines types of complex data structures"""

from typing import Union, TypeVar
from collections.abc import Mapping, MutableMapping

__V = TypeVar("__V")

Map = Mapping[str, __V]
MutMap = MutableMapping[str, __V]
NestedMap = Map[__V | "NestedMap"]
NestedMutMap = MutMap[__V | "NestedMutMap"]
