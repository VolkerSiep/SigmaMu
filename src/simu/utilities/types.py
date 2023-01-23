"""This module defines types of complex data structures"""

from typing import Union, TypeVar, Tuple, Type
from collections.abc import MutableMapping

StringDict = MutableMapping[str, str]
NestedStringDict = MutableMapping[str, Union[str, "NestedStringDict"]]

QuantityDict = MutableMapping[str, "Quantity"]
NestedQuantityDict = \
    MutableMapping[str, Union["Quantity", "NestedQuantityDict"]]

T = TypeVar("T")
TDict = MutableMapping[str, T]
NestedTDict = MutableMapping[str, Union[T, "NestedTDict"]]

ThermoContributionDict = \
    MutableMapping[str, Tuple[Type["Contribution"], MutableMapping]]
