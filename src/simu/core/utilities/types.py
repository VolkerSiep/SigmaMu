"""This module defines types of complex data structures"""

from typing import Union, TypeVar
from collections.abc import Mapping, MutableMapping

__V = TypeVar("__V")

# hope this is simpler in the future
#  for now necessary to avoid type warnings
Map = Mapping[str, __V]
__NestedMap_1 = Mapping[str, Union[__V, Map]]
__NestedMap_2 = Mapping[str, Union[__V, __NestedMap_1]]
__NestedMap_3 = Mapping[str, Union[__V, __NestedMap_2]]
__NestedMap_4 = Mapping[str, Union[__V, __NestedMap_3]]
__NestedMap_5 = Mapping[str, Union[__V, __NestedMap_4]]
__NestedMap_6 = Mapping[str, Union[__V, __NestedMap_5]]
__NestedMap_7 = Mapping[str, Union[__V, __NestedMap_6]]
__NestedMap_8 = Mapping[str, Union[__V, __NestedMap_7]]
__NestedMap_9 = Mapping[str, Union[__V, __NestedMap_8]]
NestedMap = Mapping[str, Union[__V, __NestedMap_9]]


MutMap = MutableMapping[str, __V]
__NestedMutMap_1 = MutableMapping[str, Union[__V, MutMap]]
__NestedMutMap_2 = MutableMapping[str, Union[__V, __NestedMutMap_1]]
__NestedMutMap_3 = MutableMapping[str, Union[__V, __NestedMutMap_2]]
__NestedMutMap_4 = MutableMapping[str, Union[__V, __NestedMutMap_3]]
__NestedMutMap_5 = MutableMapping[str, Union[__V, __NestedMutMap_4]]
__NestedMutMap_6 = MutableMapping[str, Union[__V, __NestedMutMap_5]]
__NestedMutMap_7 = MutableMapping[str, Union[__V, __NestedMutMap_6]]
__NestedMutMap_8 = MutableMapping[str, Union[__V, __NestedMutMap_7]]
__NestedMutMap_9 = MutableMapping[str, Union[__V, __NestedMutMap_8]]
NestedMutMap = MutableMapping[str, Union[__V, __NestedMutMap_9]]


# NestedStringDict = MutableMapping[str, Union[str, "NestedStringDict"]]

# QuantityDict = MutableMapping[str, "Quantity"]
# NestedQuantityDict = \
#     MutableMapping[str, Union["Quantity", "NestedQuantityDict"]]
