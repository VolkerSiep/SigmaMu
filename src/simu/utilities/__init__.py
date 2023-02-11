"""The utilities module contains all required data structures and utility
functions that are not specific to particular objects, but have a certain
general purpose."""

from .structures import (
    MCounter, flatten_dictionary, unflatten_dictionary)
from .quantity import (
    QFunction, Quantity, SymbolQuantity, base_magnitude, base_unit,
    conditional, jacobian, qpow, qvertcat, sqrt, sum1, unit_registry,
    extract_units_dictionary)
from .qstructures import ParameterDictionary, exp, log
from .testing import assert_reproduction, user_agree
