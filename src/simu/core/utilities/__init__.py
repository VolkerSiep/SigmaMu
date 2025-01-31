"""The utilities module contains all required data structures and utility
functions that are not specific to particular objects, but have a certain
general purpose."""

from .structures import (
    MCounter, flatten_dictionary, unflatten_dictionary, nested_map)
from .quantity import (
    QFunction, Quantity, SymbolQuantity, base_magnitude, base_unit,
    conditional, jacobian, qpow, qvertcat, qsqrt, qsum, unit_registry,
    extract_units_dictionary, simplify_quantity)
from .qstructures import (
    ParameterDictionary, QuantityDict, parse_quantities_in_struct, exp, log,
    log10, sin, cos, tan, arcsin, arccos, arctan, sinh, cosh, tanh, arcsinh,
    arccosh, arctanh, sqrt, extract_sub_structure, quantity_dict_to_strings)
from .testing import assert_reproduction, user_agree
