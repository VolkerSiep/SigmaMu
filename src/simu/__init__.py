# -*- coding: utf-8 -*-
from logging import getLogger, NullHandler

from simu.core.utilities import (
    ParameterDictionary, base_magnitude, log, qsum, conditional, jacobian,
    qpow, qvertcat, sqrt, exp, log, log10, sin, cos, tan, arcsin, arccos,
    arctan, sinh, cosh, tanh, arcsinh, arccosh, arctanh, base_magnitude,
    Quantity, QuantityDict, MCounter, SymbolQuantity, QFunction, base_unit,
    flatten_dictionary, unflatten_dictionary, extract_units_dictionary
)

from simu.core.utilities.constants import (
    PI, R_GAS, V_LIGHT, H_PLANCK, GAMMA_G, N_A, ALPHA,
    SIGMA, F, EPS_0, MU_0, E_0, K_B, STD_GRAVITY
)

from simu.core.thermo import (
    all_states, InitialState, SpeciesDefinition, SpeciesDB, StateDefinition,
    ThermoContribution, ThermoFactory, ThermoFrame
)

from simu.core.model import Model, NumericHandler

# versioning
from ._version import VERSION as __version__

# logging
logger = getLogger(__name__)
logger.addHandler(NullHandler())
