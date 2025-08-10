from logging import getLogger, NullHandler

from .core.utilities.quantity import (
    Quantity, SymbolQuantity, base_unit, base_magnitude, QFunction,
    simplify_quantity, qsum, qvertcat, qpow, conditional, jacobian,
    unit_registry, extract_units_dictionary)
from .core.utilities.qstructures import (
    QuantityDict, extract_sub_structure, parse_quantities_in_struct,
    quantity_dict_to_strings, log, log10, sqrt, exp, sin, cos, tan, arcsin,
    arccos, arctan, sinh, cosh, tanh, arcsinh, arccosh, arctanh,
    ParameterDictionary)
from .core.utilities.structures import (
    flatten_dictionary, unflatten_dictionary, MCounter)
from .core.utilities.constants import (
    PI, ALPHA, EPS_0, E_0, F, GAMMA_G, H_PLANCK, K_B, MU_0, M_E, M_N, M_P, N_A,
    R_B, R_GAS, R_INF, SIGMA, STD_GRAVITY, V_LIGHT
)

from .core.thermo.species import SpeciesDefinition, SpeciesDB
from .core.thermo.state import (
    InitialState, StateDefinition, registered_state, all_states)
from .core.thermo.contribution import (
    ThermoContribution, registered_contribution, all_contributions)
from .core.thermo.frame import ThermoFrame
from .core.thermo.factory import ThermoFactory
from .core.thermo.parameters import (
    ThermoParameterStore, AbstractThermoSource, NestedDictThermoSource,
    StringDictThermoSource)
from .core.thermo.material import MaterialDefinition, Material, MaterialSpec

from .core.model.base import Model
from .core.model.amodel import AModel
from .core.model.numeric import NumericHandler
from .core.solver.simulation import SimulationSolver

# versioning
from ._version import version as __version__

# logging
logger = getLogger(__name__)
logger.addHandler(NullHandler())
