from pathlib import Path
from collections.abc import Sequence, Callable, Mapping
from numpy import linspace, exp
from yaml import safe_load
from simu import NumericHandler, Model, InitialState

from model import MuFitModel
from thermo import lp_condensate, condensate_new, lp_steam, steam_new

THIS_DIRECTORY = Path(__file__).parent
THERMO_FIT_DEFINITION_FILE = THIS_DIRECTORY / "thermo_fit_definition.yml"

def p_sat(t: float) -> float:
    """Buck equation, https://en.wikipedia.org/wiki/Vapour_pressure_of_water"""
    return 6.1121e-3 * exp((18.678 - t / 234.5) * (t / (257.14 + t)))

def is_liquid(t: float, p: float) -> bool:
    return p_sat(t) <= p

def is_gas(t: float, p: float) -> bool:
    return p_sat(t) >= p

def gen_tp(t_space: Sequence[float], p_space: Sequence[float],
           valid: Callable[[float, float], bool]) -> Mapping:
    return {
        "columns": ["T", "p"],
        "uom": ["degC", "bar"],
        "data": [[float(t_i), float(p_i)]
                 for t_i in t_space for p_i in p_space
                 if valid(t_i, p_i)]
    }

def load_definition():
    with THERMO_FIT_DEFINITION_FILE.open() as f:
        data = safe_load(f)

    t_range = linspace(0, 150, num=20)
    p_range = [0.05, 0.1, 1.0, 2.0, 5.0, 10.0, 30.0]

    # create datasets section
    data["datasets"] = {
        "condensate": gen_tp(t_range, p_range, is_liquid),
        "steam": gen_tp(t_range, p_range, is_gas)
    }
    return data

def define_models():
    def nh(model: Model) -> NumericHandler:
        return NumericHandler(model.create_proxy().finalise())
    lp_condensate.initial_state = InitialState.from_cbar(0, 0.01, [1.0])
    # lp_steam.initial_state = InitialState.from_cbar(5, 0.1, [1.0])

    return {
        "mu_cmp_condensate": nh(MuFitModel(lp_condensate, condensate_new)),
        "mu_cmp_steam": nh(MuFitModel(lp_steam, steam_new)),
    }