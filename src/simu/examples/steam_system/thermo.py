from yaml import safe_load

from simu import (
    ThermoFactory, InitialState, MaterialDefinition, ThermoParameterStore,
    SpeciesDefinition, StringDictThermoSource, ThermoContribution, Quantity,
    qsum)
from simu.app import all_contributions, HelmholtzState
from simu.app.data import DATA_DIR
from simu.core.utilities.types import MutMap


class StreamProperties(ThermoContribution):
    def define(self, res: MutMap[Quantity]):
        res["G"] = res["n"].T @ res["mu"]
        st = res["T"] * res["S"]
        pv = res["p"] * res["V"]
        res["H"] = res["G"] + st
        res["A"] = res["G"] - pv
        res["U"] = res["A"] + res["T"] * res["S"]

        res["N"] = qsum(res["n"])
        res["m"] = res["n"] * res["mw"]
        res["M"] = qsum(res["m"])

def _iapws(phase: str):
    assert phase.lower() in ["liquid", "gas"]
    fac = ThermoFactory()  # TODO: move factory outside function!
    fac.register(*all_contributions)
    fac.register(StreamProperties)
    fac.register_state_definition(HelmholtzState)
    contributions = [
        "MolecularWeight", "ReducedStateIAPWS", "StandardStateIAPWS",
        "IdealGasIAPWS", "Residual1IAPWS", "Residual2IAPWS",
        "Residual3IAPWS", "Residual4IAPWS",
        f"{phase.capitalize()}IAPWSIdealMix", "StreamProperties"
    ]

    config = {
        "species": ["H2O"],
        "state": "HelmholtzState",
        "contributions": contributions
    }
    species_def =  {"H2O": SpeciesDefinition("H2O")}
    return fac.create_frame(species_def, config)

def _create_store():
    # load thermodynamic parameter database
    store = ThermoParameterStore()
    with open(DATA_DIR / "parameters_iapws.yml") as file:
        data = safe_load(file)
        parameter_source = StringDictThermoSource(data["data"])
    store.add_source(data["meta"]["source"], parameter_source)
    return store

_store = _create_store()
_initial_state = InitialState.from_si(600.0, 100e5, [1.0])
hp_steam = MaterialDefinition(_iapws("gas"), _initial_state, _store)

_initial_state = InitialState.from_si(600.0, 100e5, [1.0])
hp_condensate = MaterialDefinition(_iapws("liquid"), _initial_state, _store)

_initial_state = InitialState.from_si(350.0, 5e4, [1.0])
lp_steam = MaterialDefinition(_iapws("gas"), _initial_state, _store)

_initial_state = InitialState.from_si(300.0, 5e4, [1.0])
lp_condensate = MaterialDefinition(_iapws("liquid"), _initial_state, _store)
