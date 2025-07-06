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

class MaterialFactory:
    def __init__(self):
        self.factory = ThermoFactory()
        self.factory.register(*all_contributions)
        self.factory.register(StreamProperties)
        self.factory.register_state_definition(HelmholtzState)

        self.store = ThermoParameterStore()
        with open(DATA_DIR / "parameters_iapws.yml") as file:
            data = safe_load(file)
            parameter_source = StringDictThermoSource(data["data"])
        self.store.add_source(data["meta"]["source"], parameter_source)

        self.frames = {p: self._create_iapws(p) for p in ["gas", "liquid"]}

    def create(self, phase:str, t: float, p: float, n: float):
        assert phase.lower() in ["liquid", "gas"]
        initial_state = InitialState.from_si(t, p, [n])
        return MaterialDefinition(self.frames[phase], initial_state, self.store)

    def _create_iapws(self, phase:str):
        assert phase.lower() in ["liquid", "gas"]
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
        return self.factory.create_frame(species_def, config)

factory = MaterialFactory()

hp_steam = factory.create("gas", 600, 100e5, 100)
hp_condensate = factory.create("liquid", 600, 100e5, 100)
lp_steam = factory.create("gas", 350, 5e4, 100)
lp_condensate = factory.create("liquid", 300, 5e4, 100)
