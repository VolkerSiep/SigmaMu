"""This module provides materials based on publically available models and
parameters"""

# stdlib modules
from pathlib import Path

# external modules
from yaml import safe_load

# internal modules
from ..thermo import ThermoFactory, all_states, all_contributions

DATA_DIR = Path(__file__).resolve().parent / "data"


class ExampleThermoFactory(ThermoFactory):
    """This ThermoFactory subclass is capable of creating frames from the base
    SiMu installation, hence thermodynamic models that are found in open
    literature."""
    def __init__(self):
        """Default and only constructor"""
        ThermoFactory.__init__(self)
        self.register(*all_contributions)
        for state in all_states:
            self.register_state_definition(state)

        with open(DATA_DIR / "structures.yml", encoding='UTF-8') as file:
            self.__structures = safe_load(file)
        with open(DATA_DIR / "definitions.yml", encoding='UTF-8') as file:
            self.__configurations = safe_load(file)

    def create_frame(self, configuration: str):
        cfg = self.__configurations[configuration]
        cfg += self.__structures[cfg["structure"]]
        return super().create_frame(cfg)
