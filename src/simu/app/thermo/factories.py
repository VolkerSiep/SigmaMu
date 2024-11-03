# stdlib
from typing import Iterable

# external
from yaml import safe_load

# internal
from . import all_contributions, all_states
from ..data import DATA_DIR
from ... import ThermoFactory


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

    @property
    def structure_names(self) -> Iterable[str]:
        """The names of all configurations"""
        return self.__structures.keys()

    def create_frame(self, species, structure: str):
        return super().create_frame(species, self.__structures[structure])
