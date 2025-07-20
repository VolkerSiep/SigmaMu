# stdlib
from typing import Iterable

# external
from yaml import safe_load

# internal
from simu.app.data import DATA_DIR
from simu.core.thermo.factory import ThermoFactory


class RegThermoFactory(ThermoFactory):
    """This factory class already registers the current content of
    ``simu.all_contributions`` and ``simu.all_states``, making any
    further registration redundant as long as the definition of the contribution
    and states is loaded, and the entities are decorated with the appropriate
    ``register`` decorator (:class:`simu.registered_contribution` and
    :class:`simu.registered_state`).

    All states and contributions of the ``simu.app`` submodule are automatically
    imported here.
    """
    def __init__(self):
        super().__init__()
        # following imports are needed to trigger registration, even if they
        # are not used locally.
        from .state import HelmholtzState, GibbsState
        from .contributions import basic, special
        from .contributions.iapws import standard, residual
        from .contributions.cubic import core, rk
        from simu.core.thermo.contribution import all_contributions
        from simu.core.thermo.state import all_states

        self.register(*all_contributions)
        for state in all_states:
            self.register_state_definition(state)


class ExampleThermoFactory(RegThermoFactory):
    """This ThermoFactory subclass is capable of creating frames from the base
    SiMu installation, hence thermodynamic models that are found in open
    literature."""
    def __init__(self):
        """Default and only constructor"""
        super().__init__()
        from simu.core.thermo.contribution import all_contributions
        from simu.core.thermo.state import all_states

        with open(DATA_DIR / "structures.yml", encoding='UTF-8') as file:
            self.__structures = safe_load(file)

    @property
    def structure_names(self) -> Iterable[str]:
        """The names of all configurations"""
        return self.__structures.keys()

    def create_frame(self, species, structure: str):
        return super().create_frame(species, self.__structures[structure])
