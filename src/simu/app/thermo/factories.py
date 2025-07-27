# stdlib
from typing import Iterable, Any
from copy import deepcopy
from collections.abc import Iterator

# external
from yaml import safe_load

# internal
from simu.core.thermo.factory import ThermoFactory
from simu.core.thermo.species import SpeciesDefinition
from simu.core.utilities.types import Map
from simu.app.data import DATA_DIR


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
        from .contributions.augmenters import general
        from simu.core.thermo.contribution import all_contributions
        from simu.core.thermo.state import all_states

        self.register(*all_contributions)
        for state in all_states:
            self.register_state_definition(state)

# TODO: I still use this in some tests, right? Use ThermoStructure instead!

# class ExampleThermoFactory(RegThermoFactory):
#     """This ThermoFactory subclass is capable of creating frames from the base
#     SiMu installation, hence example thermodynamic models that are found in open
#     literature."""
#     def __init__(self):
#         """Default and only constructor"""
#         super().__init__()
#         with open(DATA_DIR / "structures.yml", encoding='UTF-8') as file:
#             self.__structures = safe_load(file)
#
#     @property
#     def structure_names(self) -> Iterable[str]:
#         """The names of all configurations"""
#         return self.__structures.keys()
#
#     def structure(self, name: str) -> Map:
#         """Return a copy of the contribution structure that is registered as
#         ``name``"""
#         return deepcopy(self.__structures[name])
#
#     def create_frame_from_struct(self, species: Map[SpeciesDefinition],
#                                  structure_name: str):
#         """Shortcut for
#         ``factory.create_frame(species, factory.structure(structure_name))``
#         """
#         return super().create_frame(species, self.structure(structure_name))


class ThermoStructure(Map):
    """This class represents the structure of a thermodynamic model,
    representing the list of contributions and holding the state definition.
    """

    __predefined: dict = None

    def __init__(self, state:str, contributions: list):
        self.contributions = contributions
        self.state = state

    def __add__(self, other: str):
        """Add another contribution to the model, calculating further properties
        based on the previously existing."""
        return ThermoStructure(self.state, self.contributions + [other])

    def to_dict(self):

        return {"state": self.state, "contributions": self.contributions}

    @classmethod
    def from_predefined(cls, name: str):
        """Create a pre-defined model structure based on its identifier."""
        cls.__assure_predefined()
        return cls(**cls.__predefined[name])

    @classmethod
    def predefined(cls) -> Iterable[str]:
        """Return a list of pre-defined model structures."""
        cls.__assure_predefined()
        return list(cls.__predefined.keys())

    @classmethod
    def __assure_predefined(cls):
        if cls.__predefined is None:
            with open(DATA_DIR / "structures.yml", encoding='UTF-8') as file:
                cls.__predefined = safe_load(file)

    def __getitem__(self, item) -> Any:
        """Imitate dictionary behaviour, so these objects can directly be used
        in :meth:`~simu.core.thermo.factory.ThermoFactory.create_frame'"""
        if item == "state":
            return self.state
        elif item == "contributions":
            return self.contributions
        else:
            raise KeyError(item)

    def __iter__(self) -> Iterator[str]:
        yield "state"
        yield "contributions"

    def __len__(self) -> int:
        return 2