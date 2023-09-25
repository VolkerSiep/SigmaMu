from typing import Type, Collection, Mapping, Iterable
from pathlib import Path
from yaml import safe_load

from .state import StateDefinition, all_states
from .frame import ThermoFrame
from .contribution import ThermoContribution
from .contributions import all_contributions
from ..data import DATA_DIR


class ThermoFactory:
    """The ``ThermoFactory`` class hosts the definitions for the *model
    contributions*, enabling it to create instances of thermodynamic models of
    class :class:`ThermoFrame`.

    The class is largely meant to be a singleton, but to keep doors open,
    static attributes are avoided."""

    def __init__(self):
        """Parameter-less constructor, initialising the data structure
        to host contribution definitions"""
        self.__contributions = {}
        self.__state_definitions = {}

    def register_state_definition(self, definition_cls: Type[StateDefinition]):
        """Register a new state definition with the name of its class.

        :param definition_cls: The state definition to register
        """
        name = definition_cls.__name__
        if name in self.__state_definitions:
            raise ValueError(f"State definition '{name}' already defined.")
        self.__state_definitions[name] = definition_cls

    def register(self, *contributions: Type[ThermoContribution]):
        """Registers contributions under the name of the class.
        The contributions must be a concrete subclass of
        :class:`ThermoContribution`.

        :param contributions: The contributions to register, being classes
          (not instances)
        """
        for class_ in contributions:
            name = class_.__name__
            if name in self.__contributions:
                raise ValueError(f"Contribution '{name}' already defined.")
            self.__contributions[name] = class_

    @property
    def contribution_names(self) -> Collection[str]:
        """This property contains the full names of all so long registered
        contributions"""
        return set(self.__contributions.keys())

    def create_frame(self, configuration: Mapping) -> ThermoFrame:
        """This factory method creates a :class:`ThermoFrame` object from the
        given ``configuration``.

        :param configuration: A nested dictionary with the
          following root entries:

            - ``species``: A list of strings, representing the names of the
              chemical species. These names are used as identifiers in later
              use of the model, and they define the names of thermodynamic
              parameters required by the model.
            - ``state``: An identifier that represents the type of state as
              defined by :meth:`register_state_definition`.
            - ``contributions``: A list of strings, representing the names
              of the contributions to stack. These identifiers must have been
              defined upfront by calls to :meth:`register_contribution`.
              A contribution entry can also be a dictionary with the following
              keys:

                - ``cls``: The contribution class (required). If this is the
                  only key defined, its effect is as if the sole string of the
                  class name was provided.
                - ``name``: The name of the contribution in the created
                  :class:`ThermoFrame`` object. This name must be unique within
                  the frame definition. It will also be used to define the
                  contribution parameter structure.
                  If skipped, the name will be the same as ``cls``. To be
                  unique, this doesn't work if one contribution class is used
                  multiple times.
                - ``options``: Any data structure that is accepted (and
                  hopefully documented) for the particular contribution.
                  If skipped, an empty dictionary is used.

        :return: The thermodynamic model object
        """
        contributions = {}
        for item in configuration["contributions"]:
            if isinstance(item, dict):
                class_ = item["cls"]
                name = item.get("name", class_)
                options = item.get("options", {})
            else:
                name, class_, options = item, item, {}
            class_ = self.__contributions[class_]
            if name in contributions:
                raise ValueError(f"Duplicate contribution name '{name}'")
            contributions[name] = class_, options

        species = configuration["species"]
        state_def_cls = self.__state_definitions[configuration["state"]]
        result = ThermoFrame(species, state_def_cls(), contributions)

        # set default state
        default = configuration.get("default_state", None)
        if default is not None:
            # make sure the values are float
            default = [
                float(default[0]),
                float(default[1]),
                list(map(float, default[2]))
            ]
        result.default = default
        return result


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

    @property
    def configuration_names(self) -> Iterable[str]:
        """The names of all configurations"""
        return self.__configurations.keys()

    def create_frame(self, configuration: str):
        cfg = self.__configurations[configuration]
        return super().create_frame(cfg | self.__structures[cfg["structure"]])
