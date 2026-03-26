# stdlib
from typing import Type, Union, Any
from collections.abc import Mapping, Collection, Sequence

# external
from pydantic import BaseModel, field_validator, ValidationInfo

# internal
from .state import StateDefinition
from .species import SpeciesDefinition
from .frame import ThermoFrame
from .contribution import ThermoContribution

class ContributionConfiguration(BaseModel):
    cls: str
    name: str
    options: Any

class FrameConfiguration(BaseModel):
    state: str
    contributions: Sequence[Union[str, ContributionConfiguration]]

    @field_validator("state")
    @classmethod
    def check_valid(cls, value: str, info: ValidationInfo) -> str:
        states = info.context("valid_states")
        if value not in states:
            raise ValueError(f"State {value} not registered in ThermoFactory")

    @field_validator("contributions", mode="before")
    @classmethod
    def normalize_contributions(cls, contributions: Sequence[Mapping]) \
        -> Sequence[Mapping]:
        def normalize_entry(entry: Union[str, Mapping]) -> Mapping:
            if isinstance(entry, str):
                return {
                    "cls": entry,
                    "name": entry,
                    "options": None
                }
            return entry

        return [normalize_entry(c) for c in contributions]

    @field_validator("contributions")
    @classmethod
    def convert_to_contributions(cls, contributions):
        return [ContributionConfiguration.model_validate(c)
                for c in contributions]


class ThermoFactory:
    """The ``ThermoFactory`` class hosts the definitions for the *model
    contributions*, enabling it to create instances of thermodynamic models of
    class :class:`ThermoFrame`.

    The class is largely meant to be a singleton, but to keep doors open,
    static attributes are avoided."""

    def __init__(self):
        """Parameter-less constructor, initializing the data structure
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
        """Registers contributions under the names of their classes.

        The contributions must be concrete subclasses of
        :class:`ThermoContribution`, and their names must be unique.

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
        """This property contains the full names of all registered
        contributions"""
        return set(self.__contributions.keys())

    def create_frame(self, species: Mapping[str, SpeciesDefinition],
                     configuration: Mapping) -> ThermoFrame:
        """This factory method creates a :class:`ThermoFrame` object from the
        given ``configuration``, and is the recommended way to create
        :class:`ThermoFrame` objects.

        :param species: A dictionary mapping names to species definitions
        :param configuration: A nested dictionary with the following root
          entries:

            - ``state``: A string identifier that represents the type of state
              as defined by :meth:`register_state_definition`.
            - ``contributions``: A list of strings, representing the names
              of the contributions to stack. These identifiers must have been
              defined upfront by calls to :meth:`register`.
              A contribution entry can also be a dictionary with the following
              keys:

                - ``cls``: The contribution class (required). If this is the
                  only key defined, its effect is as if the sole string of the
                  class name was provided.
                - ``name``: The name of the contribution in the created
                  :class:`ThermoFrame` object. This name must be unique within
                  the frame definition. It will also be used to define the
                  contribution parameter structure.
                  If skipped, the name will be the same as ``cls``. To be
                  unique, this doesn't work if one contribution class is used
                  multiple times.
                - ``options``: Any data structure that is accepted (and
                  hopefully documented) for the particular contribution.
                  If skipped, an empty dictionary is used.

        :return: The thermodynamic model (:class:`ThermoFrame`) object
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
