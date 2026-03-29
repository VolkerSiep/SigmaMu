# stdlib
from typing import Type, Any
from collections.abc import Collection, Sequence

# external
from pydantic import (
    BaseModel, field_validator, model_validator,
    ValidationInfo, Field, TypeAdapter)

# internal
from simu.core.utilities.types import Map
from .state import StateDefinition
from .species import SpeciesDefinition
from .frame import ThermoFrame
from .contribution import ThermoContribution


class FrameContributionConfiguration(BaseModel):
    """A data structure representing the configuration of a thermodynamic
    contribution instance."""
    cls: str
    """The identifier of the contribution class, as registered via
    :meth:`~simu.ThermoFactory.register`. or the decorator
    :func:`~simu.registered_contribution`.
    """
    name: str = Field(default=None)
    """The name of the contribution, by default equal to the class identifier,
    but to be explicitly defined if one contribution is included multiple times,
    such as a mixing rule applied on several model parts."""
    options: Any = Field(default=None)
    """Some contributions support or even require options. In such cases,
    the options are flexibly defined by the contribution and documented
    individually."""

    @model_validator(mode="after")
    def default_name(self):
        if self.name is None:
            self.name = self.cls
        return self

    @model_validator(mode="before")
    @classmethod
    def from_str_or_dict(cls, value: str|Map):
        """This method allows the class to be created from a string entity,
        which then will represent both the class identifier and contribution
        name.
        """
        return {"cls": value} if isinstance(value, str) else value

    @field_validator("cls", mode="before")
    @classmethod
    def validate_class(cls, cls_: str, info: ValidationInfo):
        contributions = info.context.get("contributions", [])
        if cls_ not in contributions:
            raise ValueError(
                f"Contribution '{cls_}' not registered in ThermoFactory")
        return cls_


FrameContributionList = TypeAdapter(Sequence[FrameContributionConfiguration])


class FrameConfiguration(BaseModel):
    """A data structure representing the configuration of a thermodynamic frame
    object.
    """
    state: str
    """A string identifier that represents the type of state as defined by 
    :meth:`~simu.ThermoFactory.register_state_definition` or via the decorator
    :func:`~simu.registered_state`."""

    contributions: Sequence[FrameContributionConfiguration]
    """Each element of this list represents the configuration of a 
    thermodynamic contribution."""

    @field_validator("state", mode="before")
    @classmethod
    def validate_state(cls, state: str, info: ValidationInfo) -> str:
        states = info.context.get("states", [])
        if state not in states:
            raise ValueError(f"State '{state}' not registered in ThermoFactory")
        return state

    @classmethod
    @field_validator("contributions", mode="before")
    def validate_contributions(cls, contributions: Sequence[Map|str]) \
            -> Sequence[FrameContributionConfiguration]:
        return FrameContributionList.validate_python(contributions)


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

    def create_frame(self, species: Map[SpeciesDefinition],
                     configuration: Map[Any]) -> ThermoFrame:
        """This factory method creates a :class:`ThermoFrame` object from the
        given ``configuration``, and is the recommended way to create
        :class:`ThermoFrame` objects.

        :param species: A dictionary mapping names to species definitions
        :param configuration:
            A nested dictionary, representing an instance
            of :class:`~simu.core.thermo.factory.FrameConfiguration`.

            A valid input (provided the registration of given entities) is

            .. code::

                {
                    "species": ["N2", "O2"],
                    "state": "HelmholtzState",
                    "contributions": [
                        "H0S0ReferenceState", "LinearHeatCapacity",
                        "StandardState", "IdealMix", "HelmholtzIdealGas"
                    ],
                }

        :return: The thermodynamic model (:class:`ThermoFrame`) object
        """
        context = {
            "states": self.__state_definitions.keys(),
            "contributions": self.__contributions.keys()
        }
        config = FrameConfiguration.model_validate(
            configuration, context=context)
        contributions = {}
        for item in config.contributions:
            class_ = self.__contributions[item.cls]
            if item.name in contributions:
                raise ValueError(f"Duplicate contribution name '{item.name}'")
            contributions[item.name] = class_, item.options

        state_def_cls = self.__state_definitions[config.state]
        return ThermoFrame(species, state_def_cls(), contributions)

