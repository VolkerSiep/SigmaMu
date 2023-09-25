from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Optional, Iterable, Collection

from . import ThermoFrame, InitialState, ThermoParameterStore

from ..utilities import Quantity
from ..utilities.types import MutMap


class MaterialSpec:
    """Representation of a requirement to a material object.

    Objects of this class are used to define the requirements to a material,
    typically in a material port of a Model. It defines which species are
    required and whether additional species are allowed.
    """

    def __init__(self, species: Optional[Iterable[str]] = None):
        #        augmenters: Optional[Iterable[Type[Augmenter]]] = None):
        self.__locked = not (species is None or "*" in species)
        self.__species = set() if species is None else set(species) - set("*")
        # self.__augmenters = set() if augmenters is None else set(augmenters)

    @property
    def species(self) -> set[str]:
        """The set of species that must be provided."""
        return set(self.__species)

    @property
    def locked(self) -> bool:
        """Whether the specifcation allows other species than the ones
        listed as :attr:`species`."""
        return self.__locked

    def is_compatible(self, material: "Material") -> bool:
        """Return true if the given material is compatible with this one.
        That is:
          - none of the specified species are missing in the material
          - if specification locks species set, none of the material
            species are missing in the specification
        """
        spe, mspe = self.species, set(material.species)
        # aug, maug = self.augmenters, material.augmenters
        locked = self.locked
        # return not ((spe - mspe) or (locked and (mspe - spe) or (aug - maug)))
        return not ((spe - mspe) or (locked and (mspe - spe)))

    # @property
    # def augmenters(self) -> set[Type[Augmenter]]:
    #     """The set of required augmentor classes"""
    #     return set(self.__augmenters)


class Material(MutMap[Quantity]):
    """This class represents a material"""
    def __init__(self,
                 definition: "MaterialDefinition",
                 flow: bool):
        self.definition = definition
        self.initial_state = definition.initial_state

        frame = definition.frame
        params = definition.store.get_symbols(frame.parameter_structure)
        state = frame.create_symbol_state()
        props = frame(state, params, squeeze_results=False, flow=flow)
        self.__properties = {n: p for n, p in props.items()
                             if not n.startswith("_")}

        # apply augmenters as required in definition

    @property
    def species(self) -> Collection[str]:
        """The species names"""
        return self.definition.frame.species

    def __getitem__(self, key: str) -> Quantity:
        return self.__properties[key]

    def __setitem__(self, key: str, value: Quantity):
        if key in self.__properties:
            raise KeyError(f"Property '{key}' already exists in material")
        self.__properties[key] = value

    def __delitem__(self, key):
        raise TypeError("Property deletion is disabled to avoid opaque "
                        "property name interpretation!")

    def __iter__(self):
        return iter(self.__properties)

    def __len__(self):
        return len(self.__properties)

    # @property
    # def augmenters(self) -> set[Type[Augmenter]]:
    #     pass


@dataclass
class MaterialDefinition:
    """A ``MaterialDefinition`` object defines a material type by its

      - frame of thermodynamic contributions,
      - initial state
      - source of thermodynamic parameters
    """
    frame: ThermoFrame
    initial_state: InitialState
    store: ThermoParameterStore

    def __post_init__(self):
        num_species = len(self.frame.species)
        num_init_species = len(self.initial_state.mol_vector.magnitude)
        if num_init_species != num_species:
            raise ValueError(
                f"Incompatible initial state with {num_init_species} "
                f"species, while {num_species} is/are expected"
            )

    @property
    def spec(self) -> MaterialSpec:
        """Return a material spec object that is implemented by this
        definition"""
        return MaterialSpec(self.frame.species)

    def create_flow(self) -> Material:
        return Material(self, True)

    def create_state(self) -> Material:
        return Material(self, False)


class Augmenter(ABC):
    """An Augmenter is a specific class to extend the physical properties
    calculated on a material instance."""
    def __init__(self, frame: ThermoFrame):
        self.species = frame.species

    @abstractmethod
    def define(self, material: "Material"):
        """Method to extend the properties of a material object"""
        ...
