"""Module containing classes to describe materials (thermodynamic phases) in
the modelling context."""

# stdlib
from typing import Optional, Type
from abc import ABC, abstractmethod
from collections.abc import Iterable, Collection

# internal
from ..utilities import Quantity
from ..utilities.types import NestedMap
from ..thermo import ThermoFrame


class Augmentor(ABC):
    """An Augmentor is a specific class to extend the physical properties
    calculated on a material instance."""
    def __init__(self, frame: ThermoFrame):
        self.species = frame.species

    @abstractmethod
    def define(self, material: "Material"):
        """Method to extend the properties of a material object"""
        ...


class MaterialSpec:
    """Representation of a requirement to a material object.

    Objects of this class are used to define the requirements to a material,
    typically in a material port of a Model. It defines which species are
    required and whether additional species are allowed.
    """

    def __init__(self, species: Optional[Iterable[str]] = None,
                 augmentors: Optional[Iterable[Type[Augmentor]]] = None):
        self.__locked = not (species is None or "*" in species)
        self.__species = set() if species is None else set(species) - set("*")
        self.__augmentors = set() if augmentors is None else set(augmentors)

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
          - none of the specified augmentors are missing in the material
        """
        # TODO: it's ok that not all augmentors are applied. The spec can do
        #  that on instantiation

        spe, mspe = self.species, set(material.species)
        aug, maug = self.augmentors, material.augmentors
        locked = self.locked
        return not ((spe - mspe) or (locked and (mspe - spe) or (aug - maug)))

    @property
    def augmentors(self) -> set[Type[Augmentor]]:
        """The set of required augmentor classes"""
        return set(self.__augmentors)


class Material:
    def __init__(self, frame: ThermoFrame, flow: bool = False):
        self.__frame = frame  # Frame is the new MaterialDefinition
        # TODO:
        #  also generate a state corresponding to the frame
        #  ask frame for parameters and ask ThermoParameterStore for symbols.



    @property
    def species(self) -> Collection[str]:
        """The species names"""
        return self.__frame.species

    @property
    def augmentors(self) -> set[Type[Augmentor]]:
        pass


class MaterialHandler:
    def __init__(self):
        self.__required = {}
        self.__materials = {}

    def require(self, name: str, spec: Optional[MaterialSpec] = None):
        """Define a material requirement of the given name and specification.
        The name must be unique in this context. If no ``spec`` is given,
        any material is accepted.
        """
        if name in self.__required:
            raise KeyError(f"Material '{name}' already required")
        self.__required[name] = MaterialSpec() if spec is None else spec

    def __getitem__(self, name: str):
        return self.__materials[name]

    def create(self, name: str, definition: MaterialDefinition):
        pass
