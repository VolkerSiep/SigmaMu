"""Module containing classes to describe materials (thermodynamic phases) in
the modelling context."""

# stdlib
from typing import Optional, Type
from abc import ABC, abstractmethod
from collections.abc import Iterable, Set, Sequence

# internal
from ..utilities.types import NestedQuantityDict
from ..thermo import ThermoFrame

class Augmentor(ABC):
    """An Augmentor is a specific class to extend the physical properties
    calculated on a material instance."""
    def __init__(self, material_spec: "MaterialSpec"):
        self.material_spec = material_spec

    @abstractmethod
    def define(self, material: "Material"):
        """Method to extend the properties of a material object"""


class MaterialSpec:
    """Representation of a requirement to a material object."""
    def __init__(self, species: Optional[Iterable[str]] = None,
                 augmentors: Optional[Iterable[Type[Augmentor]]] = None):
        self.__locked = not (species is None or "*" in species)
        self.__species = set() if species is None else (set(species) - "*")
        self.__augmentors = set() if augmentors is None else set(augmentors)

    @property
    def species(self) -> Set[str]:
        """The set of species that must be provided."""
        return set(self.__species)

    @property
    def augmentors(self) -> Set[Type[Augmentor]]:
        """The set of required augmentor classes"""
        return set(self.__augmentors)

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
        spe, mspe = self.species, set(material.species)
        aug, maug = self.augmentors, material.augmentors
        locked = self.locked
        return not ((spe - mspe) or (locked and (mspe - spe) or (aug - maug)))


class Material:
    def __init__(self,
                 thermo_frame: ThermoFrame,
                 parameter_set: NestedQuantityDict):
        pass

    @property
    def species(self) -> Sequence[str]:
        pass

    @property
    def augmentors(self) -> Set[Type[Augmentor]]:
        pass
