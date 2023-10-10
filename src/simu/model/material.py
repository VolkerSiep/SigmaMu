"""Module containing classes to describe materials (thermodynamic phases) in
the modelling context."""

# stdlib modules
from typing import Optional
from collections.abc import Iterator

# internal modules
from ..utilities.types import Map, MutMap
from ..thermo.material import MaterialSpec, Material


class MaterialHandler(Map[MaterialSpec]):
    """The MaterialHandler is the instance to define material ports of a model.
    """
    def __init__(self):
        self.__required = {}

    def port(self, name: str, spec: Optional[MaterialSpec] = None):
        """Define a material port of the given name and specification.
        The name must be unique in this context. If no ``spec`` is given,
        any material is accepted.
        """
        if name in self.__required:
            raise KeyError(f"Material port '{name}' already defined")
        self.__required[name] = MaterialSpec() if spec is None else spec

    def create_proxy(self) -> "MaterialProxy":
        """Create a proxy object for configuration in material context"""
        return MaterialProxy(self)

    def pop_definition(self, name: str) -> MaterialSpec:
        return self.__required.pop(name)

    def __getitem__(self, name: str):
        """Re-obtain the material specification, avoiding the need to keep a
        holding variable in the client scope code."""
        return self.__required[name]

    def __len__(self) -> int:
        return len(self.__required)

    def __iter__(self) -> Iterator[str]:
        return iter(self.__required)


class MaterialProxy(Map[Material]):
    def __init__(self, handler: MaterialHandler):
        self.handler = handler
        self.__materials: MutMap[Material] = {}

    def connect(self, name: str, material: Material):
        try:
            spec = self.handler.pop_definition(name)
        except KeyError:
            raise KeyError(f"Port of name {name} is not defined")
        if not spec.is_compatible(material):
            raise ValueError(f"Provided material on port {name} is "
                             "incompatible to the provided material object")
        self.__materials[name] = material

    def create(self, name: str, material: Material):
        if name in self.handler:
            raise KeyError(f"{name} is already defined as a port")
        if name in self:
            raise KeyError(f"{name} is already defined as a material")
        self.__materials[name] = material

    def __getitem__(self, name: str):
        """Re-obtain the material object, avoiding the need to keep a
        holding variable in the client scope code."""
        return self.__materials[name]

    def __len__(self) -> int:
        return len(self.__materials)

    def __iter__(self) -> Iterator[str]:
        return iter(self.__materials)
