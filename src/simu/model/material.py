"""Module containing classes to describe materials (thermodynamic phases) in
the modelling context."""

# stdlib
from typing import Optional

from ..thermo.material import MaterialSpec


# internal


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

#     def create(self, name: str, definition: MaterialDefinition):
#         pass
