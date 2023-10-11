"""Module containing classes to describe materials (thermodynamic phases) in
the modelling context."""

# stdlib modules
from typing import Optional
from collections.abc import Iterator

# internal modules
from ..thermo.material import MaterialSpec, Material
from ..utilities.types import Map, MutMap
from ..utilities.errors import DataFlowError


class MaterialHandler(Map[MaterialSpec]):
    """The material handler maintains the thermodynamic states represented as
    flows and states. When a model is created, the ``interface`` method can be
    used to define material ports. In the ``with`` context, invoked by the
    parent model, defined ports are to be connected to ``Material`` instances
    in the parent context.
    Finally, the model can create further (local) material instances.
    """
    def __init__(self):
        self.__materials: MutMap[Material] = {}
        self.__ports: MutMap[MaterialSpec] = {}

    def define_port(self, name: str, spec: Optional[MaterialSpec] = None):
        """Define a material port of the given name and specification.
        The name must be unique in this context. If no ``spec`` is given,
        any material is accepted.
        """
        if name in self.__ports:
            raise KeyError(f"Material port '{name}' already defined")
        self.__ports[name] = MaterialSpec() if spec is None else spec

    def __getitem__(self, name: str):
        """Re-obtain the material specification, avoiding the need to keep a
        holding variable in the client scope code."""
        return self.__materials[name]

    def __len__(self) -> int:
        return len(self.__materials)

    def __iter__(self) -> Iterator[str]:
        return iter(self.__materials)

    def create(self, name: str, material: Material, *, _as_port=False):
        if not _as_port and name in self.__ports:
            raise KeyError(f"{name} is already defined as a port")
        if name in self.__materials:
            raise KeyError(f"{name} is already defined as a material")
        self.__materials[name] = material

    @property
    def ports(self) -> Map[MaterialSpec]:
        return self.__ports

    def create_proxy(self) -> "MaterialProxy":
        """Create a proxy object for configuration in material context"""
        return MaterialProxy(self)


class MaterialProxy(Map[Material]):
    def __init__(self, handler: MaterialHandler):
        self.handler = handler
        self.__ports = dict(handler.ports)

    def connect(self, name: str, material: Material):
        if name not in self:
            raise KeyError(f"Port of name {name} is not defined")
        try:
            spec = self.__ports.pop(name)
        except KeyError:
            raise KeyError(f"Port of name {name} is already connected")
        if not spec.is_compatible(material):
            raise ValueError(f"Provided material on port {name} is "
                             "incompatible to the provided material object")
        self.handler.create(name, material, _as_port=True)

    def free_ports(self) -> Iterator[str]:
        return iter(self.__ports)

    def __getitem__(self, name: str):
        """Re-obtain the material object, avoiding the need to keep a
        holding variable in the client scope code."""
        return self.handler.ports[name]

    def __len__(self) -> int:
        return len(self.handler.ports)

    def __iter__(self) -> Iterator[str]:
        return iter(self.handler.ports)

    def finalise(self):
        # check that all ports are connected
        if self.__ports:
            missing = ", ".join(self.__ports.keys())
            msg = f"The following ports are not connected: {missing}"
            raise DataFlowError(msg)
