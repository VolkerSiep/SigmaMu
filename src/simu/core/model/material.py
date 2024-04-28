"""Module containing classes to describe materials (thermodynamic phases) in
the modelling context."""

# stdlib modules
from typing import Optional
from collections.abc import Iterator, Collection

# internal modules
from ..thermo.material import MaterialSpec, Material, MaterialDefinition
from ..utilities.types import Map, MutMap
from ..utilities.errors import DataFlowError


class MaterialHandler(Map[Material]):
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

    def __getitem__(self, name: str) -> Material:
        """Re-obtain the material specification, avoiding the need to keep a
        holding variable in the client scope code."""
        return self.__materials[name]

    def __len__(self) -> int:
        return len(self.__materials)

    def __iter__(self) -> Iterator[str]:
        return iter(self.__materials)

    def create_flow(self, name: str,
                    definition: MaterialDefinition) -> Material:
        """Create a Material as a flow in the local context."""
        return self.__register(name, definition.create_flow())

    def create_state(self, name: str,
                     definition: MaterialDefinition) -> Material:
        """Create a Material as a state in the local context."""
        return self.__register(name, definition.create_state())

    @property
    def ports(self) -> Map[MaterialSpec]:
        """All defined ports of the model"""
        return self.__ports

    def create_proxy(self) -> "MaterialProxy":
        """Create a proxy object for configuration in material context"""
        return MaterialProxy(self)

    def finalise(self, connections: Map[Material]):
        """Internal method called to provide connected ports."""
        self.__materials |= connections

    def __register(self, name: str, material: Material) -> Material:
        """Create a Material in the local context.
        """
        if name in self.__ports:
            raise KeyError(f"{name} is already defined as a port")
        if name in self.__materials:
            raise KeyError(f"{name} is already defined as a material")
        self.__materials[name] = material
        return material


class MaterialProxy(Map[MaterialSpec]):
    handler: MaterialHandler

    def __init__(self, handler: MaterialHandler):
        self.handler = handler
        self.__ports: MutMap[MaterialSpec] = dict(handler.ports)
        self.__connected: MutMap[Material] = {}

    def connect(self, name: str, material: Material):
        """Connect the material in local context to the port with given name
        of the child model"""
        if name not in self:
            raise KeyError(f"Port of name {name} is not defined")
        try:
            spec = self.__ports.pop(name)
        except KeyError:
            raise KeyError(f"Port of name {name} is already connected")
        if not spec.is_compatible(material):
            raise ValueError(f"Provided material on port {name} is "
                             "incompatible to the provided material object")
        self.__connected[name] = material

    def free_ports(self) -> Collection[str]:
        """Return collection of all ports that are yet free"""
        return self.__ports.keys()

    def __getitem__(self, name: str) -> MaterialSpec:
        """Re-obtain the material object, avoiding the need to keep a
        holding variable in the client scope code."""
        return self.handler.ports[name]

    def __len__(self) -> int:
        return len(self.handler.ports)

    def __iter__(self) -> Iterator[str]:
        return iter(self.handler.ports)

    def finalise(self):
        """check that all ports are connected"""
        if self.__ports:
            missing = ", ".join(self.__ports.keys())
            msg = f"The following ports are not connected: {missing}"
            raise DataFlowError(msg)
        self.handler.finalise(self.__connected)
