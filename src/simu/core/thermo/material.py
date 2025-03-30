from abc import ABC, abstractmethod
from typing import Optional
from collections.abc import Iterable, Collection, Mapping, Sequence

from ..utilities import (
    Quantity, SymbolQuantity, QuantityDict, QFunction, extract_sub_structure)
from ..utilities.types import Map, MutMap, NestedMap

from . import (
    ThermoFrame, InitialState, ThermoParameterStore, SpeciesDB,
    ThermoFactory, SpeciesDefinition)


class MaterialSpec:
    """Representation of a requirement to a material object.

    Objects of this class are used to define the requirements to a material,
    typically in a material port of a Model. It defines which species are
    required and whether additional species are allowed.
    """

    def __init__(self, species: Optional[Iterable[str]] = None,
                 flow: bool = True):
        """Create an instance based on the arguments as follows:

        :species: If ``None`` (default), allow any species. If containing
          a species name ``*``, allow any species, but demand the species
          that are element of the given iterable. If  ``species`` does not
          contain ``*``, lock the specification to the given set of species.
        :flow: ``True`` if the material needs to represent a flow, ``False`` if
          a state is to be represented.
        """
        self.__flow = flow
        self.__locked = not (species is None or "*" in species)
        self.__species = set() if species is None else set(species) - set("*")

    @classmethod
    @property
    def flow(cls):
        """A generic flow instance, allowing any species"""
        return cls()

    @classmethod
    @property
    def state(cls):
        """A generic state instance, allowing any species"""
        return cls(flow=False)

    @property
    def species(self) -> set[str]:
        """The set of species that must be provided"""
        return set(self.__species)

    @property
    def locked(self) -> bool:
        """Whether the specification allows other species than the ones
        listed as :attr:`species`."""
        return self.__locked

    def is_flow(self) -> bool:
        """Whether the specification is a flow (``True``) or
        a state (``False``)"""
        return self.__flow

    def is_compatible(self, material: "Material") -> bool:
        """Return true if the given material is compatible with this one.
        That is:

        - none of the specified species are missing in the material
        - if specification locks species set, none of the material
          species are missing in the specification
        - flows are only compatible with flows, states are only compatible
          with states.
        """
        spe, mspe = self.species, set(material.species)
        locked = self.locked
        flow_comp = self.is_flow() == material.is_flow()
        return flow_comp and not ((spe - mspe) or (locked and (mspe - spe)))


class Material(MutMap[Quantity]):
    """This class represents a material"""

    definition: "MaterialDefinition"
    """The underlying definition"""

    initial_state: InitialState
    """The provided initial state"""

    def __init__(self,
                 definition: "MaterialDefinition",
                 flow: bool):
        """Define a material based on the given definition. The ``flow`` flag
        indicates whether a *flow* or a *state* is to be generated."""
        self.definition = definition
        self.initial_state = definition.initial_state

        frame = definition.frame
        params = definition.store.get_symbols(frame.parameter_structure)
        self.__state = frame.create_symbol_state()
        props = frame(self.__state, params, squeeze_results=False, flow=flow)
        vectors = frame.vector_keys

        def convert(n: str, prop: SymbolQuantity) -> Quantity | QuantityDict:
            """If property is a registered vector, convert it to a QuantityDict,
            otherwise just return the quantity itself."""
            nonlocal vectors
            mag, unit = prop.magnitude, prop.units
            try:
                keys = vectors[n]
            except KeyError:
                return prop
            if mag.size() != (len(keys), 1):
                msg = f"Property {n} has improper shape: {mag.size()}, " \
                      f"should be ({len(keys)}, 1)"
                raise ValueError(msg)
            result = {s: Quantity(mag[i], unit) for i, s in enumerate(keys)}
            return QuantityDict(result)

        self.__properties = {n: convert(n, p) for n, p in props["props"].items()
                             if not n.startswith("_")}
        self.__bounds = props.get("bounds", {})
        self.__residuals = props.get("residuals", {})
        self.__normed_residuals = props.get("normed_residuals", {})
        self.__flow = flow

        # create a QFunction to map a state and parameters into a new initial
        # state
        args = {"state": Quantity(self.__state), "param": params}
        props = frame(self.__state, params, squeeze_results=False, flow=flow)
        res = {n: props["props"][n] for n in "Tpn"}
        self.__ini_func = QFunction(args, res, "ini_func")

    def retain_initial_state(self, state: Sequence[float],
                             parameters: NestedMap[Quantity]):
        param_struct = self.definition.frame.parameter_structure
        store_name = self.definition.store.name
        parameters = extract_sub_structure(parameters[store_name], param_struct)
        args = {"state": Quantity(state), "param": parameters}
        res = self.__ini_func(args)
        self.initial_state = InitialState(
            temperature=res["T"],
            pressure=res["p"],
            mol_vector=res["n"].reshape(-1))  # make sure it's a vector

    @property
    def species(self) -> Collection[str]:
        """The species names"""
        return self.definition.species

    @property
    def species_definitions(self) -> Map[SpeciesDefinition]:
        return self.definition.species_definitions

    @property
    def sym_state(self) -> Map[Quantity]:
        """A dictionary of quantities representing the state"""
        state = self.__state.nonzeros()
        return {f"x_{k:03d}": Quantity(x_k) for k, x_k in enumerate(state)}

    @property
    def bounds(self) -> NestedMap[Quantity]:
        """Return the bounds of the material as a nested mapping. The top
        level keys will be the names of the original contributions, while the
        second level keys will be the names of the bounds. The mapping will
        only have these two levels.
        """
        return self.__bounds

    def residuals(self, normed: bool = False) -> NestedMap[Quantity]:
        """Return the residuals of the material as a nested mapping. The top
        level keys will be the names of the original contributions, while the
        second level keys will be the names of the residuals. The mapping will
        only have these two levels.

        If ``normed`` is ``true``, return the dimensionless ratio of residual
        values and their tolerances instead.
        """
        return self.__normed_residuals if normed else self.__residuals


    def __getitem__(self, key: str) -> Quantity:
        """Return a (symbolic) property that is calculated by the underlying
        thermodynamic model or later supplements via :meth:`__setitem__`."""
        return self.__properties[key]

    def __setitem__(self, key: str, symbol: Quantity):
        """Supplement the material with another property. The ``symbol``
        argument must be a symbolic quantity, and it is highly immoral to supply
        a symbol that is a function of more than prior material properties.
        """
        if key in self.__properties:
            raise KeyError(f"Property '{key}' already exists in material")
        self.__properties[key] = symbol

    def __delitem__(self, key):
        raise TypeError("Property deletion is disabled to avoid opaque "
                        "property name interpretation.")

    def __iter__(self):
        return iter(self.__properties)

    def __len__(self):
        return len(self.__properties)

    def is_flow(self):
        """Return if the material represents a flow (``True``) or not
        (``False``)"""
        return self.__flow


class MaterialDefinition:
    """A ``MaterialDefinition`` object defines a material type by its

      - frame of thermodynamic contributions,
      - initial state
      - source of thermodynamic parameters
    """

    __frame: ThermoFrame
    __initial_state: InitialState
    __store: ThermoParameterStore

    def __init__(self, frame: ThermoFrame, initial_state: InitialState,
                 store: ThermoParameterStore):
        self.__frame = frame
        self.initial_state = initial_state
        store.get_symbols(frame.parameter_structure)  # trigger querying params
        self.__store = store

    @property
    def spec(self) -> MaterialSpec:
        """Return a material spec object that is implemented by this
        definition"""
        return MaterialSpec(self.frame.species)

    @property
    def frame(self) -> ThermoFrame:
        return self.__frame

    @property
    def store(self) -> ThermoParameterStore:
        return self.__store

    @property
    def species(self) -> Collection[str]:
        """The species names"""
        return self.frame.species

    @property
    def species_definitions(self) -> Map[SpeciesDefinition]:
        return self.frame.species_definitions

    @property
    def initial_state(self):
        return self.__initial_state

    @initial_state.setter
    def initial_state(self, new_state: InitialState):
        num_species = len(self.__frame.species)
        num_init_species = len(new_state.mol_vector.magnitude)
        if num_init_species != num_species:
            raise ValueError(
                f"Incompatible initial state with {num_init_species} "
                f"species, while {num_species} is/are expected"
            )
        self.__initial_state = new_state

    def create_flow(self) -> Material:
        return Material(self, True)

    def create_state(self) -> Material:
        return Material(self, False)


class MaterialLab:
    """A MaterialLab is an object to help design and define material
    definitions. It aggregates one or more ThermoParameterStore, the global
    SpeciesDefinition, and the ThermoFactory.

    A new material can be defined based on a species set, an initial state,
    and the contribution list.
    """
    def __init__(self, factory: ThermoFactory, species_db: SpeciesDB,
                 param_store: ThermoParameterStore):
        self.__factory = factory
        self.__species_db = species_db
        self.__param_store = param_store

    def define_material(self, species: Sequence[str],
                        initial_state: InitialState,
                        structure: Mapping) -> MaterialDefinition:
        """This method creates a material definition."""
        species_map = {s: self.__species_db[s] for s in species}
        frame = self.__factory.create_frame(species_map, structure)
        return MaterialDefinition(frame, initial_state, self.__param_store)
