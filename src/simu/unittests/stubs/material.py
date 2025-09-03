from collections.abc import Collection

from simu import Quantity, QuantityDict, SymbolQuantity
from simu.core.utilities.types import MutMap


class MaterialStub(MutMap[Quantity | QuantityDict]):
    current_instance = 0

    def __init__(self, definition: "MaterialDefinitionStub", flow: bool):
        self.__flow = flow
        self.definition = definition
        MaterialStub.current_instance += 1
        i = MaterialStub.current_instance

        def f(unit):
            return f"{unit}/s" if flow else unit

        def species_dict(name, unit, extensive=False):
            unit = f(unit) if extensive else unit
            return QuantityDict({s: SymbolQuantity(f"{name}_{s}_{i}", unit)
                                 for s in self.species})

        self.__properties = {
            "T": SymbolQuantity(f"T_{i}", "K"),
            "p": SymbolQuantity(f"p_{i}", "Pa"),
            "S": SymbolQuantity(f"S_{i}", f("J/K")),
            "V": SymbolQuantity(f"V_{i}", f("m^3")),
            "n": species_dict("n", "mol", extensive=True),
            "mu": species_dict("mu", "J/K"),
            "mw": species_dict("mw", "kg/mol")
        }

    @property
    def species(self) -> Collection[str]:
        """The species names"""
        return self.definition.species

    @property
    def species_definitions(self):
        raise NotImplemented

    @property
    def sym_state(self):
        raise NotImplemented

    @property
    def bounds(self):
        return {}

    def residuals(self, normed: bool = False):
        return {}

    def __getitem__(self, key: str) -> Quantity:
        return self.__properties[key]

    def __setitem__(self, key: str, symbol: Quantity):
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


class MaterialDefinitionStub:
    def __init__(self, species: Collection[str]):
        self.__species = species
        self.__store = ThermoParameterStoreStub

    @property
    def spec(self):
        raise NotImplemented

    @property
    def frame(self):
        raise NotImplemented

    @property
    def store(self):
        return self.__store

    @property
    def species(self) -> Collection[str]:
        return self.__species

    @property
    def species_definitions(self):
        raise NotImplemented

    @property
    def initial_state(self):
        raise NotImplemented

    @initial_state.setter
    def initial_state(self, new_state):
        pass

    def create_flow(self) -> MaterialStub:
        return MaterialStub(self, True)

    def create_state(self) -> MaterialStub:
        return MaterialStub(self, False)


class ThermoParameterStoreStub:
    pass
