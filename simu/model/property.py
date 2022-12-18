from ..utilities import SymbolQuantity


class PropertyHandler:

    def __init__(self):
        self.__symbols = {}
        self.__raw_symbols = {}
        self.__required = set()

    def provide(self, name: str, unit: str):
        if name in self.__symbols:
            raise ValueError(f"Property '{name}' already defined")
        self.__symbols[name] = SymbolQuantity(name, unit)
        self.__required.add(name)

    def __getitem__(self, name: str) -> SymbolQuantity:
        return self.__symbols[name]

    def __setitem__(self, name: str, symbol: SymbolQuantity):
        if not name in self.__required:
            raise ValueError(f"Property {name} already set or not defined")
        prototype = self.__symbols[name]
        symbol.to(prototype.units)  # DimensionalityError if not compatible
        self.__raw_symbols[name] = symbol
        self.__required.remove(name)

    def check(self):
        if self.__required:
            missing = ", ".join(self.__required)
            raise ValueError(f"Promised properties not provided: {missing}")

    @property
    def symbols(self):
        return self.__symbols

    @property
    def raw_symbols(self):
        """The raw symbols as defined in :meth:`simu.Model.define`"""
        return self.__raw_symbols

    def finalise(self, symbols):
        self.__symbols = symbols
