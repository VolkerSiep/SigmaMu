from ..utilities import SymbolQuantity


class PropertyHandler:

    def __init__(self):
        self.__defined = {}
        self.__symbols = {}

    def provide(self, name: str, unit: str):
        if name in self.__defined:
            raise ValueError(f"Property '{name}' already defined")
        self.__defined[name] = unit

    def __getitem__(self, name: str) -> SymbolQuantity:
        return self.__symbols[name]

    def __setitem__(self, name: str, symbol: SymbolQuantity):
        unit = self.__defined[name]
        symbol.to(unit)  # raises DimensionalityError if not compatible
        self.__symbols[name] = symbol

    def check(self):
        defined = set(self.__defined.keys())
        provided = set(self.__symbols.keys())

        # are all promised properties calculated?
        missing = defined - provided
        if missing:
            missing = ", ".join(missing)
            raise ValueError(f"Promised properties not provided: {missing}")

        # are
        surplus = provided - defined
        if surplus:
            surplus = ", ".join(surplus)
            raise Warning(f"Surplus properties are provided: {surplus}")

    @property
    def symbols(self):
        return self.__symbols

    @symbols.setter
    def symbols(self, symbols):
        self.__symbols = symbols