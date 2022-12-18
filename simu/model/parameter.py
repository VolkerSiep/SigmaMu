from ..utilities import Quantity, SymbolQuantity


class ParameterHandler:

    def __init__(self):
        self.__symbols = {}
        self.__values = {}
        self.__required = set()
        self.__finalised = False

    def define(self, name: str, value: float = None, unit: str = "dimless"):
        """Define a parameter in the context of the model. The name must be
        unique in the context. If value is not given (``None``), the parameter
        must be provided from the outside - and be compliant with the given
        unit of measurement. Otherwise, the value-unit pair defines the
        default value of the parameter, if not provided from parent module"""
        if name in self.__symbols:
            raise ValueError(f"Parameter '{name}' already defined")

        if value is None:
            self.__required.add(name)
        else:
            self.__values[name] = Quantity(value, unit)

        self.__symbols[name] = SymbolQuantity(name, unit)

    def update(self, **kwargs: dict[str, str]):
        """Update the default value of the parameter"""
        for name, value in kwargs.items():
            if name not in self.__symbols:
                raise KeyError(f"Parameter '{name}' undefined")
            self.__values[name] = Quantity(value)
            self.__required -= set([name])

    def provide(self, **kwargs):
        """Connect input to another symbolic quantity"""
        if self.__finalised:
            raise RuntimeError("Providing parameter to finalised model")
        for name, symbol in kwargs.items():
            if name not in self.__symbols:
                raise KeyError(f"Parameter '{name}' undefined")
            self.__symbols[name] = symbol
            self.__required -= set([name])

            # disregard, value, it's no longer relevant.
            try:
                del self.__values[name]
            except KeyError:
                pass

    def check(self):
        if self.__required:
            param = ", ".join(self.__required)
            raise ValueError(f"Non-supplied required parameters: {param}")

    def __getitem__(self, name: str) -> SymbolQuantity:
        return self.__symbols[name]

    @property
    def symbols(self):
        return self.__symbols

    @property
    def values(self):
        return self.__values

    def finalise(self):
        self.__finalised = True
