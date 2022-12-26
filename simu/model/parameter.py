from ..utilities import Quantity, SymbolQuantity


class ParameterDefinition:

    def __init__(self):
        self.__symbols = {}
        self.__values = {}
        self.__required = set()

    def define(self, name: str, value: float = None, unit: str = "dimless"):
        """Define a parameter in the context of the model. The name must be
        unique in the context. If value is not given (``None``), the parameter
        must be provided from the outside - and be compliant with the given
        unit of measurement. Otherwise, the value-unit pair defines the
        default value of the parameter, if not provided from parent module"""
        if self.defined(name):
            raise ValueError(f"Parameter '{name}' already defined")

        if value is None:
            self.__required.add(name)
        else:
            self.__values[name] = Quantity(value, unit)

        self.__symbols[name] = SymbolQuantity(name, unit)

    @property
    def symbols(self) -> dict[str, SymbolQuantity]:
        """Return dictionary of symbols representing all defined parameters,
        including those that are dependent symbols.
        """
        return self.__symbols

    def __getitem__(self, name: str) -> SymbolQuantity:
        return self.__symbols[name]

    def defined(self, name: str):
        return name in self.__symbols

    def create_handler(self):
        return ParameterHandler(self.symbols, self.__values, self.__required)


class ParameterHandler:
    """When initiated by the :class`ParameterDefinition` object, the instance
    holds the known value quantities for a subset of the parameters. The
    remaining parameters are in the ``required`` set and need to be either
    given a value or a dependent symbol.

    After all configuration has been done, all required parameters shall be
    provided a dependent symbol or a numeric quantity.

    Then, the ``_symbols`` dictionary contains the dependent symbols, and the
    ``_values`` dictionary contains the independent parameters.

    """

    def __init__(self, symbols: dict[str, SymbolQuantity],
                 values: dict[str, Quantity], required: set[str]):
        # I need to take a real copy here, as the model (ParameterDefinition)
        # instance can be reused. The _initial_ parameter situation is however
        # always the same.
        self.__symbols = dict(symbols)
        self.__values = dict(values)
        self.__required = set(required)
        self.__provided = set()
        self.__finalised = False

    def defined(self, name: str):
        return name in self.__symbols

    def update(self, **kwargs: dict[str, str]):
        """Update or provide the default value of the parameter"""
        for name, value in kwargs.items():
            if not self.defined(name):
                raise KeyError(f"Parameter '{name}' undefined")
            if name in self.__provided:
                raise KeyError(f"Parameter '{name}' not independent")
            value = Quantity(value)
            self.__values[name] = value
            self.__required -= set([name])

    def provide(self, **kwargs):
        """Connect input to another symbolic quantity"""
        if self.__finalised:
            raise RuntimeError("Providing parameter to finalised model")
        for name, symbol in kwargs.items():
            if not self.defined(name):
                raise KeyError(f"Parameter '{name}' undefined")
            self.__symbols[name] = symbol
            self.__required -= set([name])
            self.__provided.add(name)

            if name in self.__values:  # parameter no longer independent
                del self.__values[name]

    def finalise(self):
        """Check that the configuration is ready to be finalised.
        Throw appropriate exceptions if this is not the case."""
        if self.__required:
            missing = ", ".join(self.__required)
            raise ValueError(f"Non-supplied required parameters: {missing}")
        self.__finalised = True

    @property
    def values(self):
        return self.__values

    @property
    def symbols(self) -> dict[str, SymbolQuantity]:
        """Return dictionary of symbols representing all defined parameters,
        including those that are dependent symbols.
        """
        return self.__symbols

    @property
    def free_symbols(self) -> dict[str, SymbolQuantity]:
        return {name: self.__symbols[name] for name in self.__values}
