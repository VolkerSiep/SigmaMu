"""Module related to parameter handling during process modelling"""

from ..utilities import Quantity, SymbolQuantity
from ..utilities.quantity import QuantityDict, SymbolQuantityDict
from .utils import ModelStatus


class ParameterHandler:
    """When initiated by the :class:`ParameterDefinition` object, the instance
    holds the known value quantities for a subset of the parameters. The
    remaining parameters are still marked as required and need to be either
    given a value or a dependent symbol.

    After all configuration has been done, all required parameters shall be
    provided a dependent symbol or a numeric quantity.
    """

    def __init__(self, symbols: SymbolQuantityDict, values: QuantityDict,
                 required: set[str]):
        # I need to take a real copy here, as the model (ParameterDefinition)
        # instance can be reused. The _initial_ parameter situation is however
        # always the same.
        self.__symbols: SymbolQuantityDict = dict(symbols)
        self.__values: QuantityDict = dict(values)
        self.__required: set[str] = set(required)
        self.__provided: set[str] = set()
        self.status = ModelStatus.READY

    def defined(self, name: str) -> bool:
        """Check whether a parameter of given ``name`` is defined and return
        the result as a boolean."""
        return name in self.__symbols

    def update(self, **kwargs: str) -> None:
        """Update or provide the default value of the parameter"""
        self.status.assure(ModelStatus.READY, "parameter update")
        for name, value in kwargs.items():
            if not self.defined(name):
                raise KeyError(f"Parameter '{name}' undefined")
            if name in self.__provided:
                raise KeyError(f"Parameter '{name}' not independent")
            value = Quantity(value)
            self.__values[name] = value
            self.__required -= set([name])

    def provide(self, **kwargs: SymbolQuantity) -> None:
        """Connect input to another symbolic quantity"""
        self.status.assure(ModelStatus.READY, "parameter provide")
        for name, symbol in kwargs.items():
            if not self.defined(name):
                raise KeyError(f"Parameter '{name}' undefined")
            self.__symbols[name] = symbol
            self.__required -= set([name])
            self.__provided.add(name)

            if name in self.__values:  # parameter no longer independent
                del self.__values[name]

    def finalise(self) -> None:
        """Check that the configuration is ready to be finalised.
        Throw appropriate exceptions if this is not the case."""
        self.status.assure(ModelStatus.READY, "finalising")
        if self.__required:
            missing = ", ".join(self.__required)
            raise ValueError(f"Non-supplied required parameters: {missing}")
        self.status = ModelStatus.FINALISED

    @property
    def values(self) -> QuantityDict:
        """Return a dictionary with the proposed parameter values as
        ``Quantity`` objects with numerical magnitude."""
        return self.__values

    @property
    def symbols(self) -> SymbolQuantityDict:
        """Return dictionary of symbols representing all defined parameters,
        including those that are dependent symbols.
        """
        self.status.assure(ModelStatus.FINALISED, "querying symbols")
        if self.status is not ModelStatus.FINALISED:
            msg = "Symbols can only be queried when instance is prepared"
            raise RuntimeError(msg)
        return self.__symbols

    @property
    def free_symbols(self) -> SymbolQuantityDict:
        """Return dictionary of all symbols that represent free parameters.

        .. deprecated:: V0.1a1

            Don't need this, can be constructed as implemented from symbols and
            values, and it's only needed once (in prepare)::

                {name: self.symbols[name] for name in self.values}

            """
        self.status.assure(ModelStatus.FINALISED, "querying free symbols")
        return {name: self.symbols[name] for name in self.values}


class ParameterDefinition:
    """During the definition phase of a model, an instance of this class is
    used to define and connect the parameters of the module.
    This is only done once for each Model implementation.

    Once the mode instance is created, the object of this class is replaced
    with a :class:`ParameterHandler` object that handles parameter access
    from ouside the scope of the model itself.
    """

    def __init__(self):
        self.__symbols = {}
        self.__values = {}
        self.__required = set()
        self.status = ModelStatus.INTERFACE

    def define(self, name: str, value: float = None, unit: str = "dimless"):
        """Define a parameter in the context of the model. The name must be
        unique in the context. If value is not given (``None``), the parameter
        must be provided from the outside - and be compliant with the given
        unit of measurement. Otherwise, the value-unit pair defines the
        default value of the parameter, if not provided from parent module"""
        self.status.assure(ModelStatus.INTERFACE, "parameter declaration")
        if self.defined(name):
            raise ValueError(f"Parameter '{name}' already defined")

        if value is None:
            self.__required.add(name)
        else:
            self.__values[name] = Quantity(value, unit)

        self.__symbols[name] = SymbolQuantity(name, unit)

    @property
    def symbols(self) -> SymbolQuantityDict:
        """Return dictionary of symbols representing all defined parameters,
        including those that are dependent symbols.
        """
        return self.__symbols

    def __getitem__(self, name: str) -> SymbolQuantity:
        """Essential to use the defined parameters in the
        :meth:`Model.define <simu.model.Model.define>` method is to utilise
        this operator."""
        self.status.assure(ModelStatus.DEFINE, "parameter query")
        return self.__symbols[name]

    def defined(self, name: str) -> bool:
        """Check whether a parameter of given ``name`` is defined and return
        the result as a boolean."""
        return name in self.__symbols

    def create_handler(self) -> ParameterHandler:
        """After the definition phase, this method creates a
        :class:`ParameterHandler` object that handles further processing of
        the parameters in the model instance."""
        self.status.assure(ModelStatus.READY, "creating handler")
        return ParameterHandler(self.symbols, self.__values, self.__required)
