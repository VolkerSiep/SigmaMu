from ..utilities import SymbolQuantity, Quantity
from .utils import ModelStatus


class PropertyHandler:
    """After instantiating a model, this object handles the connections of the
    calculated properties to the model's client code (normally a parent model).
    """

    def __init__(self, symbols: dict[str, Quantity],
                 local_symbols: dict[str, Quantity]):
        """This constructor is called internally by the
        :meth:`PropertyDefinition.create_handler` method to transfer the
        necessary data."""
        self.__symbols = dict(symbols)  # copy required
        self.__local_symbols = dict(local_symbols)
        self.status = ModelStatus.READY

    def __getitem__(self, name: str) -> Quantity:
        """Before the model is finalised, this operator returns the empty
        quantity that represents the property as an argument to the model in
        the outer hierarchy layer.

        After the model is finalised, this operator returns the connected
        property as a result quantity for the overall model.
        """
        return self.__symbols[name]

    @property
    def symbols(self):
        """Before the model is finalised, this property holds the empty
        quantities to connect an outer model. After the model is finalised,
        it returns the properly connected properties to be calculated by an
        overall model."""
        return self.__symbols

    @property
    def local_symbols(self):
        """The local symbols as defined in :meth:`simu.Model.define`"""
        return self.__local_symbols

    def finalise(self, symbols_from_function):
        """The properties are finalised by replacing the place-holder of
        the symbols with the results of the function."""
        self.__symbols = symbols_from_function
        self.__local_symbols = None
        self.status = ModelStatus.FINALISED


class PropertyDefinition:
    """During the definition phase of a model, an object of this class is used
    to first declare the existance of the model properties, and then to receive
    their definition, such that the model then can include the parameters as
    part of its results."""

    def __init__(self):
        self.__symbols = {}
        self.__local_symbols = {}
        self.status = ModelStatus.INTERFACE

    def provide(self, name: str, unit: str):
        """During the interface definition phase, this method declares that
        a property of given ``name`` with a unit of measurement compatible to
        ``unit`` will be provided."""
        self.status.assure(ModelStatus.INTERFACE, "property declaration")
        if name in self.__symbols:
            raise ValueError(f"Property '{name}' already defined")
        self.__symbols[name] = SymbolQuantity(name, unit)

    def __setitem__(self, name: str, symbol: Quantity):
        """During the definition phase, quantities that are calculated as part
        of the implementation can be assigned.

        Each ``name`` can only be assigned once, and it must have been declared
        as to be provided (:meth:`provide`).
        """
        self.status.assure(ModelStatus.DEFINE, "property definition")
        if not name in self.__symbols:
            raise ValueError(f"Property '{name}' not defined")
        if name in self.__local_symbols:
            raise ValueError(f"Property '{name}' already assigned")
        prototype = self.__symbols[name]
        symbol.to(prototype.units)  # DimensionalityError if not compatible
        self.__local_symbols[name] = symbol

    def __getitem__(self, name: str) -> Quantity:
        """For convenience, the set property can be re-obtained from this
        object for further calculations during definition phase."""
        self.status.assure(ModelStatus.DEFINE, "property access")
        return self.__local_symbols[name]

    @property
    def local_symbols(self) -> dict[str, Quantity]:
        """To create the local function, the model requires to obtain a
        dictionary of all defined properties, as they are locally calculated.
        That dictionary is returned here."""
        self.status.assure(ModelStatus.FUNCTION, "symbol query")
        return self.__local_symbols

    def create_handler(self) -> PropertyHandler:
        """After the definition phase, this method creates a
        :class:`PropertyHandler` object that handles further processing of
        the properties in the model instance."""
        self.status.assure(ModelStatus.READY, "handler creation")
        # check if all symbols are defined
        missing = set(self.__symbols) - set(self.__local_symbols)
        if missing:
            missing = ", ".join(missing)
            raise ValueError(f"Promised properties not assigned: {missing}")
        return PropertyHandler(self.__symbols, self.__local_symbols)
