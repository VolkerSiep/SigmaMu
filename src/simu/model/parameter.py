"""This module implements functionality related to parameter handling"""

from typing import Optional, Iterator
from collections.abc import Mapping, Iterable

from ..utilities import Quantity, SymbolQuantity
from ..utilities.types import QuantityDict
from ..utilities.errors import DataFlowError


class ParameterHandler(Mapping):
    """This class, being instantiated as the :attr:`Model.parameters`
    attribute, allows to define and access process parameters.

    During :meth:`Model.initialse`, the model defines parameters with and
    without pre-defined values.

    Within the ``with`` block, the parent module connects parameters by
    providing external symbols. Here, the unit of measurement must be
    compatibe to the parameters' definition.

    During :meth:`Model.define` checks that all required symbols (the ones
    without default values) are connected. It keeps track of the non-connected
    and static parameters, as they will need to be used as arguments for the
    overall model function.
    """

    static_parameters: QuantityDict = {}
    static_values: QuantityDict = {}

    def __init__(self, static_id: str):
        self.__params: QuantityDict = {}
        self.__values: QuantityDict = {}
        self.__static_used_names: set[str] = set()
        self.__static_id = static_id

    def __len__(self) -> int:
        return len(self.__params)

    def __iter__(self) -> Iterator[str]:
        return iter(self.__params)

    def define(self,
               name: str,
               value: Optional[float] = None,
               unit: Optional[str] = "dimless") -> None:
        """Define a parameter from within the :meth:`Model.interface`
        method."""

        if name in self.__params:
            raise KeyError(f"Parameter '{name}' already defined")

        # all parameters are stored as symbols with unit
        self.__params[name] = SymbolQuantity(name, unit)
        if value is not None:
            # these have default values and don't need to be provided
            self.__values[name] = Quantity(value, unit)

    def static(self,
               name: str,
               value: float,
               unit: Optional[str] = "dimless") -> SymbolQuantity:
        """Define a static parameter over all instances of the associated
        model."""
        if name in self.__static_used_names:
            raise KeyError(f"Static Parameter '{name}' already defined")
        full_name = f"{self.__static_id}/{name}"
        cls = ParameterHandler
        if full_name not in cls.static_parameters:
            result = SymbolQuantity(name, unit)
            cls.static_parameters[full_name] = result
            cls.static_values[full_name] = Quantity(value, unit)
            return result
        else:
            return cls.static_parameters[full_name]

    def __getitem__(self, name: str) -> Quantity:
        """Return the symbol of a defined parameter. This is to be used within
        the :meth:`Model.define` method."""
        return self.__params[name]

    def get_static(self, name: str) -> SymbolQuantity:
        """Return a static parameter of given name"""
        cls = ParameterHandler
        full_name = f"{self.__static_id}/{name}"
        return cls.static_parameters[full_name]

    def static_names(self) -> Iterable[str]:
        """An iterator over all names of defined static parameters"""
        params = ParameterHandler.static_parameters
        prefix = f"{self.__static_id}/"
        return filter(lambda key: key.startswith(prefix), params.keys())

    # @property
    # def all(self) -> QuantityDict:
    #     """Return a dictionary of all defined parameters on this level"""
    #     return dict(self.__params)
    #
    # def values(self) -> QuantityDict:
    #     """Return the (default) values for the defined parameters that have
    #     such default value."""
    #     return dict(self.__values)  # create copy to preserve own set

    def create_proxy(self) -> "ParameterProxy":
        """Create a proxy object for configuration in hierarchy context"""
        return ParameterProxy(self, self.__params, self.__values)


class ParameterProxy:
    """This class is instantiated by the parent's :class:`ParameterHandler`
    to configure the parameter connections from the parent context."""

    def __init__(self, handler: ParameterHandler,
                 params: QuantityDict, values: QuantityDict):
        self.__handler = handler
        self.__model_name = "N/A"

        self.__params = params  # reference to dicts in handler
        self.__values = values  # not a copy by design
        self.__free: QuantityDict = {}

        self.__provided: set[str] = set()

    def set_name(self, name: str):
        """Set the name of the model for better error diagnostics"""
        self.__model_name = name

    def provide(self, **kwargs: Quantity) -> None:
        """Connect a parameter from parent context to child parameter."""
        for name, quantity in kwargs.items():
            self.__assure_ok(name, quantity)
            self.__params[name] = quantity
            self.__provided.add(name)

    def update(self, name: str, value: float, unit: str) -> None:
        """Provide new default value for child parameter"""
        quantity = Quantity(value, unit)
        self.__assure_ok(name, quantity)
        self.__values[name] = quantity

    def __assure_ok(self, name: str, quantity: Quantity) -> None:
        """Check that a new quantity can be processed and is compatible with
        the previously defined slot."""
        model_name = self.__model_name
        if name not in self.__params:
            msg = f"Parameter '{name}' not defined in '{model_name}'"
            raise KeyError(msg)
        if name in self.__provided:
            msg = f"Parameter '{name}' already provided in '{model_name}'"
            raise KeyError(msg)
        # check unit compatibility -> throw exception on demand:
        quantity.to(self.__params[name].units)

    @property
    def free(self):
        """Symbols representing the symbols of the parameters that have not
        been provided. These must be used to create an overall function, when
        parameter values are provided from the outside."""
        return dict(self.__free)

    @property
    def values(self):
        """Symbols representing the values of the parameters that have not
        been provided. These must be used to call the overall function."""
        return dict(self.__values)

    def finalise(self):
        """Make sure there are values for all non-provided parameters with
        no value. Remove values of provided parameters"""

        # are all required parameters connected?
        provided = self.__provided
        params = self.__params
        values = set(self.__values.keys())
        all_param = set(params.keys())
        missing = all_param - (provided | values)
        if missing:
            lst = ", ".join(map(lambda m: f"'{m}'", missing))
            name = self.__model_name
            msg = f"Model '{name}' has unresolved parameters: {lst}"
            raise DataFlowError(msg)

        # clean values of provided parameters
        self.__values = {
            name: quantity
            for name, quantity in self.__values.items() if name not in provided
        }

        # create symbols for still free variables (for later overall function)
        self.__free = {
            name: value for name, value in params.items()
            if name in self.__values
        }
