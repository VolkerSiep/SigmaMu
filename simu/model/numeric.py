"""The module defining the :class:`NumericHandler` class"""
from ..utilities import flatten_dictionary, QFunction, Quantity, SymbolQuantity
from .utils import ModelStatus

class NumericHandler:
    """The numeric handler (``NumericHandler``) is the object that represents the
    interface for all numerical calculations. It treats the ``Model`` object that
    it is associated to as the top level model, and includes all sub models
    recursively."""

    def __init__(self, parent):
        self.__parent = parent
        self.__result = None
        self.__function = None
        self.status = ModelStatus.READY


    def prepare(self) -> "NumericHandler":
        """This method assumes the associated model instance to be the top
        model and creates a numerical object to represent that model.

        This method returns the ``NumericHandler`` object itself for
        convenience.
        """
        # create overall function
        parameters = self.all_parameter_symbols
        properties = self.all_property_symbols
        # TODO: add residuals and thermodynamic states / parameters

        args = {"model_parameters": parameters}
        results = {"properties": properties}

        # TODO: do I need to make the secondary entries completely flat
        #   for casadi (to avoid thousands of arguments and results)
        self.__function = QFunction(args, results, "process_model")
        return self

    @property
    def parameters(self) -> dict[str, Quantity]:
        """After successful evaluation (:meth:`evaluate`), this method returns
        a flat dictionary of calculated properties."""
        params = flatten_dictionary({
            name: child.numerics.parameters
            for name, child in self.__parent.hierarchy.items()
        })
        params.update(self.__parent.parameters.values)
        return params

    def evaluate(self):
        """After calling :meth:`prepare`, the constructed function
        :meth:`function` is hereby called with the current state of the model,
        meaning parameters and independent variables.
        """
        args = {"model_parameters": self.parameters}
        self.__result = self.__function(args)

    @property
    def function(self) -> QFunction:
        """After calling :meth:`prepare`, this property holds a casadi
        function that maps the model parameters and independent variables
        to calculated properties and residuals."""
        return self.__function

    @property
    def properties(self) -> dict[str, Quantity]:
        """Return the calculated model properties after the model is evaluated.
        """
        return self.__result["properties"]

    @property
    def all_property_symbols(self) -> dict[str, SymbolQuantity]:
        """Return the property symbols of the associated model and all sub
        modules as a flattened dictionary with symbolic quantities as values
        This property is mainly intented to aid constructing an overall
        :class:`QFunction` object that represents the model
        """
        prop = flatten_dictionary({
            name: child.numerics.all_property_symbols
            for name, child in self.__parent.hierarchy.items()
        })
        prop.update(self.__parent.properties.symbols)
        return prop

    @property
    def all_parameter_symbols(self) -> dict[str, SymbolQuantity]:
        """Return the parameter symbols of the associated model and all sub
        modules as a flattened dictionary with symbolic quantities as values
        This property is mainly intented to aid constructing an overall
        :class:`QFunction` object that represents the model
        """
        param = flatten_dictionary({
            name: child.numerics.all_parameter_symbols
            for name, child in self.__parent.hierarchy.items()
        })
        param.update(self.__parent.parameters.symbols)
        return param
