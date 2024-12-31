"""module containing classes to obtain thermodynamic parameters from various
sources."""

from typing import Iterable, Collection
from abc import ABC, abstractmethod

from pint import DimensionalityError

from simu.core.utilities import (
    Quantity, parse_quantities_in_struct, SymbolQuantity)
from simu.core.utilities.types import NestedMap, NestedMutMap, MutMap

_RT = tuple[Quantity | NestedMap[Quantity],
            str | NestedMap[str],
            str | NestedMap[str]]


class AbstractThermoSource(ABC):
    """Any source of thermodynamic parameters is to return a Quantity object
    if the path describes an available property.

    Implementations of this class can be added to the
    :class:`ThermoParameterStore` to represent a source of thermodynamic
    parameters. It is advisable to arrange source objects by
    *bibliographic* sources.

    .. note::

        We would have loved to download large databases of thermodynamic
        properties and to publish them with this software. Unfortunate for this
        particular aspect, there are terms like *copyright* and *license* that
        prohibit this.

        While it is completely legal to obtain thermodynamic parameters for free
        from many sources for the sake of working with them, re-publishing
        compiled databases under the LGPL license is another thing.
    """
    @abstractmethod
    def __getitem__(self, path: Iterable[str]) -> Quantity:
        """
        This operator is used to extract a value from the source, addressed by
        a sequence of keys to navigate into the nested data structure:

        >>> struct = {"T": {"H2O": "100 K", "NH3": "200 K"},
        ...           "p": {"H2O": "10 bar", "NH3": "20 atm"}}
        >>> source = StringDictThermoSource(struct)
        >>> print(source["T", "H2O"])
        100 kelvin

        A ``KeyError`` is raised if either a key does not exist or if the final
        key still has sub-keys.
        """
        pass


class NestedDictThermoSource(AbstractThermoSource):
    """A source of thermodynamic parameters defined by a nested dictionary
    of Quantities"""
    def __init__(self, data: NestedMap[Quantity]):
        self.__data = data

    def __getitem__(self, path: Iterable[str]) -> Quantity:
        result = self.__data
        for key in path:
            result = result[key]
        try:
            _ = result.magnitude, result.units  # is it a Quantity?
        except AttributeError:
            raise KeyError("Key doesn't point to Quantity") from None
        return result


class StringDictThermoSource(NestedDictThermoSource):
    """Source of thermodynamic parameters defined by a nested dictionary.
    The leaf entries of the nested dictionary are strings representing
    the quantities, such as "120 K" or "32.1 kJ/(mol*K)"."""
    def __init__(self, data: NestedMap[str]):
        super().__init__(parse_quantities_in_struct(data))


class ThermoParameterStore:
    """This class connects both to the thermodynamic model instances by
    providing the parameters, and to the model's numerical interface by
    providing the symbols and values of all used parameters, allowing those
    to be altered dynamically and optimised on.

    The important aspect here is that a process model can have hundreds of
    thermodynamic model instances, and while they do not need to, most of them
    usually share the same parameters.

    Adapted to the ``SiMu`` data flow, the initial queries by the models to
    fetch parameters are of symbolic nature, using :meth:`get_symbols` and
    :meth:`get_all_symbols`. This happens during the model definition and
    assures that multiple model instances that ask for the same parameter
    receive the same symbol, such that parameter optimisation tasks can be
    accomplished easily.

    In the (numerical) evaluation phase, the store is asked for the actual
    values via :meth:`get_all_values`. In most cases these are used as constants
    in the calculations, but parameter estimation is possible by keeping
    selected parameters as free variables.
    """

    name: str

    def __init__(self):
        self.__provided_parameters: NestedMutMap[Quantity] = {}
        self.__sources: MutMap[AbstractThermoSource] = {}
        self.name = "default"

    def get_symbols(self, parameter_struct: NestedMap[str]) \
            -> NestedMap[Quantity]:
        """Query the nested structure of parameter symbols from the store.
        The ``parameter_struct`` parameter is a nested dictionary with leaf
        entries being the unit of measurement of the thermodynamic parameters
        defined by the path of keys to address them.

        The method returns a dictionary of the same structure, but with
        symbolic quantities as leaf values.

        With multiple calls to this method with varying parameter structs,
        previously defined symbols will be reused, and a
        ``DimensionalityError`` is raised if such previously defined symbol
        is incompatible with respect to the physical dimension.
        """

        def prepare(name: str, key: str, query: NestedMap[str],
                    stored: NestedMutMap[Quantity]):
            """Helper function to recursively retrieve and define symbols"""
            name = f"{name}.{key}" if name else key
            try:
                items = query.items()
            except AttributeError:
                if key in stored:
                    # compare unit for compatibility
                    try:
                        stored[key].to(query)  # convert unit, see if it works
                    except DimensionalityError as err:
                        err.extra_msg = \
                            " - Error fetching previously defined thermo " \
                            f"parameter '{name}'."
                        raise err from None
                else:
                    stored[key] = SymbolQuantity(name, query)
                return stored[key]

            if key not in stored:
                stored[key] = {}

            return {k: prepare(name, k, q, stored[key]) for k, q in items}

        return {k: prepare("", k, s, self.__provided_parameters)
                for k, s in parameter_struct.items()}

    def get_all_symbols(self) -> NestedMap[Quantity]:
        """This method returns all previously prepared symbols, cf.
        :meth:`get_symbols`, as a nested dictionary of symbolic quantities"""
        return self.__provided_parameters

    def get_all_values(self) -> NestedMap[Quantity]:
        """This method seeks in connected data sources for all previously
        prepared symbols.

        A ``KeyError`` is thrown if not all required parameters are available.
        Use :meth:`get_missing_symbols` to get the structure of missing
        parameters.

        The sources are queried sequentially in the reverse sequence of them
        being added, i.e. latest added Sources supersede previous values of
        same parameters.
        """
        found, missing, sources = self.__get_values()
        if missing:
            raise KeyError("Missing parameter values. Use " 
                           "'get_missing_symbols' to find out which")
        return found

    def get_missing_symbols(self) -> NestedMap[str]:
        """This method tries to collect values for all previously prepared
        symbols and returns a nested dictionary of quantities with those not
        found."""
        return self.__get_values()[1]

    def get_sources(self) -> NestedMap[str]:
        """This method returns the name of the data source in which the
        individual values are found"""
        return self.__get_values()[2]

    def add_source(self, name: str, source: AbstractThermoSource):
        """Add a given source of thermodynamic parameters to the store.
        """
        if name in self.__sources:
            raise KeyError(f"Source '{name}' already defined")
        self.__sources[name] = source

    def __get_values(self) -> _RT:
        """Return a tuple of

          a. values found for previously defined symbols, and
          b. symbols for those entries where values are not found.
          c. names of the sources where the values are found
        """
        def get_value(path: Collection[str],
                      qty: Quantity) -> tuple[Quantity, str]:
            """Query the sources for value"""
            for source_name, source in reversed(self.__sources.items()):
                try:
                    result = source[*path]
                except KeyError:
                    continue
                try:
                    result.to(qty.units)
                except DimensionalityError as err:
                    err.extra_msg = \
                        " - Error fetching thermodynamic property " + \
                        ".".join(path)
                    raise err from None
                else:
                    return result, source_name
            else:
                raise KeyError("Parameter not found")

        def extract(path: list[str],
                    struct: Quantity | NestedMutMap[Quantity]) -> _RT:
            """Recursive helper function to extract values"""
            try:
                items = struct.items()
            except AttributeError:  # found a leaf node
                try:
                    value, source = get_value(path, struct)
                    return value, {}, source
                except KeyError:
                    # the space after the unit is important here, for the
                    # string not to evaluate to false in case of dimensionless
                    # parameters!
                    return {}, f"{struct.units:~} ", {}

            # found a sub structure
            found: NestedMap[Quantity] = {}
            missing: NestedMutMap[str] = {}
            source: NestedMap[str] = {}
            for key, value in items:
                found_i, miss, source_i = extract(path + [key], value)
                found[key], source[key] = found_i, source_i
                if miss:
                    missing[key] = miss

            return found, missing, source

        return extract([], self.__provided_parameters)
