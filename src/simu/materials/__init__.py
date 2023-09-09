"""This module provides materials based on publically available models and
parameters"""

# stdlib modules
from pathlib import Path
from typing import Iterable, Tuple, List, Union, Mapping
from abc import ABC, abstractmethod

from ..utilities import SymbolQuantity, Quantity
from ..utilities.qstructures import parse_quantities_in_struct
from ..utilities.errors import DimensionalityError

# external modules
from yaml import safe_load

# internal modules
from ..thermo import ThermoFactory, all_states, all_contributions


NestedQuantityMapping = Union[Quantity, Mapping[str, "NestedQuantityMapping"]]
NestedStringMapping = Union[str, Mapping[str, "NestedStringMapping"]]

_RT = Tuple[NestedQuantityMapping, NestedStringMapping, NestedStringMapping]

DATA_DIR = Path(__file__).resolve().parent / "data"


class AbstractThermoSource(ABC):
    """Any source of thermodynamic parameters is to return a Quantity object
    if the path describes an available property.
    """
    @abstractmethod
    def __getitem__(self, path: Iterable[str]) -> Quantity:
        pass


class NestedDictThermoSource(AbstractThermoSource):
    """A source of thermodynamic parameters defined by a nested dictionary
    of Quantities"""
    def __init__(self, data: NestedQuantityMapping):
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
    def __init__(self, data: NestedStringMapping):
        super().__init__(parse_quantities_in_struct(data))


class ThermoPropertyStore:
    """This class connects both to the thermodynamic model instances by
    providing the parameters, and to the model's numerical interface by
    providing the symbols and values of all used parameters, allowing those
    to be altered dynamically and optimised on."""
    def __init__(self):
        self.__provided_parameters: NestedQuantityMapping = {}
        self.__sources: Mapping[str, AbstractThermoSource] = {}

    def get_symbols(self, parameter_struct: NestedStringMapping) \
            -> NestedQuantityMapping:
        """Query the nested structure of parameter symbols from the store.
        The ``param_struct`` parameter is a nested dictionary with leaf entries
        being the unit of measurement of the thermodynamic parameters defined
        by the path of keys to address them.

        The method returns a dictionary of the same structure, but with
        symbolic quantities as leaf values.

        With multiple calls to this method with varying parameter structs,
        previously defined symbols will be reused, and a
        ``DimensionialityError`` is raised if such previously defined symbol
        is incompatible with respect to the physical dimension.
        """

        def prepare(name: str, key: str, query: NestedStringMapping,
                    stored: NestedQuantityMapping):
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

    def get_all_symbols(self) -> NestedQuantityMapping:
        """This method returns all previously prepared symbols, cf.
        :meth:`get_symbols`, as a nested dictionary of symbolic quantities"""
        return self.__provided_parameters

    def get_all_symbol_values(self) -> NestedQuantityMapping:
        """This method seeks in connected data sources for all previously
        prepared symbols.

        A ``KeyError`` is thrown if not all required parameters are available.
        Use
        """
        found, missing, sources = self.__get_values()
        if missing:
            raise KeyError("Missing parameter values. Use " 
                           "'get_missing_symbols' to find out which")
        return found

    def get_missing_symbols(self) -> NestedStringMapping:
        """This method tries to collect values for all previously prepared
        symbols and returns a nested dictionary of quantities with those not
        found."""
        return self.__get_values()[1]

    def get_sources(self) -> NestedStringMapping:
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
        def get_value(path: List[str], qty: Quantity) -> Tuple[Quantity, str]:
            """Query the sources for value"""
            for source_name, source in self.__sources.items():
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

        def extract(path: List[str], struct: NestedQuantityMapping) -> _RT:
            """Recursive helper function to extract values"""
            try:
                items = struct.items()
            except AttributeError:  # found a leaf node
                try:
                    value, source = get_value(path, struct)
                    return value, {}, source
                except KeyError:
                    return {}, f"{struct.units:~}", {}

            # found a sub structure
            found: NestedQuantityMapping = {}
            missing: NestedStringMapping = {}
            source: NestedStringMapping = {}
            for key, value in items:
                found[key], miss, source[key] = extract(path + [key], value)
                if miss:
                    missing[key] = miss

            return found, missing, source

        return extract([], self.__provided_parameters)


class ExampleThermoFactory(ThermoFactory):
    """This ThermoFactory subclass is capable of creating frames from the base
    SiMu installation, hence thermodynamic models that are found in open
    literature."""
    def __init__(self):
        """Default and only constructor"""
        ThermoFactory.__init__(self)
        self.register(*all_contributions)
        for state in all_states:
            self.register_state_definition(state)

        with open(DATA_DIR / "structures.yml", encoding='UTF-8') as file:
            self.__structures = safe_load(file)
        with open(DATA_DIR / "definitions.yml", encoding='UTF-8') as file:
            self.__configurations = safe_load(file)

    @property
    def configuration_names(self) -> Iterable[str]:
        return self.__configurations.keys()

    def create_frame(self, configuration: str):
        cfg = self.__configurations[configuration]
        return super().create_frame(cfg | self.__structures[cfg["structure"]])
