"""This module provides materials based on publically available models and
parameters"""

# stdlib modules
from pathlib import Path
from typing import Iterable, Tuple, List
from abc import ABC, abstractmethod

from ..utilities import SymbolQuantity, Quantity
from ..utilities.errors import DimensionalityError
from ..utilities.types import NestedStringDict, NestedQuantityDict

# external modules
from yaml import safe_load

# internal modules
from ..thermo import ThermoFactory, all_states, all_contributions

DATA_DIR = Path(__file__).resolve().parent / "data"


class AbstractThermoSource(ABC):
    @abstractmethod
    def __getitem__(self, path: List[str]) -> Quantity:
        pass


class ThermoPropertyStore:
    def __init__(self):
        self.__provided_parameters: NestedQuantityDict = {}
        self.__sources: List[AbstractThermoSource] = []

    def get_symbols(self, parameter_struct: NestedStringDict) \
            -> NestedQuantityDict:
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

        def prepare(name: str, key: str, query: NestedStringDict,
                    stored: NestedQuantityDict):
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

    def get_all_symbols(self) -> NestedQuantityDict:
        """This method returns all previously prepared symbols, cf.
        :meth:`get_symbols`, as a nested dictionary of symbolic quantities"""
        return self.__provided_parameters

    def get_all_symbol_values(self) -> NestedQuantityDict:
        """This method seeks in connected data sources for all previously
        prepared symbols.

        A ``KeyError`` is thrown if not all required parameters are available.
        Use
        """
        found, missing = self.__get_values()
        if missing:
            raise KeyError("Missing parameter values. Use " 
                           "'get_missing_symbols' to find out which")
        return found

    def get_missing_symbols(self) -> NestedQuantityDict:
        """This method tries to collect values for all previously prepared
        symbols and returns a nested dictionary of quantities with those not
        found."""
        return self.__get_values()[1]

    def __get_values(self) -> Tuple[NestedQuantityDict, NestedQuantityDict]:
        """Return a tuple of

          a. values found for previously defined symbols, and
          b. symbols for those entries where values are not found.
        """
        def get_value(path: List[str], qty: Quantity) -> Quantity:
            """Query the sources for value"""
            for source in self.__sources:
                try:
                    result = source[path]
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
                    return result
            else:
                raise KeyError("Parameter not found")

        def extract(path: List[str], struct: NestedQuantityDict) \
                -> Tuple[NestedQuantityDict, NestedQuantityDict]:
            """Recursive helper function to extract values"""
            try:
                items = struct.items()
            except AttributeError:  # found a leaf node
                try:
                    return get_value(path, struct), {}
                except KeyError:
                    return {}, struct

            # found a sub structure
            found: NestedQuantityDict = {}
            missing: NestedQuantityDict = {}
            for key, value in items:
                found[key], miss = extract(path + [key], value)
                if miss:
                    missing[key] = miss

            return found, missing

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
