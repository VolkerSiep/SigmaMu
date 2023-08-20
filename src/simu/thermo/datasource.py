"""This module defines an interface and a simple implementation class to
obtain sets of thermodynamic parameters. These data sources can be chained
for retrieving data with decreasing priority from different sources."""

from abc import ABC, abstractmethod

from ..utilities import Quantity
from ..utilities.types import (
    QuantityDict, NestedQuantityDict, NestedStringDict)


class ThermoDataSource(ABC):
    """This abstract class represents a data source for thermodynamic data
    sets."""
    def __int__(self):
        pass

    @abstractmethod
    def extract(self, structure: NestedStringDict) -> NestedQuantityDict:
        """Given a recursive structure, of which the leaf values are strings
        representing the units of measurements, extract the parameters from the
        database."""

    def extract_all_used_symbols(self) -> QuantityDict:
        """Return a flat dictionary with all symnbolic Quantities that have
        been extracted before"""

    def extract_all_used_values(self) -> QuantityDict:
        """Return a flat dictionary with all numeric Quantities that have
        been extracted before"""

    # need also something to update values

    def __setitem__(self, key: str, value: Quantity):
        """Set the value of a parameter to a new numerical value"""
