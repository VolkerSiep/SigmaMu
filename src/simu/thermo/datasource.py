"""This module defines an interface and a simple implementation class to
obtain sets of thermodynamic parameters. These data sources can be chained
for retrieving data with decreasing priority from different sources."""

from abc import ABC, abstractmethod


class ThermoDataSource(ABC):
    """This abstract class represents a data source for thermodynamic data
    sets."""
    def __int__(self):
        pass

    @abstractmethod
    def extract(self, parameter_structure: dict):
        """Given a recursive structure, of which the leaf values are strings
        representing the units of measurements, extract the parameters from the
        database."""
