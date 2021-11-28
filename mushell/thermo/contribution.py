# -*- coding: utf-8 -*-

"""This module defines the class :cls:`ThermoContribution`, which defines
the building blocks of a :cls:`ThermoFrame` function object."""

from abc import ABC, abstractmethod

# external modules
from casadi import vertcat


class ThermoContribution(ABC):
    """This abstract class defines the interface of a contribution to a
    thermodynamic state function, as collected in :cls:`ThermoFrame` objects.

    The definition is based on the ``casadi`` library, and its definition
    is required to build a ``casadi`` evaluation structure as the
    implementation of the belonging equations.

    The usage of this class is mainly indirect by instantiation via the
    :cls`ThermoFactory` objects and parametrisation via the provided parameter
    structures.
    """
    def __init__(self, species, options):
        self.species = species
        self.options = options

    @abstractmethod
    def define(self, res, par):
        """Abstract method that implements the ``casadi`` expressions
        that make up this contribution."""
        pass

    @property
    def parameter_structure(self):
        return {}

    def relax(self, state, delta_state, parameters):
        """Virtual function to report the maximal allowable step size in
        the state variables"""
        # default implementation
        return 100  # a number greater than 1 / gamma for practical gamma

    @staticmethod
    def _tensor_structure(*keys):
        if not keys:
            return None
        current, *rest = keys
        func = ThermoContribution._tensor_structure
        return {k: func(*rest) for k in current}

    def _vector(self, dictionary, keys=None):
        if keys is None:
            keys = self.species
        return vertcat(*[dictionary[key] for key in keys])

    def initial_state(self, T, p, n, parameters):
        """To be overwritten if a contribution can initialise the state"""
        return None
