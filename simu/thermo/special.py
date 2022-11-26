# -*- coding: utf-8 -*-

# external modules
from casadi import jacobian

# internal modules
from .contribution import ThermoContribution


class Derivative(ThermoContribution):
    r"""This auxiliary contribution provides the derivative of an arbitrary
    property with respect to an independent (state) variable. This can for
    instance be used to equip a thermodynamic model with extra temperature
    derivatives for calculating heat capacity or partial molar enthalpy
    as canonical properties.

    The contribution requires an ``option`` dictionary with the following
    entries:

        - ``x``: The name of the independent property :math:`x`
        - ``y``: The name of the dependent property :math:`y`

    The derivative :math:`\partial y/\partial x` will be provied as
    ``f"d{options['y']}_d{options['x']}"``.
    """

    def define(self, res, par):
        independent = self.options["x"]
        dependent = self.options["y"]
        name = f"d{dependent}_d{independent}"
        if name not in res:
            res[name] = jacobian(res[dependent], res[independent])
