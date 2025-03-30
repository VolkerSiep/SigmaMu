# -*- coding: utf-8 -*-

# internal modules
from simu import jacobian, ThermoContribution


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

    The derivative :math:`\partial y/\partial x` will be provided as
    ``f"d{options['y']}_d{options['x']}"``.
    """

    def define(self, res):
        independent = self.options["x"]
        dependent = self.options["y"]
        name = f"d{dependent}_d{independent}"
        if name not in res:
            res[name] = jacobian(res[dependent], res[independent])

        # TODO: what do I do if one or even both of these are vectors?
        #  if one is a vector, I can still use the normal declare_vector_keys,
        #  but if both are vectors, I need to generalize that routine.