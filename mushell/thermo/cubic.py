
# -*- coding: utf-8 -*-

# external modules
from casadi import dot, log, vertsplit, vertcat, sum1

# internal modules
from .contribution import ThermoContribution
from ..constants import R_GAS


class RedlichKwongEOS(ThermoContribution):
    r"""This contribution implements a general Redlich-Kwong equation of state
    with Peneloux volume translation:

    .. math::

        p = \frac{N\,R\,T}{V - B + C} - \frac{A}{(V + C)\,(V + B + C)}

    The following properties need to be procided by upstream contributions:

    ======== ===================================================
    Property Symbol
    ======== ===================================================
    RK_A     :math:`A`
    RK_A_T   :math:`A_T:=\partial A/\partial T|_{V,n}`
    RK_A_N   :math:`A_n:=\partial A/\partial n|_{T, V}` (vector)
    RK_B     :math:`B`
    RK_B_T   :math:`B_T:=\partial B/\partial T|_{V,n}`
    RK_B_N   :math:`B_n:=\partial B/\partial n|_{T, V}` (vector)
    RK_C     :math:`C`
    RK_C_T   :math:`C_T:=\partial C/\partial T|_{V,n}`
    RK_C_N   :math:`C_n:=\partial C/\partial n|_{T, V}` (vector)
    ======== ===================================================

    Care is to be taken when utilising a temperature-dependent :math:`C`
    contribution, as doing so can have significant effects on the calorimetric
    properties.

    As such, there are no further model parameters to be provided at this
    point. The residual Helmholtz function is

    .. math::
        A^\mathrm{res} = \int\limits_V^\infty
           p - \frac{N\,R\,T}{V} \mathrm{d}V
         = N\,R\,T\,\ln \frac{V}{V + C - B} +
           \frac{A}{B}\,\ln \frac{V + C}{V + C + B}

    The implemented temperature derivative is

    .. math::

        -S^\mathrm{res} &= N\,R\,\left [
            \ln\frac{V}{V + B + C} + T\,\frac{B_T + C_T}{V + B + C}
            \right ]\\& +
            \frac1B\left (A_T + \frac{A}{B}\,B_T\right )\,
            \ln \frac{V + C}{V + B + C} +
            \frac{A}{B}\left  [
                \frac{C_T}{V + C} - \frac{B_T + C_T}{V + B + C}
            \right ]

    The volume derivative is the negative residual pressure:

    .. math::

        -p^\mathrm{res} =
          N\,R\,T\, \left [ \frac1V - \frac1{V - B + C}\right ] +
          \frac{A}{(V + C)\,(V + B + C)}

    The derivative with respect to molar quantities is

    .. math::

        \mu_i^\mathrm{res} &= R\,T\,\left [
            \ln\frac{V}{V + B + C} + N\,\frac{B_{n,i} + C_{n,i}}{V + B + C}
            \right ]\\& +
            \frac1B\left (A_{n,i} + \frac{A}{B}\,B_{n,i}\right )\,
            \ln \frac{V + C}{V + B + C} +
            \frac{A}{B}\left  [
                \frac{C_{n,i}}{V + C} - \frac{B_{n,i} + C_{n,i}}{V + B + C}
            \right ]

    The contribution updates are

    .. math::

        S &\leftarrow S + S^\mathrm{res}\\
        p &\leftarrow p + p^\mathrm{res}\\
        \mu_i &\leftarrow \mu_i + \mu_i^\mathrm{res}\\
    """
    @property
    def parameter_structure(self):
        t_s = ThermoContribution._tensor_structure
        # TODO: implement

    def define(self, res, par):
        # TODO: implement
        pass

    def relax(self, current_result, delta_state):
        """A proper relaxation strategy for the cubic equation of state is
        possible, but not trivial. A series of constraints apply. Firstly,
        the volume cannot become less than :math:`B - C`
        """
        # TODO: implement (V limits)
        pass




# TODO: derive for liquid and gas, different initialise and relax method
