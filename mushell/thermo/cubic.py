
# -*- coding: utf-8 -*-

# external modules
from casadi import dot, log, vertsplit, vertcat, sum1, SX

# internal modules
from .contribution import ThermoContribution
from ..constants import R_GAS


class RedlichKwongEOS(ThermoContribution):
    r"""This contribution implements a general Redlich-Kwong equation of state
    with Peneloux volume translation:

    .. math::

        p = \frac{N\,R\,T}{V - B + C} - \frac{A}{(V + C)\,(V + B + C)}

    The following properties need to be provided by upstream contributions:

    ======== ===================================================
    Property Symbol
    ======== ===================================================
    RK_A     :math:`A`
    RK_A_T   :math:`A_T:=\partial A/\partial T|_{V,n}`
    RK_A_N   :math:`A_n:=\partial A/\partial n|_{T, V}` (vector)
    RK_B     :math:`B`
    RK_B_T   :math:`B_T:=\partial B/\partial T|_{V,n}`
    RK_B_N   :math:`B_n:=\partial B/\partial n|_{T, V}` (vector)
    CEOS_C   :math:`C`
    CEOS_C_T :math:`C_T:=\partial C/\partial T|_{V,n}`
    CEOS_C_N :math:`C_n:=\partial C/\partial n|_{T, V}` (vector)
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
    def __init__(self, *args, **kwargs):
        """The constructor of this class is disabled, as an instance would
        not have a decent initialiation nor relaxation functionality and
        would therefore be without value. Instantiate objects of
        :class:`RedlichKwongEOSLiquid` or :class:`RedlichKwongEOSLiquid`
        instead."""
        raise NotImplementedError("Don't do this!")

    @property
    def parameter_structure(self):
        t_s = ThermoContribution._tensor_structure
        # TODO: implement

    def define(self, res, par):
        # TODO: implement
        pass


class RedlichKwongEOSLiquid(RedlichKwongEOS):
    """As a sub-class of :class:`RedlichKwongEOS`, this entity specialises
    on describing liquid (and super-critical) phases. The distinct elements
    are the initialisation and the relaxation, ensuring the state to be within
    the correct domain with respect to volume."""
    def relax(self, current_result, delta_state):
        """A proper relaxation strategy for the cubic equation of state is
        possible, but not trivial. A series of constraints apply. Firstly,
        the volume cannot become less than :math:`B - C`. Further, the
        calculated pressure shall not become negative, and thirdly, the
        derivative of pressure with respect to volume shall not become
        positive.
        """
        # TODO: implement (V limits)
        pass

    def initial_state(self, T, p, n, parameters):
        # same problem as with relax.
        # Maybe I need to calculate for T, V, n where V = NRT/p.
        # Then I have the EOS equation and can find the roots.

        # giving T, p, and n as initial results, try to evaluate all
        # contributions. Some will failes, as they might need V or x, but
        # the ones that succeed will fill some of the results.
        #
        # Then, start to ask the contributions top down to initialise.
        # The first one not returning None determines the result.
        pass

class RedlichKwongEOSGas(RedlichKwongEOS):
    """As a sub-class of :class:`RedlichKwongEOS`, this entity specialises
    on describing gas (and super-critical) phases. The distinct elements
    are the initialisation and the relaxation, ensuring the state to be within
    the correct domain with respect to volume."""
    def relax(self, current_result, delta_state):
        """A proper relaxation strategy for the cubic equation of state is
        possible, but not trivial. A series of constraints apply. Firstly,
        the volume cannot become less than :math:`B - C`. Further, the
        derivative of pressure with respect to volume shall not become
        positive.
        """
        # TODO: implement (V limits)
        pass

    def initial_state(self, T, p, n, parameters):
        pass


class LinearPenelouxVolumeShift(ThermoContribution):
    r"""The Peneloux volume shift provides the ``CEOS_C`` parameter and its
    temperature and compositional derivatives ``CEOS_C_T`` and ``CEOS_C_N``.
    As such, this contribution determines the `C` parameter in cubic equations
    of state, representing the volume translation.

    This class implements only a volume shift linear in molar quantities and
    not dependent on temperature - hereby not triggering impact on the
    caliometric properties of the EOS:

    .. math::

        C = \sum_i c_i\,n_i
    """
    @property
    def parameter_structure(self):
        t_s = ThermoContribution._tensor_structure
        return {"c_i": t_s(self.species)}

    def define(self, res, par):
        c_i = self._vector(par["c_i"])
        res["CEOS_C"] = dot(c_i, res["n"])
        res["CEOS_C_T"] = SX.zeros(1,1)
        res["CEOS_C_N"] = c_i


class NonSymmmetricMixingRule(ThermoContribution):
    r"""The mixing rule combines the pure species a-contributions (:math:`a_i`)
    into the ``RK_A`` parameter of the :class:`RedlichKwongEOS` class,
    including the derivatives ``RK_A_T`` and ``RK_A_N``.

    The following properties need to be provided by upstream contributions:

    ======== ===================================================
    Property Symbol
    ======== ===================================================
    RK_a     :math:`a` (vector)
    RK_a_T   :math:`a_T:=\partial a/\partial T|_{V,n}` (vector)
    RK_a_N   :math:`a_n:=\partial a/\partial n|_{T, V}` (matrix)
    ======== ===================================================

    Implemented with sparse binary interaction coefficient structure, the
    contribution includes a symmetric interaction (:math:`k`) and an
    anti-symmetric contribution (:math:`l`).

    .. math::

        A &= \sum_{ij} \sqrt{a_i(T)\,a_j(T)} \left [
               x_i\,x_j - 2\,x_i\,x_j\,k_{ij}(T) -
               2\,(x_i^2\,x_j - x_i\,x_j^2)\,l_{ij}(T)
               \right ]\\
          &= \frac1{N}\sum_{ij} \sqrt{a_i(T)\,a_j(T)} \left [
               n_i\,n_j - 2\,n_i\,n_j\,k_{ij}(T) -
               \frac2N\,(n_i^2\,n_j - n_i\,n_j^2)\,l_{ij}(T)
               \right ]\\

    The following properties will then be generated:

    ======== ===================================================
    Property Symbol
    ======== ===================================================
    RK_A     :math:`A`
    RK_A_T   :math:`A_T:=\partial A/\partial T|_{V,n}`
    RK_A_N   :math:`A_n:=\partial A/\partial n|_{T, V}` (vector)
    ======== ===================================================

    The binary interaction coefficients can be a function of temperature
    according to

    .. math::

        k_{ij}(T) &:= k_{0,ij} +
                      k_{1,ij}\,\left (1 - \frac{T}{T_\mathrm{ref}} \right ) +
                      k_{2,ij}\,\left (1 - \frac{T_\mathrm{ref}}{T} \right )\\
        l_{ij}(T) &:= l_{0,ij} +
                     l_{1,ij}\,\left (1 - \frac{T}{T_\mathrm{ref}} \right ) +
                     l_{2,ij}\,\left (1 - \frac{T_\mathrm{ref}}{T} \right )

    Each pair of species can only be define once. The pre-factor 2 is used to
    yield the same parameter values as if :math:`k` was a symmetric matrix and
    :math:`l` was anti-symmetric.

    Their temperature derivatives are

    .. math::

        k_{T,ij}(T) &:= \frac{1}{T_\mathrm{ref}} \left ( k_{1,ij} +
                        k_{2,ij}\,\frac{T_\mathrm{ref}^2}{T^2} \right )\\
        l_{T,ij}(T) &:= \frac{1}{T_\mathrm{ref}} \left ( l_{1,ij} +
                        l_{2,ij}\,\frac{T_\mathrm{ref}^2}{T^2} \right )

    For :math:`A_T`, we firstly define

    .. math::
        \mathcal A_{ij}   &:= \sqrt{a_i(T)\,a_j(T)}\\
        \mathcal A_{T,ij} &:= \hat a + \hat a^\mathrm{T}\ \text{with}
          \ \hat a = \sqrt{\frac{a_i}{2\,a_j}}\,a_{T,j}\\
        \mathcal B_{ij} &:= -2\,\left [ n_i\,n_j\,k_{T,ij} +
                           \frac1N\,(n_i^2\,n_j - n_i\,n_j^2)\,l_{T,ij}(T)
                           \right ]




    """
    @property
    def parameter_structure(self):
        result = {"T_ref": None}

        def assert_binary(name, binaries, pair):
            """Check that no binary is defined twice"""
            new_pair = sorted(pair)
            msg = f"duplicate definition of {name}: {new_pair}"
            assert new_pair not in binaries, msg
            binaries.add(sorted([s_i, s_j]))

        for name in "k_1 k_2 k_3 l_1 l_2 l_3".split():
            try:
                options = self.options[name]
            except KeyError:
                pass
            else:
                res_i = {}
                binaries = set()
                for s_i, s_j in options:
                    assert_binary(name, binaries, [s_i, s_j])
                    res_i[s_i] = res_i.get(s_i, {})
                    res_i[s_i][s_j] = None
                result[name] = res_i
        return result

    def define(self, res, par):
        pass
