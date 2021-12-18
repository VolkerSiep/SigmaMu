
# -*- coding: utf-8 -*-

# external modules
from casadi import dot, jacobian, sqrt, sum1

# internal modules
from ..utilities import iter_binary_parameters
from ..constants import R_GAS
from .contribution import ThermoContribution


class RedlichKwongEOS(ThermoContribution):
    r"""This contribution implements a general Redlich-Kwong equation of state
    with Peneloux volume translation:

    .. math::

        p = \frac{N\,R\,T}{V - B + C} - \frac{A}{(V + C)\,(V + B + C)}

    The following properties need to be provided by upstream contributions:

    ======== ========= ======
    Property Symbol    UOM
    ======== ========= ======
    RK_A     :math:`A` J * m3
    RK_B     :math:`B` m3
    CEOS_C   :math:`C` m3
    ======== ========= ======

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

    The explicitly implemented temperature derivative is

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

    def define(self, res, par):
        abc_names = ["RK_A", "RK_B", "CEOS_C"]
        T, n, A, B, C = [res[i] for i in ["T", "n"] + abc_names]
        for i in abc_names:
            res[f"{i}_T"] = jacobian(res[i], T)
            res[f"{i}_n"] = jacobian(res[i], n)
        A_t, B_t, C_t = [res[f"{i}_T"] for i in abc_names]
        A_n, B_n, C_n = [res[f"{i}_n"] for i in abc_names]

        # TODO:
        #  - TEST THIS AND ALL OTHER CUBIC CONTRIBUTIONS SO FAR!!!!
        #  - update S, p and mu

        pass


class RedlichKwongEOSLiquid(RedlichKwongEOS):
    """As a sub-class of :class:`RedlichKwongEOS`, this entity specialises
    on describing liquid (and super-critical) phases. The distinct elements
    are the initialisation and the relaxation, ensuring the state to be within
    the correct domain with respect to volume."""
    def __init__(self, species, options):
        # reenable constructor
        ThermoContribution.__init__(self, species, options)

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
    r"""The Peneloux volume shift provides the ``CEOS_C`` parameter. As such,
    this contribution determines the `C` parameter in cubic equations of state,
    representing the volume translation.

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
        res["CEOS_C"] = dot(self._vector(par["c_i"]), res["n"])


class NonSymmmetricMixingRule(ThermoContribution):
    r"""The mixing rule combines the pure species a-contributions ``RK_A_I``
    (:math:`a_i`) into the ``RK_A`` parameter of the :class:`RedlichKwongEOS`
    class.

    Implemented with sparse binary interaction coefficient structure, the
    contribution includes a symmetric interaction (:math:`k`) and an
    anti-symmetric contribution (:math:`l`).

    .. math::

        A = \sum_{ij} \sqrt{a_i(T)\,a_j(T)} \left [
              n_i\,n_j - 2\,n_i\,n_j\,k_{ij}(T) -
              \frac2N\,(n_i^2\,n_j - n_i\,n_j^2)\,l_{ij}(T)
            \right ]

    These dimensionless interaction coefficients can be a function of
    temperature according to

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

    Despite what above equation suggests based on the double sum, the
    complexity of this contribution in terms of both memory and runtime is
    linear in the number of species and linear in the number of non-zero
    interaction parameters. This is achieved using the relationship

    .. math::
        \sum_{ij} \sqrt{a_i(T)\,a_j(T)} n_i\,n_j
          = \left (\sum_i \sqrt{a_i(T)}\,n_i \right )^2
    """
    @property
    def parameter_structure(self):
        result = {"T_ref": None}

        def assert_binary(name, binaries, pair):
            """Check that no binary is defined twice"""
            new_pair = tuple(sorted(pair))
            msg = f"duplicate definition of {name}: {new_pair}"
            assert new_pair not in binaries, msg
            binaries.add(new_pair)

        for name in "k_1 k_2 k_3 l_1 l_2 l_3".split():
            try:
                options = self.options[name]
            except KeyError:
                pass
            else:
                res_i = {}
                binaries = set()
                for s_i, s_j in options:
                    assert_binary(name, binaries, (s_i, s_j))
                    res_i[s_i] = res_i.get(s_i, {})
                    res_i[s_i][s_j] = None
                result[name] = res_i
        return result


    def define(self, res, par):
        T, n, a_i = res["T"], res["n"], res["RK_A_I"]
        tau = T / par["T_ref"]
        tau_1, tau_2 = 1 - tau, 1 - 1 / tau

        # calculate first term
        an = sqrt(a_i) * n
        result = sum1(an)
        result *= result  # = sum_ij an_i * an_j

        # calculate second term (symmetric interaction)
        # pre-factors can be reused to minimise graph size
        cache = {}
        for name, factor in [("k_1", 1), ("k_2", tau_1), ("k_3", tau_2)]:
            coefficients = iter_binary_parameters(self.species, par, name)
            for i, j, symbol in coefficients:
                if (i, j) not in cache:
                    cache[(i, j)] = 2* an[i] * an[j]
                result -= cache[(i, j)] * symbol * factor
        # calculate second term (symmetric interaction)
        two_by_n = 2.0 / sum1(n)
        cache = {}
        for name, factor in [("l_1", 1), ("l_2", tau_1), ("l_3", tau_2)]:
            coefficients = iter_binary_parameters(self.species, par, name)
            for i, j, symbol in coefficients:
                if (i, j) not in cache:
                    cache[(i, j)] = two_by_n * an[i] * an[j] * (n[i] - n[j])
                result -= cache[(i, j)] * symbol * factor
        res["RK_A"] = result


# TODO:
#   - implement a_i contribution
#   - implement alpha function
#   - implement m-factor
#   - implement b_i contribution


class CriticalParameters(ThermoContribution):
    r"""This class does not perform any calculations, but provides the basic
    critical parameters as a basis for the typical equation of state
    contributions.

    The following parameters need to be provided (all as species vector):

    ======== ============== =========================
    Property Symbol         Description
    ======== ============== =========================
    T_C      :math:`T_c`    Critical temperatures [K]
    P_C      :math:`p_c`    Critical pressure[Pa]
    OMEGA    :math:`\omega` Acentric factor [-]
    ======== ============== =========================

    The same symbols will just be published as intermediate results for the
    actual model contributions to be consumed.
    """
    @property
    def parameter_structure(self):
        t_s = ThermoContribution._tensor_structure
        return t_s(["T_C", "P_C", "OMEGA"], self.species)

    def define(self, res, par):
        for name in ["T_C", "P_C", "OMEGA"]:
            res[name] = self._vector(par[name])
