# -*- coding: utf-8 -*-

# external modules
from casadi import conditional, dot, exp, sqrt, sum1, vertcat

# internal modules
from ...utilities import iter_binary_parameters
from ..contribution import ThermoContribution


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

    category = "cubic_eos_c"
    name = "linear"

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

    name = "non_symmetric"
    category = "cubic_eos_a"
    requires = ["cubic_eos_a_function"]

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

    category = "critical_parameters"

    @property
    def parameter_structure(self):
        t_s = ThermoContribution._tensor_structure
        return t_s(["T_C", "P_C", "OMEGA"], self.species)

    def define(self, res, par):
        for name in ["T_C", "P_C", "OMEGA"]:
            res[name] = self._vector(par[name])


class BostonMathiasAlphaFunction(ThermoContribution):
    r"""This contribution represents the Mathias alpha function with the
    Boston-Mathias extrapolation for supercritical temperatures.

    The following properties need to be provided upstream:

    ======== ============== ===========================================
    Property Symbol         Description
    ======== ============== ===========================================
    T        :math:`T`      Actual temperatures [K]
    T_C      :math:`T_c`    Critical temperatures [K]
    MFAC     :math:`m_i`    m-factor as function of acentric factor [-]
    ======== ============== ===========================================

    Additionally, the contribution requires a polar parameter :math:`\eta_i`,
    named ``ETA``. We define the root of the reduced temperature as
    :math:`\tau_i := \sqrt{T/T_{c,i}}`. Then, we define for
    :math:`\tau_i \le 1`:

    .. math::

        \alpha_i^{\frac12} =
          1 + m_i\,(1 - \tau_i) - \eta_i\,(1 - \tau_i)(0.7 - \tau_i^2)

    As described in Appendix (:ref:`alpha_extensions`), the extrapolation into
    the super-critical region is implemented as

    .. math::

        \alpha_i^{\frac12} =
          \left [\frac{c}{d}(1-\tau^{d})\right ]
       \quad\text{with}\quad
    	 c = m + 0.3\eta
           \quad\text{and}\quad
         d = 1 + \frac{4\,\eta}{c} + c

    The calculated vector is provided as a property called ``ALPHA``

    .. warning::
        This alpha-function is implemented to reproduce results of
        thermodynamic models that have been parameterised in
        `Aspen Plus by AspenTech <http://aspentech.com>`_. The models differ
        however in the continuity of the second derivative (see
        :ref:`alpha_extensions`) and therefore slightly in the calculated
        properties for temperatures higher than critical temperatures of the
        involved species.
    """

    name = "Boston_Mathias"
    category = "alpha_function"
    requires = ["critical_parameters", "m_factor"]

    @property
    def parameter_structure(self):
        t_s = ThermoContribution._tensor_structure
        return t_s(["ETA"], self.species)

    def define(self, res, par):
        eta = self._vector(par["ETA"])
        T, T_c, m = res["T"], res["T_C"], res["MFAC"]
        tau = T / T_c
        stau = sqrt(tau)

        # define sub and super-critical expression
        alpha_sub = 1 + m * (1 - stau) - eta * (1 - stau) * (0.7 - tau)

        c = m + 0.3 * eta
        d = 1 + c + 4 * eta / c
        alpha_sup = exp(c / d * (1 - stau ** d))

        # implement switch with scalar conditional function
        alpha = [conditional(tau[i] > 1, [alpha_sub[i]], alpha_sup[i])
                 for i in range(len(self.species))]
        alpha = vertcat(*alpha)
        # result is square of above
        res["ALPHA"] = alpha * alpha
