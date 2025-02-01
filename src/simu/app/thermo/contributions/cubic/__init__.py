# -*- coding: utf-8 -*-

# internal modules
from simu import ThermoContribution, exp, conditional, sqrt, qpow, qsum
from simu.core.utilities import SymbolQuantity

from .rk import (RedlichKwongEOSLiquid, RedlichKwongEOSGas,
                 RedlichKwongAFunction, RedlichKwongBFunction,
                 RedlichKwongMFactor)


class LinearMixingRule(ThermoContribution):
    r"""This linear mixing rule represents any contribution that computes a
    lumped quantity as weighted sum over the molar quantities:

    .. math::

        x = \sum_i x_i\,n_i

    The contribution requires an ``option`` dictionary with the following
    entries:

      - ``target``: The name of the property :math:`x`.
      - ``source``: The name of the property :math:`x_i`. If omitted, it
        will be generated as the target name with ``_i`` appended.
    """

    def define(self, res):
        target = self.options["target"]
        source = self.options.get("source", target + "_i")
        res[target] = res[source].T @ res["n"]


class NonSymmetricMixingRule(ThermoContribution):
    r"""The mixing rule combines the pure species a-contributions ``a_i``
    (:math:`a_i`) into a lumped ``a`` property.

    The contribution requires an ``option`` dictionary with the following
    entries:

      - ``target``: The name of the property :math:`a`.
      - ``source``: The name of the property :math:`a_i`. If omitted, it
        will be generated as the target name with ``_i`` appended.

    Implemented with sparse binary interaction coefficient structure, the
    contribution includes a symmetric interaction (:math:`k`) and an
    antisymmetric contribution (:math:`l`).

    .. math::

        a = \sum_{ij} \sqrt{a_i(T)\,a_j(T)} \left [
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

    Each pair of species can only be defined once. The pre-factor 2 is used to
    yield the same parameter values as if :math:`k` was a symmetric matrix and
    :math:`l` was antisymmetric.

    Despite what above equation suggests based on the double sum, the
    complexity of this contribution in terms of both memory and runtime is
    linear in the number of species and linear in the number of non-zero
    interaction parameters. This is achieved using the relationship

    .. math::
        \sum_{ij} \sqrt{a_i(T)\,a_j(T)} n_i\,n_j
          = \left (\sum_i \sqrt{a_i(T)}\,n_i \right )^2
    """

    def define(self, res):
        target = self.options["target"]

        def extract():
            source = self.options.get("source", target + "_i")

            temp, n, a_i = res["T"], res["n"], res[source]
            tau = temp / self.par_scalar("T_ref", "K")
            tau_1, tau_2 = 1 - tau, 1 - 1 / tau
            a_n = sqrt(a_i) * n
            N = qsum(n)
            return n, N, a_n, tau_1, tau_2

        n, N, a_n, tau_1, tau_2 = extract()

        # pre-factors can be reused to minimise graph size
        cache = [{}, {}]

        def c_1(i: int, j: int) -> SymbolQuantity:
            if (i, j) not in cache[0]:
                idx_i = self.species.index(i)
                idx_j = self.species.index(j)
                cache[0][(i, j)] = a_n[idx_i] * a_n[idx_j]
            return cache[0][(i, j)]

        def c_2(i: int, j: int) -> SymbolQuantity:
            if (i, j) not in cache[1]:
                idx_i = self.species.index(i)
                idx_j = self.species.index(j)
                cache[1][(i, j)] = c_1(i, j) * (n[idx_i] - n[idx_j])
            return cache[1][(i, j)]

        def symmetric():
            contributions = []
            for name, factor in [("k_1", 1), ("k_2", tau_1), ("k_3", tau_2)]:
                try:
                    pairs = self.options[name]
                except KeyError:
                    continue
                coeff = self.par_sparse_matrix(name, pairs, "dimless")
                if coeff:
                    term = sum(c_1(i, j) * c for i, j, c in coeff.pair_items())
                    contributions.append(factor * term)
            return -2 * sum(contributions) if contributions else None

        def asymmetric():
            contributions = []
            for name, factor in [("l_1", 1), ("l_2", tau_1), ("l_3", tau_2)]:
                try:
                    pairs = self.options[name]
                except KeyError:
                    continue
                coeff = self.par_sparse_matrix(name, pairs, "dimless")
                if coeff:
                    term = sum(c_2(i, j) * c for i, j, c in coeff.pair_items())
                    contributions.append(factor * term)
            return -2 / N * sum(contributions) if contributions else None

        sym = symmetric()
        asym = asymmetric()
        res[target] = qsum(a_n) ** 2
        if sym is not None:
            res[target] += sym
        if asym is not None:
            res[target] += asym

        self.add_bound("T", res["T"])


class CriticalParameters(ThermoContribution):
    r"""This class does not perform any calculations, but provides the basic
    critical parameters as a basis for the typical equation of state
    contributions.

    The following parameters need to be provided (all as species vector):

    ======== ============== =========================
    Property Symbol         Description
    ======== ============== =========================
    T_c      :math:`T_c`    Critical temperatures [K]
    p_c      :math:`p_c`    Critical pressure[Pa]
    omega    :math:`\omega` Acentric factor [-]
    ======== ============== =========================

    The same symbols will just be published as intermediate results for the
    actual model contributions to be consumed.
    """

    provides = ["_T_c", "_p_c", "_omega"]

    def define(self, res):
        res["_T_c"] = self.par_vector("T_c", self.species, "K")
        res["_p_c"] = self.par_vector("p_c", self.species, "bar")
        res["_omega"] = self.par_vector("omega", self.species, "dimless")


class VolumeShift(ThermoContribution):
    r"""This class does not perform any calculations, but provides volume
    shift parameters to be used via mixing rules as the C-parameter in
    equations of state. The following parameter needs to be provided as a
    species vector:

    ======== ============== =================================
    Property Symbol         Description
    ======== ============== =================================
    c_i      :math:`c_i`    Volume shift parameter [m**3/mol]
    ======== ============== =================================
    """

    provides = ["_ceos_c_i"]

    def define(self, res):
        res["_ceos_c_i"] = self.par_vector("c_i", self.species, "m**3/mol")


class BostonMathiasAlphaFunction(ThermoContribution):
    r"""This contribution represents the Mathias alpha function with the
    Boston-Mathias extrapolation for supercritical temperatures.

    The following properties need to be provided upstream:

    ======== ============== ===========================================
    Property Symbol         Description
    ======== ============== ===========================================
    T        :math:`T`      Actual temperatures [K]
    T_c      :math:`T_c`    Critical temperatures [K]
    m_factor :math:`m_i`    m-factor as function of acentric factor [-]
    ======== ============== ===========================================

    Additionally, the contribution requires a polar parameter :math:`\eta_i`,
    named ``eta``. We define the root of the reduced temperature as
    :math:`\tau_i := \sqrt{T/T_{c,i}}`. Then, we define for
    :math:`\tau_i \le 1`:

    .. math::

        \alpha_i^{\frac12} =
          1 + m_i\,(1 - \tau_i) - \eta_i\,(1 - \tau_i)(0.7 - \tau_i^2)

    .. important::

        While the paper :cite:p:`Mathias1983` is currently unavailable to me, I
        must recognize that all but :cite:p:`AspenTech2001` quote the expression
        differently, namely:

        .. math::

            \alpha_i^{\frac12} =
              1 + m_i\,(1 - \tau_i) - \eta_i\,
              (1 - \boxed{\tau_i^2})(0.7 - \tau_i^2)

        Likely, there was a misprint in the original paper
        :cite:p:`Mathias1983`, as the parameter :math:`\eta_i` shall not alter
        the impact of the acentric factor :math:`\omega_i` at
        :math:`T = 0.7\, T_{c,i}`, as this point is defined to match the
        saturation pressure via

        .. math::

            \omega_i = -\log_{10}
              \frac{p^\mathrm{sat}_i(0.7\cdot T_{c,i})}{p_{c, i}}

        The reason why AspenTech could and did use the correct form is that the
        author P. Mathias was working for the company at the time of these
        developments, and the co-workers were not dependent on the journal
        publication.

    As described in Appendix (:ref:`alpha_extensions`), the extrapolation into
    the super-critical region is implemented as

    .. math::

        \alpha_i^{\frac12} = \left [\frac{c}{d}(1-\tau^{d})\right ]
        \quad\text{with}\quad c = m + 0.3\eta
        \quad\text{and}\quad d = 1 + \frac{4\,\eta}{c} + c

    The calculated vector is provided as a property called ``alpha``
    """

    provides = ["_alpha"]

    def define(self, res):
        eta = self.par_vector("eta", self.species, "dimless")
        temp, critical_temp, m_fac = res["T"], res["_T_c"], res["_m_factor"]
        tau = temp / critical_temp
        stau = sqrt(tau)

        # define sub and super-critical expression
        alpha_sub = 1 + m_fac * (1 - stau) - eta * (1 - stau) * (0.7 - tau)

        bm_c = m_fac + 0.3 * eta
        bm_d = 1 + bm_c + 4 * eta / bm_c
        alpha_sup = exp(bm_c / bm_d * (1 - qpow(stau, bm_d)))

        alpha = conditional(tau > 1, alpha_sub, alpha_sup)
        # result is square of above
        res["_alpha"] = alpha * alpha

        res["T"] = temp