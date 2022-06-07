# -*- coding: utf-8 -*-

# external modules
from casadi import conditional, dot, exp, sqrt, sum1, vertcat

# internal modules
from ...utilities import iter_binary_parameters
from ..contribution import ThermoContribution

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
      - ``src_mode``: If set to ``child`` (default), the contribution aquires
         the quantities :math:`x_i` as a result of previous calculations.
         If set to ``parameter``, the quantity is required as a vector
         parameter to the contribution.
    """

    CHILD = "CHILD"  # canot be enum as I wish to serialise readable with json
    PARAMETER = "PARAMETER"

    @property
    def parameter_structure(self):
        if self.options.get("src_mode", "child") == "child":
            return {}
        else:
            cts = ThermoContribution.create_tensor_structure
            name = self.options.get("source", self.options["target"] + "_i")
            return {name: cts(self.species)}

    def define(self, res, par):
        target = self.options["target"]
        source = self.options.get("source", target + "_i")
        src_mode = self.options.get("src_mode", self.CHILD).upper()
        assert src_mode in (self.CHILD, self.PARAMETER), \
            f"Invalid src_mode: '{src_mode}'"
        if src_mode == self.CHILD:
            res[target] = dot(res[source], res["n"])
        else:
            res[target] = dot(self.create_vector(par[source]), res["n"])


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
    anti-symmetric contribution (:math:`l`).

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
                pass  # all interaction matrices are optional
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
        target = self.options["target"]
        source = self.options.get("source", target + "_i")

        temp, quant, a_i = res["T"], res["n"], res[source]
        tau = temp / par["T_ref"]
        tau_1, tau_2 = 1 - tau, 1 - 1 / tau

        # calculate first term
        a_n = sqrt(a_i) * quant
        result = sum1(a_n)
        result *= result  # = sum_ij an_i * an_j

        # calculate second term (symmetric interaction)
        # pre-factors can be reused to minimise graph size
        cache = {}
        for name, factor in [("k_1", 1), ("k_2", tau_1), ("k_3", tau_2)]:
            coefficients = iter_binary_parameters(self.species, par, name)
            for i, j, symbol in coefficients:
                if (i, j) not in cache:
                    cache[(i, j)] = 2 * a_n[i] * a_n[j]
                result -= cache[(i, j)] * symbol * factor
        # calculate second term (symmetric interaction)
        two_by_n = 2.0 / sum1(quant)
        cache = {}
        for name, factor in [("l_1", 1), ("l_2", tau_1), ("l_3", tau_2)]:
            coefficients = iter_binary_parameters(self.species, par, name)
            for i, j, symbol in coefficients:
                if (i, j) not in cache:
                    cache[(i, j)] = two_by_n * a_n[i] * a_n[j] * (quant[i] - quant[j])
                result -= cache[(i, j)] * symbol * factor
        res[target] = result


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

    provides = ["T_c", "p_c", "omega"]

    @property
    def parameter_structure(self):
        cts = ThermoContribution.create_tensor_structure
        return cts(["T_c", "p_c", "omega"], self.species)

    def define(self, res, par):
        for name in ["T_c", "p_c", "omega"]:
            res[name] = self.create_vector(par[name])


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

    As described in Appendix (:ref:`alpha_extensions`), the extrapolation into
    the super-critical region is implemented as

    .. math::

        \alpha_i^{\frac12} = \left [\frac{c}{d}(1-\tau^{d})\right ]
        \quad\text{with}\quad c = m + 0.3\eta
        \quad\text{and}\quad d = 1 + \frac{4\,\eta}{c} + c

    The calculated vector is provided as a property called ``alpha``

    .. warning::
        This alpha-function is implemented to reproduce results of
        thermodynamic models that have been parameterised in
        `Aspen Plus by AspenTech <http://aspentech.com>`_. The models differ
        however in the continuity of the second derivative (see
        :ref:`alpha_extensions`) and therefore slightly in the calculated
        properties for temperatures higher than critical temperatures of the
        involved species.
    """

    provides = ["alpha"]

    @property
    def parameter_structure(self):
        cts = ThermoContribution.create_tensor_structure
        return cts(["eta"], self.species)

    def define(self, res, par):
        eta = self.create_vector(par["eta"])
        temp, critical_temp, m_fac = res["T"], res["T_c"], res["m_factor"]
        tau = temp / critical_temp
        stau = sqrt(tau)

        # define sub and super-critical expression
        alpha_sub = 1 + m_fac * (1 - stau) - eta * (1 - stau) * (0.7 - tau)

        bm_c = m_fac + 0.3 * eta
        bm_d = 1 + bm_c + 4 * eta / bm_c
        alpha_sup = exp(bm_c / bm_d * (1 - stau ** bm_d))

        # implement switch with scalar conditional function
        alpha = [conditional(tau[i] > 1, [alpha_sub[i]], alpha_sup[i])
                 for i in range(len(self.species))]
        alpha = vertcat(*alpha)
        # result is square of above
        res["alpha"] = alpha * alpha
