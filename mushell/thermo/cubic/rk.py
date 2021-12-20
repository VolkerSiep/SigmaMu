# -*- coding: utf-8 -*-

# external modules
from casadi import log, jacobian, sum1

# internal modules
from ..contribution import ThermoContribution
from ...constants import R_GAS

class RedlichKwongEOS(ThermoContribution):
    r"""This contribution implements a general Redlich-Kwong equation of state
    with Peneloux volume translation:

    .. math::

        p = \frac{N\,R\,T}{V - B + C} - \frac{A}{(V + C)\,(V + B + C)}

    The following properties need to be provided by upstream contributions:

    ======== ========= ====== ==========
    Property Symbol    UOM
    ======== ========= ====== ==========
    RK_A     :math:`A` J * m3
    RK_B     :math:`B` m3
    CEOS_C   :math:`C` m3     (optional)
    ======== ========= ====== ==========

    Care is to be taken when utilising a temperature-dependent :math:`C`
    contribution, as doing so can have significant effects on the calorimetric
    properties. If ``CEOS_C`` is not provided, the contribution is assumed
    zero.

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
            \ln\frac{V}{V - B + C} + T\,\frac{B_T - C_T}{V - B + C}
            \right ]\\& +
            \frac1B\left (A_T - \frac{A}{B}\,B_T\right )\,
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
            \ln\frac{V}{V - B + C} + N\,\frac{B_{n,i} - C_{n,i}}{V - B + C}
            \right ]\\& +
            \frac1B\left (A_{n,i} - \frac{A}{B}\,B_{n,i}\right )\,
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

    category = "cubic_eos"
    name = "Redlich_Kwong"
    requires = ["cubic_eos_a", "cubic_eos_b",
                ["ideal_gas", "Helmholtz"],
                ["a_function", "Redlich_Kwong"],
                ["b_function", "Redlich_Kwong"]]

    def __init__(self, *args, **kwargs):
        """The constructor of this class is disabled, as an instance would
        not have a decent initialiation nor relaxation functionality and
        would therefore be without value. Instantiate objects of
        :class:`RedlichKwongEOSLiquid` or :class:`RedlichKwongEOSLiquid`
        instead."""
        raise NotImplementedError("Don't do this!")

    def define(self, res, par):
        abc_names = ["RK_A", "RK_B", "CEOS_C"]
        T, V, n, A, B = [res[i] for i in ["T", "V", "n"] + abc_names[:-1]]
        C = res.get("CEOS_C", 0)  # optional
        for i in abc_names:
            res[f"{i}_T"] = jacobian(res[i], T)
            res[f"{i}_n"] = jacobian(res[i], n).T  # jacobian transposes
        A_t, B_t, C_t = [res[f"{i}_T"] for i in abc_names]
        A_n, B_n, C_n = [res[f"{i}_n"] for i in abc_names]

        # common terms
        N = sum1(n)
        NR, RT = N * R_GAS, T * R_GAS
        VC = V + C
        VmBC, VpBC = VC - B, VC + B
        AB = A / B

        # entropy contribution
        m_dS = NR * (log(V / VmBC) + T * (B_t - C_t) / VmBC)
        m_dS += (A_t - AB * B_t) / B * log(VC / VpBC)
        m_dS += AB * (C_t / VC - (B_t + C_t) / VpBC)
        res["S"] -= m_dS

        # pressure contribution
        res["p"] -= NR * T * (1 / V - 1 / VmBC) + A / (VC * VpBC)

        # chemical potential contribution
        dmu = RT * (log(V / VmBC) + N * (B_n - C_n) / VmBC)
        dmu += (A_n - AB * B_n) / B * log(VC / VpBC)
        dmu += AB * (C_n / VC - (B_n + C_n) / VpBC)
        res["mu"] += dmu


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


class RedlichKwongAFunction(ThermoContribution):
    r"""Given critical temperature :math:`T_{c,i}` and pressure
    :math:`T_{c,i}`, this contribution scales the :math:`\alpha`-function
    (``ALPHA_I``) to define the :math:`a`-contribution (``RK_A_I``) for the
    individual species. It is

    .. math::

        a_i = \alpha_i\,\Omega_a\,\frac{R^2\,T_{c,i}^2}{p_{c,i}}
        \quad\text{with}\quad
        \Omega_a = \frac19\,(2^{1/3} - 1)^{-1}
    """

    category = "a_function"
    name = "Redlich_Kwong"
    requires = ["critical_parameters", "alpha_function"]

    def define(self, res, par):
        omega_r2 = R_GAS * R_GAS / (9 * (2 ** (1 / 3) - 1))
        alpha, T_c, p_c = [res[i] for i in "ALPHA_I T_C P_C".split()]
        res["RK_A_I"] = omega_r2 * alpha * (T_c * T_c) / p_c

class RedlichKwongBFunction(ThermoContribution):
    r"""Given critical temperature :math:`T_{c,i}` and pressure
    :math:`T_{c,i}`, this contribution calculates the
    :math:`b`-contribution (``RK_B_I``) for the individual species. It is

    .. math::

        b_i = \Omega_b\,\frac{R\,T_{c,i}}{p_{c,i}}
        \quad\text{with}\quad
        \Omega_b = \frac13\,(2^{1/3} - 1)
    """

    category = "b_function"
    name = "Redlich_Kwong"
    requires = ["critical_parameters"]

    def define(self, res, par):
        omega_r = R_GAS * (2 ** (1 / 3) - 1) / 3
        T_c, p_c = [res[i] for i in "T_C P_C".split()]
        res["RK_B_I"] = omega_r * T_c / p_c
