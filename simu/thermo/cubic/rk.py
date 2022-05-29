# -*- coding: utf-8 -*-

# stdlib modules
from abc import abstractmethod

# external modules
from numpy import roots
from casadi import log, jacobian, sum1, SX

# internal modules
from ..contribution import ThermoContribution
from ...constants import R_GAS


class RedlichKwongEOS(ThermoContribution):
    r"""This contribution implements a general Redlich-Kwong equation of state
    with Peneloux volume translation:

    .. math::

        p = \frac{N\,R\,T}{V - B + C} - \frac{A}{(V + C)\,(V + B + C)}

    The following properties need to be provided by upstream contributions:

    ======== ========= ======== ==========
    Property Symbol    UOM
    ======== ========= ======== ==========
    ceos_a   :math:`A` J * |m3|
    ceos_b   :math:`B` |m3|
    ceos_c   :math:`C` |m3|     (optional)
    ======== ========= ======== ==========

    Care is to be taken when utilising a temperature-dependent :math:`C`
    contribution, as doing so can have significant effects on the calorimetric
    properties. If ``ceos_c`` is not provided, the contribution is assumed
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

    .. todo::

        The original publication on the volume correction of CEOS (Peneloux?)
        suggests to shift both V and B, not only V.
    """

    provides = ["VCB", "VCB_x", "p_x", "p_V", "p_V_x"]

    def define(self, res, par):
        abc_names = ["ceos_a", "ceos_b", "ceos_c"]
        T, V, n, A, B = [res[i] for i in ["T", "V", "n"] + abc_names[:-1]]
        C = res.get("ceos_c", SX(1, 1))  # C is optional
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

        # quantities specific for relaxation
        state = res["state"]
        res["VBC"] = VmBC
        res["VBC_x"] = jacobian(VmBC, state)
        res["p_x"] = jacobian(res["p"], state)
        res["p_V"] = jacobian(res["p"], V)
        res["p_V_x"] = jacobian(res["p_V"], state)


    def relax(self, current_result, delta_state):
        # V - B + C > 0 ?
        y = current_result["VBC"]
        d_y = current_result["VBC_x"].T @ delta_state
        beta = -y / d_y if d_y < 0 else 100

        # V > 0 ?
        if delta_state[1] < 0:
            beta = min(beta, -current_result["V"] / delta_state[1])

        #  dp/dV < 0 ?
        y = current_result["p_V"]
        d_y = current_result["p_V_x"].T @ delta_state
        if d_y > 0:
            beta = min(beta, -y / d_y)

        # p > 0 for liquid phase? (cannot happen in gas phase)
        beta_more = self.relax_more(current_result, delta_state)
        beta = beta if beta_more is None else min(beta, beta_more)
        return beta

    @abstractmethod
    def relax_more(self, current_result, delta_state):
        """For the Redlich Kwong equation of state, both phases are required
        to keep the EOS denominators positive and the volume derivative of
        pressure negative. More specific constraints are to be implemented
        by the sub-classes to represent in particular the liquid phase"""

    @staticmethod
    def find_zeros(T, p, n, properties):
        A, B, C = [float(properties[f"ceos_{i}"]) for i in "abc"]
        NRT = sum(n) * R_GAS * T
        alpha = A * p / (NRT * NRT)
        beta = B * p / NRT
        zeros = roots([1, -1, alpha - beta * (1 + beta), -alpha * beta])
        zeros = zeros[abs(zeros.imag) < 1e-7 * abs(zeros)]
        return zeros[zeros > beta].real * NRT / p - C


class RedlichKwongEOSLiquid(RedlichKwongEOS):
    """As a sub-class of :class:`RedlichKwongEOS`, this entity specialises
    on describing liquid (and super-critical) phases. The distinct elements
    are the initialisation and the relaxation, ensuring the state to be within
    the correct domain with respect to volume."""

    def relax_more(self, current_result, delta_state):
        """Additionally to the main relaxation method, this implementation
        also assures positive pressures"""
        y = current_result["p"]
        d_y = current_result["p_x"].T @ delta_state
        return -y / d_y if d_y < 0 else None

    def initial_state(self, temperature, pressure, quantities, properties):
        V = min(self.find_zeros(temperature, pressure, quantities,
                                properties))
        return [temperature, V] + quantities

class RedlichKwongEOSGas(RedlichKwongEOS):
    """As a sub-class of :class:`RedlichKwongEOS`, this entity specialises
    on describing gas (and super-critical) phases. The distinct elements
    are the initialisation and the relaxation, ensuring the state to be within
    the correct domain with respect to volume."""

    def relax_more(self, current_result, delta_state):
        """For the gas phase, no more constraints apply then the ones that
        apply for both phases."""
        pass

    def initial_state(self, temperature, pressure, quantities, properties):
        V = max(self.find_zeros(temperature, pressure, quantities,
                                properties))
        return [temperature, V] + quantities


class RedlichKwongAFunction(ThermoContribution):
    r"""Given critical temperature ``T_c`` (:math:`T_{c,i}`) and pressure
    ``p_c`` (:math:`p_{c,i}`), this contribution scales the
    :math:`\alpha`-function ``alpha`` (:math:`\alpha_i`) to define the
    :math:`a`-contribution ``rk_a_i`` (:math:`a_i`) for the individual species.
    It is

    .. math::

        a_i = \alpha_i\,\Omega_a\,\frac{R^2\,T_{c,i}^2}{p_{c,i}}
        \quad\text{with}\quad
        \Omega_a = \frac19\,(2^{1/3} - 1)^{-1}
    """

    provides = ["ceos_a_i"]

    def define(self, res, par):
        omega_r2 = R_GAS * R_GAS / (9 * (2 ** (1 / 3) - 1))
        alpha, T_c, p_c = [res[i] for i in "alpha T_c p_c".split()]
        res["ceos_a_i"] = omega_r2 * alpha * (T_c * T_c) / p_c


class RedlichKwongBFunction(ThermoContribution):
    r"""Given critical temperature ``T_c`` (:math:`T_{c,i}`) and pressure
    ``p_c`` (:math:`p_{c,i}`), this contribution calculates the
    :math:`b`-contribution ``ceos_b_i`` for the individual species. It is

    .. math::

        b_i = \Omega_b\,\frac{R\,T_{c,i}}{p_{c,i}}
        \quad\text{with}\quad
        \Omega_b = \frac13\,(2^{1/3} - 1)
    """

    provides = ["ceos_b_i"]

    def define(self, res, par):
        omega_r = R_GAS * (2 ** (1 / 3) - 1) / 3
        T_c, p_c = [res[i] for i in "T_c p_c".split()]
        res["ceos_b_i"] = omega_r * T_c / p_c


class RedlichKwongMFactor(ThermoContribution):
    r"""This contribution calculates the Redlich Kwong m-factor that is used
    in various alpha-functions. Based on provided acentric factors ``omega``
    (:math:`\omega_i`), it calculates ``m_factor`` (:math:`m_i`) as

    .. math::

        m_i = 0.48508 + (1.55171 - 0.15613\,\omega_i)\,\omega_i
    """

    provides = ["m_factor"]

    def define(self, res, par):
        omega = res["omega"]
        res["m_factor"] = 0.48508 + (1.55171 - 0.15613 * omega) * omega
