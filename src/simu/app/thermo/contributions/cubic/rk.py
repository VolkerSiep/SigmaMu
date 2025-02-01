# stdlib modules
from abc import ABC, abstractmethod

# external modules
from casadi import SX
from numpy import roots

# internal modules
from simu import (
    InitialState, base_magnitude, jacobian, log, qsum, R_GAS,
    ThermoContribution, Quantity
)


class RedlichKwongEOS(ThermoContribution, ABC):
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

    .. important::

        The original publication on the volume correction of EOS
        :cite:p:`Peneloux82` suggests to shift both V and B, not only V. The
        derivation of this proposal is not given, and a numerical experiment
        shows that doing so drastically impacts the equilibrium and the density
        predictions in a negative way. We therefore stick to the approach taken
        by :cite:p:`AspenTech2001`, namely a pure volume translation. For this
        translation, it can easily be shown, as by :cite:p:`Peneloux82` that
        equilibrium conditions are not affected:

        Consider the chemical potential derived from the Helmholtz function

        .. math::
            A = A(T, V, \vec n)\quad\Rightarrow\quad
            \vec \mu = \left . \frac{\partial A}{\partial \vec n}
              \right |_{T, V}

        If we translate :math:`\tilde V = V - \vec c\cdot\vec n`, we can derive

        .. math::
            \tilde A = \tilde A(T, \tilde V, \vec n)\quad\Rightarrow\quad
             \vec{\tilde \mu} = \left . \frac{\partial \tilde A}{\partial \vec n}
              \right |_{T, \tilde V}
              = \left . \frac{\partial \tilde A}{\partial \vec n}
              \right |_{T, V} +
              \left . \frac{\partial \tilde A}{\partial \tilde V}
              \right |_{T, n}\cdot
              \left . \frac{\partial \tilde V}{\partial \vec n}
                \right |_{T}
              =  \vec \mu - p\,\vec c

        Hence at given temperature and pressure, the chemical potential only
        shifts by a constant, not affecting any equilibrium constraints.
    """

    provides = ["_VCB", "_VCB_x", "_p_x", "_p_V", "_p_V_x"]

    def define(self, res):
        ab_names = ["_ceos_a", "_ceos_b"]
        T, V, n, A, B = [res[i] for i in ["T", "V", "n"] + ab_names]

        for i in ab_names:
            res[f"{i}_T"] = jacobian(res[i], T)
            res[f"{i}_n"] = jacobian(res[i], n).T  # jacobian transposes
        A_t, B_t = [res[f"{i}_T"] for i in ab_names]
        A_n, B_n = [res[f"{i}_n"] for i in ab_names]

        if (c_name := "_ceos_c") in res:  # C is optional
            C = res[c_name]
            C_t = res[f"{c_name}_T"] = jacobian(res[c_name], T)
            C_n = res[f"{c_name}_n"] = jacobian(res[c_name], n).T  # jacobian transposes
        else:
            C = res[c_name] = Quantity(SX(1, 1), "m**3/mol") * n.units
            C_t = res[f"{c_name}_T"] = Quantity(SX(1, 1), "m**3/mol/K") * n
            C_n = res[f"{c_name}_n"] = Quantity(SX(*n.shape), "m**3/mol")

        # common terms
        N = qsum(n)
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

        self.add_bound("RK_VmBC", VmBC)  # V - B + C > 0
        self.add_bound("V", V)
        self.add_bound("dp_dV", -jacobian(res["p"], V))  # dp/dv < 0
        self.add_bound("p", res["p"])

    @staticmethod
    def find_zeros(state: InitialState, properties):
        """Given conditions and properties (A, B, C contribution of RK-EOS),
        calculate the compressibility factor roots analytically and return
        the vector quantity of volumes based on these"""
        T, p, n = state.temperature, state.pressure, state.mol_vector
        A, B, C = [properties[f"_ceos_{i}"] for i in "abc"]
        NRT = sum(n) * R_GAS * T
        alpha = float(A * p / (NRT * NRT))  # must be dimensionless
        beta = float(B * p / NRT)  # dito
        zeros = roots([1, -1, alpha - beta * (1 + beta), -alpha * beta])
        zeros = zeros[abs(zeros.imag) < 1e-7 * abs(zeros)]
        return zeros[zeros > beta].real * NRT / p - C


    @abstractmethod
    def initial_state(self, state, properties):
        """Force implementation of this method"""
        ...


class RedlichKwongEOSLiquid(RedlichKwongEOS):
    """As a subclass of
    :class:`~simu.app.thermo.contributions.cubic.rk.RedlichKwongEOS`, this
    entity specialises on describing liquid (and super-critical) phases.
    The distinct element is the initialisation.
    """

    def initial_state(self, state, properties):
        volume = min(self.find_zeros(state, properties))
        return ([base_magnitude(state.temperature), base_magnitude(volume)] +
                list(base_magnitude(state.mol_vector)))


class RedlichKwongEOSGas(RedlichKwongEOS):
    """As a subclass of
    :class:`~simu.app.thermo.contributions.cubic.rk.RedlichKwongEOS`, this
    entity specialises on describing gas (and super-critical) phases. The
    distinct elements is the initialisation."""

    def initial_state(self, state, properties):
        zeros = self.find_zeros(state, properties)
        volume = max(zeros)
        return ([base_magnitude(state.temperature), base_magnitude(volume)] +
                list(base_magnitude(state.mol_vector)))


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

    provides = ["_ceos_a_i"]

    def define(self, res):
        omega_r2 = R_GAS * R_GAS / (9 * (2**(1 / 3) - 1))
        alpha, T_c, p_c = [res[i] for i in "_alpha _T_c _p_c".split()]
        res["_ceos_a_i"] = omega_r2 * alpha * (T_c * T_c) / p_c


class RedlichKwongBFunction(ThermoContribution):
    r"""Given critical temperature ``T_c`` (:math:`T_{c,i}`) and pressure
    ``p_c`` (:math:`p_{c,i}`), this contribution calculates the
    :math:`b`-contribution ``ceos_b_i`` for the individual species. It is

    .. math::

        b_i = \Omega_b\,\frac{R\,T_{c,i}}{p_{c,i}}
        \quad\text{with}\quad
        \Omega_b = \frac13\,(2^{1/3} - 1)
    """

    provides = ["_ceos_b_i"]

    def define(self, res):
        omega_r = R_GAS * (2**(1 / 3) - 1) / 3
        T_c, p_c = [res[i] for i in "_T_c _p_c".split()]
        res["_ceos_b_i"] = omega_r * T_c / p_c


class RedlichKwongMFactor(ThermoContribution):
    r"""This contribution calculates the Redlich Kwong m-factor that is used
    in various alpha-functions. Based on provided acentric factors ``omega``
    (:math:`\omega_i`), it calculates ``m_factor`` (:math:`m_i`) as

    .. math::

        m_i = 0.48508 + (1.55171 - 0.15613\,\omega_i)\,\omega_i
    """

    provides = ["_m_factor"]

    def define(self, res):
        omega = res["_omega"]
        m = 0.48508 + (1.55171 - 0.15613 * omega) * omega
        res["_m_factor"] = Quantity(m)
