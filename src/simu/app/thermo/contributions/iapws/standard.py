# internal modules
from simu import ThermoContribution, R_GAS, log, exp, qsum


class ReducedStateIAPWS(ThermoContribution):
    r"""In the IAPWS EOS, all contributions are expressed in terms of

    .. math::

        \tau_i = \frac{T_{c,i}}{T} \quad\text{and}\quad
        \varrho_i = \frac{\rho_i}{\rho_{c,i}}
          = \frac{M_i\,n_i}{V\,\rho_{c,i}}

    There is no limitation at this point to pure species models. The
    contribution simply calculates the dimensionless temperatures and
    mole-densities for all species in vector form.

    This contribution is part of the EOS as described in :cite:p:`Wagner_2002`,
    section 6.2.

    It relies on prior definition of molecular weight ``mw``, as the volumetric
    terms are implemented as density. The following parameters need to be
    provided as species vectors:

    ======== ============== =========================
    Property Symbol         Description
    ======== ============== =========================
    T_c      :math:`T_c`    Critical temperatures [K]
    rho_c    :math:`\rho_c` Critical density [kg/m3]
    ======== ============== =========================

    .. note::

        Despite numerous occasions in literature where :math:`\tau := T / T_c`,
        the reciprocal value is used here.

    """
    provides = ["_tau", "_rho"]

    def define(self, res):
        temp, vol, n, mw = res["T"], res["V"], res["n"], res["mw"]
        rho_c = self.par_vector("rho_c", self.species, "kg/m^3")
        t_c = self.par_vector("T_c", self.species, "K")

        res["_tau"] = t_c / temp
        res["_rho"] = mw * n / vol / rho_c

        self.add_bound("T", temp)  # it is divided by T and V
        self.add_bound("V", vol)


class StandardStateIAPWS(ThermoContribution):
    r"""
    With previously defined :math:`\tau` (``_tau``) and  :math:`\varrho`
    (``_rho``), the standard state of the IAPWS model according to
    :cite:`Wagner_2002` is defined via the chemical potential. With

    .. math::

      \eta_i = \exp(-\gamma_{k,i}^0\,\tau_i)

    we write

    .. math::

        \phi_i = n_{1, i}^0 + n_{2, i}^0\,\tau_i + n_{3, i}^0\,\ln \tau_i
          + \sum_{k=4}^{8}\,n_{k, i}^0\,\ln(1 - \eta_i)

    As such, the chemical potential is

    .. math:: \mu_i = R\,T\,\phi_i.

    The entropy contribution is calculated as

    .. math::

        S = \sum_i n_i\, \left .
              \frac{\partial \mu_i}{\partial T} \right |_p
          = R\,\sum_i n_i\, \left [
              \phi_i - \tau_i\, \left .
                \frac{\partial \phi_i}{\partial \tau_i} \right |_p
            \right ]

    The derivative is implemented as

    .. math::

        \tau_i\, \left . \frac{\partial \phi_i}{\partial \tau_i} \right |_p
        = n_{2, i}^0\,\tau_i + n_{3, i}^0 +
          \sum_{k=4}^{8}\,n_{k, i}^0\,\gamma_{k, i}^0\,\tau_i\,
          \frac{\eta_i}{1-\eta_i}


    As such, the contribution requires the following standard state parameters

    ========= ============================== ========= ======================
    Parameter Symbol                         Parameter Symbol
    ========= ============================== ========= ======================
    ``n_1``   :math:`n_{1,i}^0\hspace{3cm}`  ``g_4``   :math:`\gamma_{4,i}^0`
    ...                                      ...
    ``n_8``   :math:`n_{8,i}^0`              ``g_8``   :math:`\gamma_{8,i}^0`
    ========= ============================== ========= ======================

    The standard state is technically an ideal gas pure component standard state
    at the critical point of the individual species.
    """

    provides = ["mu", "S"]

    def define(self, res):
        species = self.species
        T, tau, n = res["T"], res["_tau"], res["n"]
        pn = [self.par_vector(f"n_{i:d}", species, "") for i in range(1, 9)]
        pg = [self.par_vector(f"g_{i:d}", species, "") for i in range(4, 9)]
        rt = R_GAS * T

        # construct chemical potential as d(NRT*phi)/dn) = RT * phi
        phi_0 = pn[0] + pn[1] * tau + pn[2] * log(tau)
        t_exp = [exp(-pg_i * tau) for pg_i in pg]
        t_log = [log(1 - t_exp_i) for t_exp_i in t_exp]
        phi_0 += sum(pn_i * t_log_i for pn_i, t_log_i in zip(pn[3:], t_log))
        res["mu"] = rt * phi_0

        # construct entropy as -dmu/dT * n
        phi_tau = pn[1] * tau + pn[2]
        phi_tau += sum(pn_i * g_i * tau * e_i / (1 - e_i)
                       for pn_i, g_i, e_i in zip(pn[3:], pg, t_exp))
        res["S"] = -R_GAS * (phi_0 - phi_tau).T @ n

        self.declare_vector_keys("mu")


class IdealGasIAPWS(ThermoContribution):
    r"""The standard state :class:`StandardStateIAPWS` is defined with a
    reference to the critical density of the individual species. This
    contribution adds the ideal gas and ideal mix contribution based on this
    reference. Due to the individual nature, it is not efficient to distinguish
    between the ideal mix and ideal gas contribution.

    The IAPWS ideal gas contribution in terms of Helmholtz energy is

    .. math::

        \Delta A^\mathrm{ig} =
        \sum_k n_k\,R\,T\,\ln \varrho_k
        = \sum_k n_k\,R\,T\,\ln \frac{M_k\,n_j}{V\,\rho_{c, k}}

    This yields a chemical potential update

    .. math::

        \mu_i \leftarrow \mu_i +
          R\,T\,\left (\ln \frac{M_i\,n_i}{V\,\rho_{c, i}} + 1 \right )
          = \mu_i + R\,T\,(\ln \varrho_i + 1)

    And for the entropy:

    .. math::

        S \leftarrow S -
          \sum_k n_k\,R\,\ln \varrho_k

    Finally, pressure is defined here as

    .. math:: p = \frac{\sum_k n_k\,R\,T}{V}

    As such, this contribution does not require further parameters.
    """

    provides = ["p"]

    def define(self, res):
        temperature, n, rho, volume = res["T"], res["n"], res["_rho"], res["V"]
        ln_rho = log(rho)
        n_total = qsum(n)
        r_t = R_GAS * temperature
        n_r_t = n_total * r_t
        res["mu"] += r_t * (ln_rho + 1)
        res["S"] -= R_GAS * (n.T @ ln_rho)
        res["p"] = n_r_t / volume