# internal modules
from simu import ThermoContribution, R_GAS, log, exp


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
        T, V, n, mw = res["T"], res["V"], res["n"], res["mw"]
        rho_c = self.par_vector("_rho_c", self.species, "bar")
        t_c = self.par_vector("_t_c", self.species, "K")

        res["_tau"] = t_c / T
        res["_rho"] = mw * n / V / rho_c

        self.add_bound("T", T)  # it is divided by T and V
        self.add_bound("V", V)


class StandardStateIAPWS(ThermoContribution):
    """
    - requires _tau and mw
    - does not separate out reference state
    - requires parameters n_1 to n_8 and g_4 to g_8
    - defines standard state as ideal gas pure component standard state
      at the critical pressure of the individual species.
      This is why we need a special ideal gas contribution, to pick up each
      species at that pressure.
    """

    provides = ["mu", "S"]

    def define(self, res):
        species = self.species
        T, tau, n = res["T"], res["_tau"], res["n"]
        pn = [self.par_vector(f"n_{i:d}", species, "-") for i in range(1, 9)]
        pg = [self.par_vector(f"g_{i:d}", species, "-") for i in range(4, 9)]
        rt = R_GAS * T

        # construct chemical potential as d(NRT*phi)/dn) = RT * phi
        phi_0 = pn[0] + pn[1] * tau + pn[2] * log(tau)
        t_exp = [exp(-pg_i * tau) for pg_i in pg]
        t_log = [log(1 - t_exp_i) for t_exp_i in t_exp]
        phi_0 += [pn_i * t_log_i for pn_i, t_log_i in zip(pn[4:], t_log)]
        res["mu"] = rt * phi_0

        # construct entropy as -dmu/dT * n
        mu_tau = pn[1] + pn[2] / tau
        mu_tau += [pn_i * g_i * e_i / (1 - e_i)
                   for pn_i, g_i, e_i in zip(pn[4:], pg, t_exp)]
        res["S"] = R_GAS * (phi_0 - tau * mu_tau).T @ n

        self.declare_vector_keys("mu")
