from abc import abstractmethod

from casadi import DM
from simu import (ThermoContribution, registered_contribution, Quantity,
                  N_A, E_0, EPS_0, K_B, R_GAS, PI, sqrt, log, qvertcat)

_M0 = Quantity(1.0, "mol/kg")


@registered_contribution
class ElectrolyteBasics(ThermoContribution):
    r"""

    """
    def define(self, res):
        n = res["n"]

        # find index of water
        species_def = self.species_definitions
        solvent_name = self.options.get("solvent_name", "H2O")
        solvent_idx = self.species.index(solvent_name)

        # unity vector in solvent direction
        delta = DM.zeros(len(self.species))
        delta[solvent_idx] = 1
        res["_delta_i_s"] = Quantity(delta)

        res["_mw_solvent"] = species_def[solvent_name].molecular_weight
        m_sol = n[solvent_idx] * res["_mw_solvent"]

        res["m_solvent"] = m_sol
        res["charge"] = c = qvertcat(*[s.charge for s in species_def.values()])
        res["molality"] = m = n / (m_sol * _M0)
        res["I"] = m.T @ (c ** 2) / 2


class ExcessBasePitzer(ThermoContribution):
    r"""The Pitzer model is formulated in terms of a reduced excess Gibbs
    energy as follows:

    .. math::

        \frac{G^\mathrm{ex}}{R\,T\,M_s\,n_s\,m_0} = \chi(T, m_i)

    Here, :math:`M_s` is the molecular weight of the solvent, :math:`n_s` the
    molar quantity or flow of the solvent, :math:`m_0 = 1 mol/kg`, and
    :math:`m_i = n_i/(M_s\,n_s\,m_0)` the dimensionless molalities.

    Both the long-range contribution (Pitzer-Debye-Hückel) and the short-range
    contribution expressed by binary and ternary parameters are part of the
    function :math:`\chi(T, m_i)`.

    The chemical potential is then

    .. math::

        \Delta \mu_i = R\,T\,\left [\frac{\partial \chi}{\partial m_i} +
         \left (M_s\,n_s\,\chi - \frac{\partial \chi}{\partial m_j}\,
           \frac {n_j}{n_s} \right )\right ]

    The entropy is

    .. math::

        \Delta S = R\,M_s\,n_s\,m_0\,\left [
            \chi + T\,\frac{\partial \chi}{\partial T} \right ]
    """
    def define(self, res):
        temp, n, ionic_strength = res["T"], res["n"], res["I"]
        molality, d_is = res["molality"], res["_delta_i_s"]
        m_solvent, mw_solvent = res["m_solvent"], res["mw_solvent"]
        chi = self._get_chi(temp, molality)
        s_res = m_solvent * _M0 * R_GAS * (chi["chi"] + temp * chi["chi_T"])

        mu_res = (1 - d_is) * R_GAS * temp * chi["chi_m"]
        mu_res += (mw_solvent * _M0 * chi["chi"]) * d_is

        res["S"] += s_res
        res["mu"] += mu_res

    def _get_chi(self, temp: Quantity, molality: Quantity):
        # use chain rule on fresh independent variables, like
        # _get_phi in ResidualBaseIAPWS

        return {"chi": None,
                "chi_T": None,
                "chi_m": None}

    @staticmethod
    @abstractmethod
    def define_chi(self, temp: Quantity, ionic_strength: Quantity,
                   molality: Quantity):
        pass


@registered_contribution
class PitzerDebyeHueckel(ExcessBasePitzer):
    r"""

    """

    # TODO: into define_chi
    def define_chi(self, temp, ionic_strength, molality):
        rho_sol = self.par_scalar("RHO_SOLVENT", "kg/m**3")
        eps_r = self.par_scalar("EPS_SOLVENT", "dimless")
        b = self.par_scalar("B", "dimless")

        a_gamma = (
            sqrt(2 * PI * N_A * _M0 * rho_sol)
            * (E_0 ** 2 / (4 * PI * EPS_0 * eps_r * K_B * temp)) ** (3 / 2))
        b_sqi = b * sqrt(ionic_strength)
        f = -4 / 3  * a_gamma * ionic_strength * log(1 + b_sqi) / b

        # TODO: can I still have a mechanism to save f in res?
        return f
