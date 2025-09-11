from abc import abstractmethod

from casadi import DM
from simu import (ThermoContribution, registered_contribution, Quantity,
                  N_A, E_0, EPS_0, K_B, R_GAS, PI, sqrt, log, exp, qvertcat)
from simu.core.utilities.types import Map, MutMap

_M0 = Quantity(1.0, "mol/kg")
_C0 = Quantity(1.0, "e/mol")


@registered_contribution
class ElectrolyteBasics(ThermoContribution):
    r"""This contribution prepares some basic properties relevant for
    electrolyte systems.

    An input option ``solvent_name`` can be used to specify the name of the
    solvent species, which is by default ``H2O``.

    Charges and molalities are considered as dimensionless by normalizing with
    :math:`m_0 = 1\ \mathrm{mol/kg}` and :math:`c_0 = 1\ \mathrm{e / mol}`.

    ``_delta_i_s`` (:math:`\delta_{is}`)
        The Kronecker operator, being unity if :math:`i = s` and otherwise zero

    ``_mw_solvent`` (:math:`M_s`)
        The molecular weight of the solvent [kg/mol]

    ``_m_solvent`` (:math:`\hat m_s`)
        The mass (flow) of the solvent [kg]/[kg/s]

    ``charge`` (:math:`c_i`)
        The charge vector [e/mol]

    ``molality`` (:math:`m_i`)
       The molality vector: :math:`m_i = n_i / (\hat m_s\cdot m_0)` [-]

    ``I`` (:math:`I`)
       Ionic strength: :math:`I = \sum_i m_i\,(c_i / c_0)^2` [-]

    """

    provides = ["_delta_i_s", "_mw_solvent", "m_solvent",
                "charge", "molality", "I"]

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
        res["charge"] = qvertcat(*[s.charge for s in species_def.values()])
        res["molality"] = m = n / (m_sol * _M0)
        res["I"] = m.T @ ((res["charge"] / _C0) ** 2) / 2


class ExcessBasePitzer(ThermoContribution):
    r"""The Pitzer model is formulated in terms of a reduced excess Gibbs
    energy as follows:

    .. math::

        \frac{G^\mathrm{ex}}{R\,T\,M_s\,n_s\,m_0} = \chi(T, m_i)

    Here, :math:`M_s` is the molecular weight of the solvent, :math:`n_s` the
    molar quantity or flow of the solvent, :math:`m_0 = 1\ \mathrm{mol/kg}`,
    and :math:`m_i = n_i/(M_s\,n_s\,m_0)` the dimensionless molalities.

    Both the long-range contribution (Pitzer-Debye-Hückel) and the short-range
    contribution expressed by binary and ternary parameters are part of the
    function :math:`\chi(T, m_i)`.

    The chemical potential is then

    .. math::

        \Delta \mu_i = R\,T\,\left [\frac{\partial \chi}{\partial m_i} +
         \left (M_s\,n_s\,\chi - \frac{\partial \chi}{\partial m_j}\,
           \frac {n_j}{n_s} \right )\delta_{is}\right ]

    Here, :math:`\delta_{is}` is the Kronecker operator, being unity if
    :math:`i = s` and otherwise zero. The entropy is

    .. math::

        \Delta S = R\,M_s\,n_s\,m_0\,\left [
            \chi + T\,\frac{\partial \chi}{\partial T} \right ]
    """
    def define(self, res):
        temp, n, ionic_strength = res["T"], res["n"], res["I"]
        molality, d_is = res["molality"], res["_delta_i_s"]
        charge = res["charge"] / _C0
        m_solvent, mw_solvent = res["m_solvent"], res["mw_solvent"]
        chi_res = self.define_chi(res)
        chi, chi_t = chi_res["chi"], chi_res["chi_t"]
        chi_m = chi_res["chi_m"] + chi_res["chi_i"] * charge ** 2 / 2

        s_res = m_solvent * _M0 * R_GAS * (chi + temp * chi_t)
        mu_res = (1 - d_is) * R_GAS * temp * chi_m
        mu_res += (mw_solvent * _M0 * chi) * d_is

        res["S"] += s_res
        res["mu"] += mu_res

    @abstractmethod
    def define_chi(self, res: MutMap[Quantity]) -> Map[Quantity]:
        """
        Provide dimensionless residual contribution :math:`\chi(T, I, m_i)`
        and the partial derivatives :math:`\chi_T` (``chi_t``,
        :math:`\chi_I` (``chi_i``) and :math:`\chi_{m_i}` (``chi_m``).
        """
        ...


@registered_contribution
class PitzerDebyeHueckel(ExcessBasePitzer):
    r"""The Pitzer-Debye-Hückel rule defines the long-range interaction as
    described in :cite:`Pitzer_1980`.

    A Debye-Hückel parameter is solely a function of temperature and defined as

    .. math::

       A_\gamma =\sqrt{2\pi\, N_A\, m_0\,\rho_{\rm s}}\cdot
         \left [\frac{e^2}{4\pi\,\epsilon_0\epsilon_r\,k_B\,T} \right ]^{3/2}

    Here, :math:`N_A` is the Avogadro constant, :math:`e` the electron charge,
    :math:`\epsilon_0` the vacuum permittivity, and :math:`k_B` the Boltzmann
    constant.
    Further, :math:`\rho_s` the density and :math:`\epsilon_r` the
    dielectric constant of the solvent.

    One could consider the solvent properties as being temperature-dependent,
    but this is omitted, as the dominant impact is captured by the
    :math:`T^{-3/2}` dependency, and the model was developed with the intention
    to keep the property values at 25 |degC|. Secondary impact of temperature
    is to be caught by parametrising the short-range binary and ternary
    interaction.

    The dimensionless long-range contribution is then additionally a function
    of ionic strength :math:`I`.

    .. math::

       f(I)=-\frac43A_\gamma I^{3/2}\,\frac{\ln (1+b\sqrt{I})}{b\sqrt{I}}
            \quad {\rm with}\quad b=1.2

    The function :math:`f` is a compatible contribution :math:`\chi` as defined
    for the :class:`ExcessBasePitzer` base-class.

    The required derivatives are

    .. math::
        f_T = -\frac32\,\frac{f}{T}\qquad
        f_I = \frac{f}{I} -
          \frac23\,A_\gamma\,\frac{\sqrt{I}}{1 + b\,\sqrt{I}}\qquad
        f_m = 0

    The contribution expects the following parameters:

    =========== ================== =========================== =======
    Name        Symbol             Description                 Unit
    =========== ================== =========================== =======
    RHO_SOLVENT :math:`\rho_s`     Standard density of solvent [kg/m3]
    EPS_SOLVENT :math:`\epsilon_r` Dielectric constant solvent [-]
    B           :math:`b`          Ion size parameter          [-]
    =========== ================== =========================== =======

    Generally, :math:`b = 1.2` is used universally. For water, normally
    :math:`\rho_s = 980\ \mathrm{kg/m3}` and :math:`\epsilon_r = 80`.
    """
    def define_chi(self, res):
        temp, ionic_strength = res["T"], res["I"]
        rho_sol = self.par_scalar("RHO_SOLVENT", "kg/m**3")
        eps_r = self.par_scalar("EPS_SOLVENT", "dimless")
        b = self.par_scalar("B", "dimless")

        a_gamma = (
            sqrt(2 * PI * N_A * _M0 * rho_sol)
            * (E_0 ** 2 / (4 * PI * EPS_0 * eps_r * K_B * temp)) ** 1.5)
        b_sqi = b * (sqi := sqrt(ionic_strength))
        f = -4 / 3  * a_gamma * ionic_strength * log(1 + b_sqi) / b
        res["_pdh_f"] = f
        return {
            "chi": f,
            "chi_t": -1.5 * f / temp,
            "chi_i": f / ionic_strength - a_gamma * sqi / (1 + b_sqi) / 1.5,
            "chi_m": Quantity(0.0)}


class PitzerBinaryInteraction(ExcessBasePitzer):
    r"""The binary interaction in the Pitzer model is dependent on temperature
    and ionic strength as:

    .. math::

        \lambda_{ij}(T, I) = \beta^{(0)}_{ij}(T) +
          \frac{\beta^{(1)}_{ij}(T)}{2I}\left [
            1-(1+2\sqrt{I})\,\exp(-2\sqrt{I}) \right ]

    with

    .. math::

         \beta^{(k)}_{ij}(T) = \beta^{(k)}_{ij,1} +
           \beta^{(k)}_{ij,2}(T-\Theta) +
           \beta^{(k)}_{ij,3}\left (\frac1T-\frac1{\Theta}\right ) +
         \beta^{(k)}_{ij,4}\ln\frac{T}{\Theta} +
          \beta^{(k)}_{ij,5}\left (T^2-\Theta^2\right )

    In other implementations, the interaction coefficients are assumed
    symmetric and hence accounted for twice. This is not the case in this
    fully sparse implementation. For each species pair :math:`i,j`, a
    contribution is defined as

    .. math::

        \Delta\chi = m_i\,m_j\,\lambda_{ij}(T, I)

    The derivatives are coded manually with

    .. math::

        \beta^{(k)}_{ij,T} = \beta^{(k)}_{ij,2} -
           \beta^{(k)}_{ij,3} \frac1{T^2} + \beta^{(k)}_{ij,4}\,\frac{1}{T} +
          2\,\beta^{(k)}_{ij,5}\,T

    as

    .. math::
       :nowrap:

       \begin{align*}
        \Delta\chi_T &= m_i\,m_j\,\lambda_{ij,T}(T, I)\quad\text{with}\quad
            \lambda_{ij,T}(T, I) = \beta^{(0)}_{ij,T} +
            \frac{\beta^{(1)}_{ij,T}}{2I}\left [
            1-(1+2\sqrt{I})\,\exp(-2\sqrt{I}) \right ]\\
        \Delta\chi_I &= m_i\,m_j\,\lambda_{ij,I}(T, I)\quad\text{with}\quad
          \lambda_{ij,I}(T, I) = \beta^{(1)}_{ij,T}\,
          \frac{1 + (2\,I + 2\,\sqrt{I} - 1)\,\exp(-2\sqrt{I})}{2\,I^2}\\
        \Delta\boldsymbol{\chi}_m &= \lambda_{ij,T}(T, I)\,(
          m_i\,\mathbf{e}_j + m_j\,\mathbf{e}_i)
       \end{align*}

    """
    def define_chi(self, res):
        temp, molality, ionic_strength = res["T"], res["molality"], res["I"]
        t_ref = self.par_scalar("T_ref", "K")
        tsi = 2 * sqrt(ionic_strength)
        i_factor = (1 - (1 + tsi) * exp(-tsi))/ (2 * ionic_strength)
        i_factor_i = ((1 + (2 * ionic_strength + tsi - 1) * exp(-tsi)) /
                      (2 * ionic_strength ** 2))

        # pre-factors for binary interactions
        factors = [1, temp - t_ref, 1 / temp - 1 / t_ref, log(temp / t_ref),
                   temp ** 2 - t_ref ** 2]
        factors = [[factors], [f_i * i_factor for f_i in factors]]

        factors_t = [0, 1, -1 / temp ** 2, 1 / temp, 2 * temp]
        factors_t = [[factors_t], [f_i * i_factor for f_i in factors]]
        factors_i = [[0.0] * 5, [f_i * i_factor_i for f_i in factors[0]]]

        units = ["dimless", "1/K", "K", "dimless", "K**-2"]
        cache = {}

        def pair(idx_i: int, idx_j: int) -> Quantity:
            if (idx_i, idx_j) not in cache:
                cache[(idx_i, idx_j)] = molality[idx_i] * molality[idx_j]
            return cache[(idx_i, idx_j)]

        chi, chi_t, chi_i = Quantity(0), Quantity(0, "1/K"), Quantity(0)
        chi_m = Quantity(DM.zeros(len(self.species)))

        for k in (0, 1):
            for m in range(5):
                p_name = f"beta_{k}{m+1}"
                try:
                    pairs = self.options[p_name]
                except KeyError:
                    continue
                coefficients = self.par_sparse_matrix(p_name, pairs, units[m])
                f, f_t, f_i = factors[k][m], factors_t[k][m], factors_i[k][m]
                for i, j, c in coefficients.pair_items():
                    ii, ij = self.species.index(i), self.species.index(j)
                    term = pair(ii, ij) * c
                    chi += term * f
                    chi_t += term * f_t
                    chi_i += term * f_i
                    chi_m[ii] += term * molality[ij]
                    chi_m[ij] += term * molality[ii]

        return {"chi": chi, "chi_t": chi_t, "chi_i": chi_i, "chi_m": chi_m}






