from abc import abstractmethod

from casadi import SX, DM
from simu import (ThermoContribution, registered_contribution, Quantity,
                  N_A, E_0, EPS_0, K_B, R_GAS, PI, sqrt, log, exp, qvertcat)
from simu.core.utilities.types import Map, MutMap

_M0 = Quantity(1.0, "mol/kg")
_C0 = Quantity(1.0, "e/mol")


@registered_contribution
class ElectrolyteBasics(ThermoContribution):
    r"""This contribution prepares some basic properties relevant for
    electrolyte systems.

    **Molality** (:math:`b_i`) is defined as the molar quantities per mass of
    solvent. Generalized, we consider all non-ionic species as solvent. With
    the charge vector :math:`c_i`, the Kronecker symbol for marking solvent
    components is

    .. math::

        \delta_{si} = \begin{cases}
            1\quad\text{for}\ c_i = 0\\
            0\quad\text{else} \end{cases}

    As such, :math:`m_s = \sum_i n_i\,M_i\,b_0\,\delta_{si}` and molalities are
    defined as :math:`b_i = n_i / m_s`. Here, :math:`b_0 = 1` mol/kg
    is a common factor used to yield dimensionless molalities.

    Based on molality, **ionic strength** :math:`I` is defined as

    .. math:: I = \frac12\,\sum_i b_i\,\left (\frac{c_i}{c_0}\right )^2

    Again, to support the dimensionless mind of electro-chemists for empirical
    freedom, the charge is normalized by :math:`c_0 = 1` e/mol.

    To derive molality-based Gibbs excess contributions with respect to molar
    quantities, we pre-calculate the following derivatives:

    .. math::

       \frac{\mathrm{d} m_s}{\mathrm{d} n_k} = M_k\,b_0\,\delta_{sk}\qquad
       \frac{\mathrm{d} I}{\mathrm{d} b_k} =
         \frac12 \left (\frac{c_k}{c_0} \right )^2
    """

    provides = ["b", "I", "_m_s", "_di_db", "_dms_dn", "charge"]

    def define(self, res):
        n, mw = res["n"], res["mw"]
        species_def = self.species_definitions
        charge_list = [s.charge for s in species_def.values()]
        kron_s = qvertcat(*[Quantity(1 if c == 0 else 0) for c in charge_list])
        res["charge"] = charge = qvertcat(*charge_list)  # [e/mol]

        # define molality, ionic strength, and relevant derivatives
        res["_di_db"] = di_db = (charge / _C0) ** 2 / 2
        res["_dms_dn"] = d_ms_dn = mw * kron_s * _M0 # [-]
        res["_m_s"] = m_s = n.T @ d_ms_dn  # [mol (mol/s)]
        res["b"] = b = n / m_s  # [-]
        res["I"] = b.T @ di_db


class ExcessBasePitzer(ThermoContribution):
    r"""The Pitzer model is formulated in terms of a reduced excess Gibbs
    energy contributions follows:

    .. math::

        \frac{\Delta \chi G^\mathrm{ex}}{R\,T} = m_s\,\chi(T, b_i)

    Both the long-range contribution (Pitzer-Debye-Hückel) and the short-range
    contribution expressed by binary and ternary parameters are expressed in
    terms of :math:`\chi(T, b_i)`.

    The chemical potential is then

    .. math::

        \frac{\Delta_\chi \mu_i}{R\,T} = \left [
            \chi - \left .
              \sum_k\frac{\partial \chi}{\partial b_k}\right |_T\,\,b_k
          \right ]\,\frac{\mathrm{d} m_s}{\mathrm{d} n_i} +
          \left . \frac{\partial \chi}{\partial b_i} \right |_T

    The entropy is

    .. math::

        \Delta_\chi S = R\,m_s\,\left [
            \chi + T\,\left .\frac{\partial \chi}{\partial T} \right |_{n}
            \right ]

    From the subclasses, the :meth:`define_chi` method is to return for
    convenience partial derivatives with respect to molalities and ionic
    strength separately. The required combined derivative is then

    .. math::

        \left . \frac{\partial \chi}{\partial b_i} \right |_T
         = \left . \frac{\partial \chi}{\partial b_i} \right |_{T, I} +
         \left . \frac{\partial \chi}{\partial I} \right |_{T, b}\,
        \frac{\mathrm{d} I}{\mathrm{d} b_k}

    """
    def define(self, res):
        names = ["T", "b", "I", "_m_s", "_di_db", "_dms_dn"]
        temp, b, ios, m_s, didb, dmsdn = [res[n] for n in names]

        chi_res = self.define_chi(res)
        chi, chi_t = chi_res["chi"], chi_res["chi_t"]
        chi_b = chi_res["chi_b"] +  chi_res["chi_i"] * didb

        res["S"] += R_GAS * m_s * (chi + temp * chi_t)
        res["mu"] += R_GAS * temp * (chi_b + (chi - chi_b.T @ b) * dmsdn)

    @abstractmethod
    def define_chi(self, res: MutMap[Quantity]) -> Map[Quantity]:
        """
        Provide dimensionless residual contribution :math:`\chi(T, I, b_i)`
        and the partial derivatives :math:`\chi_T` (``chi_t``),
        :math:`\chi_I` (``chi_i``) and :math:`\chi_{b}` (``chi_b``).
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

       \chi^\mathrm{PDH}(I) = -\frac43A_\gamma I^{3/2}\,
         \frac{\ln (1+b\sqrt{I})}{b\sqrt{I}} \quad {\rm with}\quad b=1.2

    The function :math:`\chi^\mathrm{PDH}(I)` is a compatible dimensionless
    contribution  as defined for the :class:`ExcessBasePitzer` base-class.

    The required derivatives are

    .. math::
        \chi^\mathrm{PDH}_T = -\frac32\,\frac{\chi^\mathrm{PDH}}{T}\qquad
        \chi^\mathrm{PDH}_I = \frac{\chi^\mathrm{PDH}}{I} -
          \frac23\,A_\gamma\,\frac{\sqrt{I}}{1 + b\,\sqrt{I}}\qquad
        \chi^\mathrm{PDH}_b = 0

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
        chi = -4 / 3  * a_gamma * ionic_strength * log(1 + b_sqi) / b
        res["chi_pdh_f"] = chi
        return {
            "chi": chi,
            "chi_t": -1.5 * chi / temp,
            "chi_i": chi / ionic_strength - a_gamma * sqi / (1 + b_sqi) / 1.5,
            "chi_b": Quantity(0.0)}


@registered_contribution
class PitzerBinaryInteraction(ExcessBasePitzer):
    r"""The binary interaction in the Pitzer model is dependent on temperature
    and ionic strength as:

    .. math::

        \lambda_{ij}(T, I) = \beta^{(0)}_{ij}(T) +
          \frac{1-(1+2\sqrt{I})\,\exp(-2\sqrt{I})}{2\,I}\,\beta^{(1)}_{ij}(T)

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

        \Delta_\lambda\chi = b_i\,b_j\,\lambda_{ij}(T, I)

    The derivatives are coded manually with

    .. math::

        \beta^{(k)}_{ij,T} = \beta^{(k)}_{ij,2} -
           \beta^{(k)}_{ij,3} \frac1{T^2} + \beta^{(k)}_{ij,4}\,\frac{1}{T} +
          2\,\beta^{(k)}_{ij,5}\,T

    as

    .. math::
       :nowrap:

       \begin{align*}
        \Delta_\lambda\chi_T &= b_i\,b_j\,\lambda_{ij,T}(T, I)\quad\text{with}\quad
            \lambda_{ij,T}(T, I) = \beta^{(0)}_{ij,T} +
            \frac{1-(1+2\sqrt{I})\,\exp(-2\sqrt{I})}{2\,I}\,\beta^{(1)}_{ij,T}\\
        \Delta_\lambda\chi_I &= b_i\,b_j\,\lambda_{ij,I}(T, I)\quad\text{with}\quad
          \lambda_{ij,I}(T, I) =
            \frac{1 + (2\,I + 2\,\sqrt{I} - 1)\,\exp(-2\sqrt{I})}{2\,I^2}\,
            \beta^{(1)}_{ij}\\
        \Delta_\lambda\boldsymbol{\chi}_b &= \lambda_{ij,T}(T, I)\,(
          b_i\,\mathbf{e}_j + b_j\,\mathbf{e}_i)
       \end{align*}

    """
    def define_chi(self, res):
        temp, b, ios = res["T"], res["b"], res["I"]
        t_ref = self.par_scalar("T_ref", "K")
        tsi = 2 * sqrt(ios)
        i_factor = (1 - (1 + tsi) * exp(-tsi))/ (2 * ios)
        i_factor_i = ((1 + (2 * ios + tsi - 1) * exp(-tsi)) / (2 * ios ** 2))

        # pre-factors for binary interactions
        # Operations involving factors of zero are required to provide the
        # correct unit of measurement, e.g. 0 / temp = 0 1/K.
        factors = [1, temp - t_ref, 1 / temp - 1 / t_ref, log(temp / t_ref),
                   temp ** 2 - t_ref ** 2]
        factors = [factors, [f_i * i_factor for f_i in factors]]

        factors_t = [0 / temp, 1, -1 / temp ** 2, 1 / temp, 2 * temp]
        factors_t = [factors_t, [f_i * i_factor for f_i in factors_t]]
        factors_i = [[0 * f_i for f_i in factors[0]],
                     [f_i * i_factor_i for f_i in factors[0]]]

        units = ["dimless", "1/K", "K", "dimless", "K**-2"]
        cache = {}

        def pair(idx_i: int, idx_j: int) -> Quantity:
            if (idx_i, idx_j) not in cache:
                cache[(idx_i, idx_j)] = b[idx_i] * b[idx_j]
            return cache[(idx_i, idx_j)]

        chi, chi_t, chi_i = Quantity(0), Quantity(0, "1/K"), Quantity(0)
        chi_b = Quantity(SX.zeros(len(self.species)))

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
                    d_chi = term * f
                    chi += d_chi
                    chi_t += term * f_t
                    chi_i += term * f_i
                    chi_b[ii] += d_chi * b[ij]
                    chi_b[ij] += d_chi * b[ii]

        res["_pitzer_bin_chi"] = chi
        return {"chi": chi, "chi_t": chi_t, "chi_i": chi_i, "chi_b": chi_b}

@registered_contribution
class PitzerTernaryInteraction(ExcessBasePitzer):
    r"""The ternary interaction in the Pitzer model is dependent only on
    temperature:

    .. math:: \Delta_\gamma \chi_{ijk} = \gamma_{ijk}(T)\,b_i\,b_j\,b_k

    The interaction coefficients are parameterized as

    .. math::

        \gamma_{ijk}(T) = \gamma_{ijk,1} + \gamma_{ijk,2}(T-\Theta) +
            \gamma_{ijk,3}\left (\frac1T-\frac1{\Theta}\right ) +
            \gamma_{ijk,4}\ln\frac{T}{\Theta} +
            \gamma_{ijk,5}\left (T^2-\Theta^2\right )

    The derivatives are coded manually with

    .. math::

        \gamma_{ijk,T} = \gamma_{ijk,2} - \gamma_{ijk,3} \frac1{T^2} +
            \gamma_{ijk,4}\,\frac{1}{T} + 2\,\gamma_{ijk,5}\,T

    The required derivatives are provided analytically:

    .. math::

         \Delta_{\gamma} \chi_{ijk, T} = \gamma_{ijk_T}\,b_i\,b_j\,b_k\qquad
         \Delta_{\gamma} \chi_{ijk, I} = 0\qquad
         \Delta_{\gamma} \chi_{ijk, m} = \gamma_{ijk}\,\left (
            b_i\,b_j\,\mathbf{e}_k + b_i\,b_k\,\mathbf{e}_j
            + b_j\,b_k\,\mathbf{i}_k
         \right )
    """
    def define_chi(self, res):
        temp, b = res["T"], res["b"]
        t_ref = self.par_scalar("T_ref", "K")

        factors = [1, temp - t_ref, 1 / temp - 1 / t_ref, log(temp / t_ref),
                   temp ** 2 - t_ref ** 2]
        factors_t = [0 / temp, 1, -1 / temp ** 2, 1 / temp, 2 * temp]

        units = ["dimless", "1/K", "K", "dimless", "K**-2"]
        cache = {}

        def pair(idx_i: int, idx_j: int, idx_k: int) -> Quantity:
            if (idx_i, idx_j, idx_k) not in cache:
                cache[(idx_i, idx_j, idx_k)] = \
                    b[idx_i] * b[idx_j] * b[idx_k]
            return cache[(idx_i, idx_j, idx_k)]

        chi, chi_t, chi_i = Quantity(0), Quantity(0, "1/K"), Quantity(0)
        chi_b = Quantity(SX.zeros(len(self.species)))

        for m in range(5):
            p_name = f"gamma_{m+1}"
            try:
                pairs = self.options[p_name]
            except KeyError:
                continue
            coefficients = self.par_sparse_3d(p_name, pairs, units[m])
            f, f_t = factors[m], factors_t[m]
            for i, j, k, c in coefficients.pair_items():
                ii, ij, ik = map(self.species.index, (i, j, k))
                term = pair(ii, ij, ik) * c
                d_chi = term * f
                chi += d_chi
                chi_t += term * f_t
                chi_b[ii] += d_chi * b[ij] * b[ik]
                chi_b[ij] += d_chi * b[ii] * b[ik]
                chi_b[ik] += d_chi * b[ii] * b[ij]

        res["_pitzer_ternary_chi"] = chi
        return {"chi": chi, "chi_t": chi_t, "chi_i": chi_i, "chi_b": chi_b}
