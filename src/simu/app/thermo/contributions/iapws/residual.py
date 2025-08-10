from abc import abstractmethod
from collections.abc import Sequence
from numpy import isnan

# internal modules
from simu.core.thermo.contribution import ThermoContribution, registered_contribution
from simu.core.utilities.constants import R_GAS
from simu.core.utilities.quantity import (
    Quantity, qsum, qpow, qvertcat, QFunction, SymbolQuantity,
    base_magnitude, jacobian)
from simu.core.utilities.qstructures import exp
from simu.core.utilities.types import Map

# TODO:
#  - document parameters required for each contribution


class ResidualBaseIAPWS(ThermoContribution):
    r"""All IAPWS residual contributions define a molar residual contribution
    :math:`\phi_i^\mathrm{res}(\tau_i, \varrho_i)` being the sum of :math:`m`
    terms, and applying for a subset of species :math:`i`.

    The contribution can be parameterized by two options:

      - ``species``: If provided, the contribution is only applied to a sub-set
        of species. This allows to combine the rigorous residual contribution
        with simplified models for secondary species.
      - ``number_of_terms``: The default number of terms for this contribution
        is :math:`m = 7`. This parameter allows to change the number to
        accommodate different model structures.

    The Helmholtz function is then implemented as

    .. math::

        A^\mathrm{res} = R\,T\,\sum_i n_i\,\phi_i^\mathrm{res}

    `Casadi`_ is used to obtain the derivatives as

    .. math::

        S \leftarrow S - \left .
            \frac{\partial A^\mathrm{res}}{\partial T}
          \right |_{V, n_i}\\
        p \leftarrow p - \left .
            \frac{\partial A^\mathrm{res}}{\partial V}
          \right |_{T, n_i}\\
        \mu_i \leftarrow \mu_i + \left .
            \frac{\partial A^\mathrm{res}}{\partial n_i}
          \right |_{T, V}

    """
    def define(self, res):
        species = self.options.get("species", self.species)
        # Only take species that are in self.species
        species = [s for s in species if s in self.species]
        number_of_terms = self.options.get("number_of_terms",
                                           self.default_number_of_terms())
        props = "T V _tau _rho n".split()
        temp, vol, tau, rho, n = [res[i] for i in props]

        # select active species from tau, rho and n
        idx = [self.species.index(s) for s in species]
        tau = qvertcat(*[tau[i] for i in idx])
        rho = qvertcat(*[rho[i] for i in idx])
        n_sub = qvertcat(*[n[i] for i in idx])

        p_names = self.parameter_names()
        num_terms_rng = range(1, number_of_terms + 1)
        params = {n: [self.par_vector(f"{n}_{i:02d}", species, "")
                      for i in num_terms_rng]
                  for n in p_names}
        phi = self._get_phi(tau, rho, params, species, number_of_terms)

        # use chain rule to construct derivatives
        # Doesn't make a big difference to brute force derivatives, however :(
        s_res = -R_GAS * (n_sub.T @ (phi["phi"] - tau.T @ phi["phi_tau"]))
        p_res = R_GAS * temp / vol * (n_sub.T @ phi["phi_rho"] @ rho)
        mu_res = R_GAS * temp * (phi["phi"] + rho @ phi["phi_rho"])

        res["mu"] += mu_res
        res["S"] += s_res
        res["p"] += p_res

    def _get_phi(self, tau: Quantity, rho: Quantity,
                 params: Map[Sequence[Quantity]], species: Sequence[str],
                 number_of_terms: int) -> Map[Quantity]:
        """Define phi as a casadi function and evaluate the derivatives
        with respect to tau and rho."""
        tau_independent = SymbolQuantity("tau", "", species)
        rho_independent = SymbolQuantity("rho", "", species)
        phi = self.define_phi(tau_independent, rho_independent, params)
        phi_tau = jacobian(phi, tau_independent)
        phi_rho = jacobian(phi, rho_independent)
        p_struct = {n: {f"{i:02d}": p_n[i] for i in range(number_of_terms)}
                    for n, p_n in params.items()}
        args = {"tau": tau_independent, "rho": rho_independent,
                "params": p_struct}
        result = {"phi": phi, "phi_rho": phi_rho, "phi_tau": phi_tau}
        f = QFunction(args, result)
        result = f({"tau": tau, "rho": rho, "params": p_struct})
        return result

    @staticmethod
    @abstractmethod
    def parameter_names() -> Sequence[str]:
        """Return a sequence of parameter names to be required by the
        contribution."""
        ...

    @staticmethod
    @abstractmethod
    def default_number_of_terms() -> int:
        """Return the default number of terms for this contribution.
        This can be over-defined by the option data given to the contribution
        object on instantiation."""
        ...

    @staticmethod
    @abstractmethod
    def define_phi(tau: Quantity, rho: Quantity,
                   parameters: Map[Sequence[Quantity]]) -> Quantity:
        r"""Based on the map of parameters, whereas the keys are the names
        defined in :meth:`parameter_names`, define the residual
        :math:`\phi_i^\mathrm{res}(\tau_i, \varrho_i)` vector function, such
        that the Helmholtz function contribution can be constructed as

        .. math::

            A^{\mathrm{res}, j} = R\,T\,\sum_i n_i\,\phi_i^\mathrm{res}
        """
        ...


@registered_contribution
class Residual1IAPWS(ResidualBaseIAPWS):
    r"""The first contribution of the IAPWS residual Helmholtz function, as
    represented by this contribution, is formulated as in :cite:p:`Wagner_2002`:

    .. math::

        \phi_i^{\mathrm{res}, 1} =
        \sum_{k=1}^{m} n_{k, i}\,
          \varrho_i^{d_{k, i}}\,
          \tau_i^{t_{k, i}}\,

    The parameters are defined as follows, the indices formatted with two
    digits. The default number of terms :math:`m` for this contribution is 44.

    ==================== ==================================
    Parameter            Symbol
    ==================== ==================================
    ``n_01`` to ``n_mm`` :math:`n_{1,i}^\mathrm{res,1}` to
                         :math:`n_{m,i}^\mathrm{res,1}`
    ``d_01`` to ``d_mm`` :math:`d_{1,i}` to :math:`d_{m,i}`
    ``t_01`` to ``t_mm`` :math:`t_{1,i}` to :math:`t_{m,i}`
    ==================== ==================================
    """
    @staticmethod
    def parameter_names():
        return "d t n".split()

    @staticmethod
    def default_number_of_terms():
        return 7

    @staticmethod
    def define_phi(tau, rho, parameters):
        param = [parameters[i] for i in Residual1IAPWS.parameter_names()]
        return sum(
            n_i * qpow(rho, d_i) * qpow(tau, t_i)
            for d_i, t_i, n_i in zip(*param)
        )

@registered_contribution
class Residual2IAPWS(ResidualBaseIAPWS):
    r"""This contribution defines the second group of terms in the residual
    Helmholtz energy, defined as in :cite:p:`Wagner_2002`:

    .. math::

        \phi_i^{\mathrm{res}, 2} =
           \sum_{k=1}^{m} n_{k, i}\,
              \varrho_i^{d_{k, i}}\,
              \tau_i^{t_{k, i}}\,
              \exp \left ( -\varrho_i^{c_{k, i}}\right )

    The parameters are defined as follows, the indices formatted with two
    digits. The default number of terms :math:`m` for this contribution is 44.

    ==================== ==================================
    Parameter            Symbol
    ==================== ==================================
    ``n_01`` to ``n_mm`` :math:`n_{1,i}` to :math:`n_{m,i}`
    ``d_01`` to ``d_mm`` :math:`d_{1,i}` to :math:`d_{m,i}`
    ``t_01`` to ``t_mm`` :math:`t_{1,i}` to :math:`t_{m,i}`
    ``c_01`` to ``c_mm`` :math:`c_{1,i}` to :math:`c_{m,i}`
    ==================== ==================================
    """
    @staticmethod
    def parameter_names():
        return "c d t n".split()

    @staticmethod
    def default_number_of_terms():
        return 44

    @staticmethod
    def define_phi(tau, rho, parameters):
        param = [parameters[i] for i in Residual2IAPWS.parameter_names()]
        return sum(
            n_i * qpow(rho, d_i) * qpow(tau, t_i) * exp(-qpow(rho, c_i))
            for c_i, d_i, t_i, n_i in zip(*param)
        )


@registered_contribution
class Residual3IAPWS(ResidualBaseIAPWS):
    r"""This contribution defines the third group of terms in the residual
    Helmholtz energy, defined as in :cite:p:`Wagner_2002`:

    .. math::

        \phi_i^{\mathrm{res}, 3} =
           \sum_{k=1}^{m} n_{k, i}^\mathrm{res}\,
              \varrho_i^{d_{k, i}}\,
              \tau_i^{t_{k, i}}\,
              \exp \left [
                -\alpha_{k,i}\,(\varrho_i - \epsilon_{k,i})^2
                -\beta_{k,i}\,(\tau_i - \gamma_{k,i})^2
              \right ]

    The parameters are defined as follows, the indices formatted with two
    digits. The default number of terms :math:`m` for this contribution is 3.

    ==================== ================================================
    Parameter            Symbol
    ==================== ================================================
    ``n_01`` to ``n_mm`` :math:`n_{1,i}` to :math:`n_{m,i}`
    ``d_01`` to ``d_mm`` :math:`d_{1,i}` to :math:`d_{m,i}`
    ``t_01`` to ``t_mm`` :math:`t_{1,i}` to :math:`t_{m,i}`
    ``a_01`` to ``a_mm`` :math:`\alpha_{1,i}` to :math:`\alpha_{m,i}`
    ``b_01`` to ``b_mm`` :math:`\beta_{1,i}` to :math:`\beta_{m,i}`
    ``g_01`` to ``g_mm`` :math:`\gamma_{1,i}` to :math:`\gamma_{m,i}`
    ``e_01`` to ``e_mm`` :math:`\epsilon_{1,i}` to :math:`\epsilon_{m,i}`
    ==================== ================================================
    """
    @staticmethod
    def parameter_names():
        return "d t n a b g e".split()

    @staticmethod
    def default_number_of_terms():
        return 3

    @staticmethod
    def define_phi(tau, rho, parameters):
        param = [parameters[i] for i in Residual3IAPWS.parameter_names()]
        return sum(
            n_i * qpow(rho, d_i) * qpow(tau, t_i) * exp(
                -a_i * (rho  - e_i) ** 2 - b_i * (tau - g_i) ** 2
            ) for d_i, t_i, n_i, a_i, b_i, g_i, e_i in zip(*param)
        )


@registered_contribution
class Residual4IAPWS(ResidualBaseIAPWS):
    r"""This contribution defines the forth group of terms in the residual
    Helmholtz energy, defined as in :cite:p:`Wagner_2002`:

    .. math::

        \phi_i^{\mathrm{res}, 4} =
           \sum_{k=1}^{m} n_{k, i}^\mathrm{res}\,
            \Delta_{k,i}^{b_{k,i}}\,\varrho_i\,\psi_{k,i}

    with

    .. math::

        \Delta_{k,i} &= \theta_{k,i}^2 + B_i\,\hat\varrho_i^{a_{k,i}}\\
        \theta_{k,i} &= 1-\tau_i + A_{k,i}\,\hat\varrho_i^{1/(2\,\beta_{k,i})}\\
        \psi_i &= \exp \left [
          -C_{k,i}\,\hat\varrho_i - D_{k,i}\,(\tau_i - 1)^2
          \right ]\\
        \hat \varrho_i &= (\varrho_i - 1)^2

    The parameters are defined as follows, the indices formatted with two
    digits. The default number of terms :math:`m` for this contribution is 2.

    ========================== ================================================
    Parameter                  Symbol
    ========================== ================================================
    ``n_01`` to ``n_mm``       :math:`n_{1,i}` to :math:`n_{m,i}`
    ``a_01`` to ``a_mm``       :math:`a_{1,i}` to :math:`a_{m,i}`
    ``b_01`` to ``b_mm``       :math:`b_{1,i}` to :math:`b_{m,i}`
    ``B_01`` to ``B_mm``       :math:`B_{1,i}` to :math:`B_{m,i}`
    ``C_01`` to ``C_mm``       :math:`C_{1,i}` to :math:`C_{m,i}`
    ``D_01`` to ``D_mm``       :math:`D_{1,i}` to :math:`D_{m,i}`
    ``A_01`` to ``A_mm``       :math:`A_{1,i}` to :math:`A_{m,i}`
    ``beta_01`` to ``beta_mm`` :math:`\beta_{1,i}` to :math:`\beta_{m,i}`
    ========================== ================================================

    """
    @staticmethod
    def parameter_names():
        return "a b B n C D A beta".split()

    @staticmethod
    def default_number_of_terms():
        return 2

    @staticmethod
    def define_phi(tau, rho, parameters):
        a, b, B, n, C, D, A, beta = \
            [parameters[i] for i in Residual4IAPWS.parameter_names()]
        rho_hat = (rho - 1) ** 2
        psi = [exp(-C_i * rho_hat - D_i * (tau - 1) ** 2)
               for C_i, D_i in zip(C, D)]
        theta = [1 - tau + A_i * qpow(rho_hat, 1 / (2 * beta_i))
                 for A_i, beta_i in zip(A, beta)]
        delta = [theta_i ** 2 + B_i * qpow(rho_hat, a_i)
                 for theta_i, B_i, a_i in zip(theta, B, a)]
        return sum(
            n_i * qpow(delta_i, b_i) * rho * psi_i
            for n_i, delta_i, b_i, psi_i  in zip(n, delta, b, psi)
        )


@registered_contribution
class GasIAPWSIdealMix(ThermoContribution):
    r"""This contribution does not add any terms to the state function itself,
    but prepares and provides the initialization if the IAPWS model gas phase
    based on empirical expressions for the pure species, paired with ideal
    mixing rules for multi-component mixtures.

    To find the gas volume, we calculate the saturation volume according to the
    empirical formula in :cite:p:`Wagner_2002`:

    .. math::

        \ln \varrho_i^* = \sum_{j=1}^m c_{j,i}\,(1-\tau_i)^{f_{j,i}}

    The total saturation volume according to Raoult's law is then

    .. math:: V^* = \sum_i \frac{n_i\, M_i}{\varrho_i^*\,\rho_{c,i}}

    Instead of initialising with the saturation volume itself, we compensate
    proportionally with any pressure departure. For that, the saturation
    pressure is evaluated according to :cite:p:`Wagner_2002` as

    .. math::

        \ln \frac{p^*_i}{p_{c,i}} = \tau_i\,\sum_{j=1}^{m}
          a_{j,i}\,(1-\tau_i)^{e_{j,i}}

    The total saturation pressure is the mole-fraction weighted sum of partial
    pressures:

    .. math:: p^* = \frac{1}{N}\sum_i n_i\,p^*_i

    As such, the volume estimate is calculated to

    .. math:: V^\mathrm{est} = \frac{p^*}{p}\,V^*

    This contribution requires the following parameters

    ========================== ================================================
    Parameter                  Symbol
    ========================== ================================================
    ``a_01`` to ``a_mm``       :math:`a_{1,i}` to :math:`a_{m,i}`
    ``e_01`` to ``e_mm``       :math:`e_{1,i}` to :math:`e_{m,i}`
    ``c_01`` to ``c_mm``       :math:`c_{1,i}` to :math:`c_{m,i}`
    ``f_01`` to ``f_mm``       :math:`f_{1,i}` to :math:`f_{m,i}`
    ========================== ================================================

    By default, the number of terms is :math:`m = 6`, but this can be altered
    by the ``number_of_terms`` option.
    """
    def define(self, res):
        props = "T n _tau _p_c _rho_c mw".split()
        temp, n, tau, p_c, rho_c, mw = [res[i] for i in props]
        number_of_terms = self.options.get("number_of_terms", 6)
        ntr = range(1, number_of_terms + 1)
        a, e, c, f = [[self.par_vector(f"{n}_{i:02d}", self.species, "")
                       for i in ntr] for n in ("a", "e", "c", "f")]
        theta = 1 - 1 / tau
        x = n / qsum(n)  # Raoult's law
        expr_p = tau * sum(a_i * qpow(theta, e_i) for a_i, e_i in zip(a, e))
        res["_p_sat"] = p_c * (x.T @ exp(expr_p))

        expr_r = sum(c_i * qpow(theta, f_i) for c_i, f_i in zip(c, f))
        rho = rho_c * exp(expr_r)
        res["_v_sat"] = qsum(n * mw / rho)  # no mixing volume

    def initial_state(self, state, properties):
        p_sat, v_sat = properties["_p_sat"], properties["_v_sat"]
        # check if super-critical, as then there is no v_sat or p_sat
        #  then, just initialize with ideal gas (there should only be one root)
        if isnan(p_sat) or isnan(v_sat):
            volume = (qsum(state.mol_vector) * R_GAS * state.temperature /
                      state.pressure)
        else:
            volume = p_sat / state.pressure * v_sat
        return ([base_magnitude(state.temperature),
                 base_magnitude(volume)] +
                list(base_magnitude(state.mol_vector)))


@registered_contribution
class LiquidIAPWSIdealMix(ThermoContribution):
    r"""This contribution does not add any terms to the state function itself,
    but prepares and provides the initialization if the IAPWS model liquid phase
    based on empirical expressions for the pure species, paired with ideal
    mixing rules for multi-component mixtures.

    To find the liquid volume, we calculate the saturation volume according to
    the empirical formula in :cite:p:`Wagner_2002`:

    .. math::

        \varrho_i^* = 1 + \sum_{j=1}^m b_{j,i}\,(1-\tau_i)^{f_{j,i}}

    The total saturation volume is then

    .. math::

       V^\mathrm{est} = V^* = \sum_i \frac{n_i\, M_i}{\varrho_i^*\,\rho_{c,i}}

    Unlike for initialising the gas phase (:class:`GasIAPWSIdealMix`), we do not
    pressure-compensate the saturation volume.

    ========================== ================================================
    Parameter                  Symbol
    ========================== ================================================
    ``b_01`` to ``b_mm``       :math:`b_{1,i}` to :math:`b_{m,i}`
    ``f_01`` to ``f_mm``       :math:`f_{1,i}` to :math:`f_{m,i}`
    ========================== ================================================

    By default, the number of terms is :math:`m = 6`, but this can be altered
    by the ``number_of_terms`` option.
    """
    def define(self, res):
        props = "T n _tau _p_c _rho_c mw".split()
        temp, n, tau, p_c, rho_c, mw = [res[i] for i in props]
        number_of_terms = self.options.get("number_of_terms", 6)
        ntr = range(1, number_of_terms + 1)
        b, f = [[self.par_vector(f"{n}_{i:02d}", self.species, "")
                for i in ntr] for n in ("b", "f")]
        theta = 1 - 1 / tau
        expr_r = 1 + sum(b_i * qpow(theta, f_i) for b_i, f_i in zip(b, f))
        rho = rho_c * expr_r
        res["_v_sat"] = qsum(n * mw / rho)

    def initial_state(self, state, properties):
        v_sat = properties["_v_sat"]
        # check if super-critical, as then there is no v_sat
        #  then, just initialize with ideal gas (there should only be one root)
        if isnan(v_sat):
            volume = (qsum(state.mol_vector) * R_GAS * state.temperature /
                      state.pressure)
        else:
            volume = v_sat
        return ([base_magnitude(state.temperature),
                 base_magnitude(volume)] +
                list(base_magnitude(state.mol_vector)))
