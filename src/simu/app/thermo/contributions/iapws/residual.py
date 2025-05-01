from abc import abstractmethod
from collections.abc import Sequence

# internal modules
from simu import ThermoContribution, R_GAS, qvertcat, Quantity, exp
from simu.core.utilities.types import Map
from simu.core.utilities.quantity import jacobian

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
        num_range = range(1, number_of_terms + 1)

        props = "T V _tau _rho n".split()
        temp, vol, tau, rho, n = [res[i] for i in props]

        # select active species from tau, rho and n
        idx = [self.species.index(s) for s in species]
        tau = qvertcat(*[tau[i] for i in idx])
        rho = qvertcat(*[rho[i] for i in idx])
        n_sub = qvertcat(*[n[i] for i in idx])

        p_names = self.parameter_names()
        params = {n: [self.par_vector(f"{n}_{i:02d}", species, "")
                      for i in num_range]
                  for n in p_names}
        phi = self.define_phi(tau, rho, params)
        a_res = R_GAS * temp * (n_sub.T @ phi)

        res["mu"] += jacobian(a_res, n)
        res["S"] -= jacobian(a_res, temp)
        res["p"] -= jacobian(a_res, vol)

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
        on instantiation."""
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


class Residual1IAPWS(ResidualBaseIAPWS):
    r"""The first contribution of the IAPWS residual Helmholtz function, as
    represented by this contribution, is formulated as in :cite:p:`Wagner_2002`:

    .. math::

        \phi_i^{\mathrm{res}, 1} =
        \sum_{k=1}^{7} n_{k, i}^\mathrm{res}\,
          \varrho_i^{d_{k, i}}\,
          \tau_i^{t_{k, i}}\,

    Here, :math:`d_{k, i}` and :math:`t_{k, i}` are mostly integer coefficients
    or simple fractions, while :math:`n_{k, i}^\mathrm{res}` are the real tuning
    parameters.

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
            n_i * rho ** d_i * tau ** t_i
            for d_i, t_i, n_i in zip(*param)
        )

class Residual2IAPWS(ResidualBaseIAPWS):
    r"""This contribution defines the second group of terms in the residual
    Helmholtz energy, defined as in :cite:p:`Wagner_2002`:

    .. math::

        \phi_i^{\mathrm{res}, 2} =
           \sum_{k=1}^{44} n_{k, i}^\mathrm{res}\,
              \varrho_i^{d_{k, i}}\,
              \tau_i^{t_{k, i}}\,
              \exp \left ( -\varrho_i^{c_{k, i}}\right )
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
            n_i * rho ** d_i * tau ** t_i * exp(-rho ** c_i)
            for c_i, d_i, t_i, n_i in zip(*param)
        )

class Residual3IAPWS(ResidualBaseIAPWS):
    r"""This contribution defines the third group of terms in the residual
    Helmholtz energy, defined as in :cite:p:`Wagner_2002`:

    .. math::

        \phi_i^{\mathrm{res}, 2} =
           \sum_{k=1}^{3} n_{k, i}^\mathrm{res}\,
              \varrho_i^{d_{k, i}}\,
              \tau_i^{t_{k, i}}\,
              \exp \left [
                -\alpha_i\,(\varrho_i - \epsilon_i)^2
                -\beta_i\,(\tau_i - \gamma_i)^2
              \right ]
    """
    @staticmethod
    def parameter_names():
        return "d t n a b g e".split()

    @staticmethod
    def default_number_of_terms():
        return 3

    @staticmethod
    def define_phi(tau, rho, parameters):
        param = [parameters[i] for i in Residual2IAPWS.parameter_names()]
        return sum(
            n_i * rho ** d_i * tau ** t_i * exp(
                -a_i * (rho  - e_i) ** 2 - b_i * (tau - g_i) ** 2
            ) for d_i, t_i, n_i, a_i, b_i, g_i, e_i in zip(*param)
        )

class Residual4IAPWS(ResidualBaseIAPWS):
    r"""This contribution defines the forth group of terms in the residual
    Helmholtz energy, defined as in :cite:p:`Wagner_2002`:

    .. math::

        \phi_i^{\mathrm{res}, 2} =
           \sum_{k=1}^{3} n_{k, i}^\mathrm{res}\,
            \Delta_i^{b_i}\,\varrho\,\psi_i

    with

    .. math::

        \Delta_i &= \theta_i^2 + B_i\,\hat\varrho_i^{a_i}\\
        \theta &= 1 - \tau_i + A_i\,\hat\varrho_i^{1/(2\,\beta_i)}\\
        \psi_i &= \exp \left [
          -C_i\,\hat\varrho_i - D_i\,(\tau_i - 1)^2
          \right ]\\
        \hat \varrho_i &= (\varrho_i - 1)^2
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
            [parameters[i] for i in Residual2IAPWS.parameter_names()]
        rho_hat = (1 - rho) ** 2
        psi = [exp(-C_i * rho_hat - D_i * (tau - 1) ** 2)
               for C_i, D_i in zip(C, D)]
        theta = [1 - tau + A_i * rho_hat ** (1 / (2 * beta_i))
                 for A_i, beta_i in zip(A, beta)]
        delta = [theta_i ** 2 + B_i * rho_hat ** a_i
                 for theta_i, B_i, a_i in zip(theta, B, a)]
        return sum(
            n_i * delta_i ** b_i * rho * psi_i
            for n_i, delta_i, b_i, psi_i  in zip(n, delta, b, psi)
        )
