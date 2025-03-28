# -*- coding: utf-8 -*-

# stdlib modules
from copy import copy
from typing import Sequence

# internal modules
from simu import ThermoContribution, R_GAS, log, qsum, base_magnitude, qvertcat
from simu.core.utilities.errors import DimensionalityError


class H0S0ReferenceState(ThermoContribution):
    r"""This contribution defines the reference state based on enthalpy of
    formation and standard entropy. It requires the following parameters:

    =========== ================================== ======================
    Parameter   Description                        Symbol
    =========== ================================== ======================
    ``dh_form`` Molar enthalpy of formation vector :math:`\Delta_f h`
    ``s_0``     Vector of standard entropies       :math:`s_0`
    ``T_ref``   Reference temperature              :math:`T_\mathrm{ref}`
    ``p_ref``   Reference pressure                 :math:`T_\mathrm{ref}`
    =========== ================================== ======================

    Based on this, the reference state is defined in therms of entropy
    and chemical potentials:

    ========= ========================================= ======================
    Property  Description                               Symbol
    ========= ========================================= ======================
    ``S``     Reference state entropy of the system     :math:`S`
    ``mu``    Reference state chemical potential vector :math:`\mu`
    ``T_ref`` Reference temperature (same as parameter) :math:`T_\mathrm{ref}`
    ``p_ref`` Reference pressure (same as parameter)    :math:`p_\mathrm{ref}`
    ========= ========================================= ======================

    The contribution does not support any options.

    .. math::

        S &= \sum_i s_{0,i}\, n_i\\
        \mu_i &= \Delta_f h_i - T\,s_{0,i}
    """

    provides = ["T_ref", "p_ref", "S", "mu"]

    def define(self, res):
        species = self.species

        s_0 = self.par_vector("s_0", species, "J/(mol*K)")
        dh_form = self.par_vector("dh_form", species, "J/mol")
        res["S"] = s_0.T @ res["n"]
        res["mu"] = dh_form - res["T"] * s_0
        res["T_ref"] = self.par_scalar("T_ref", "K")
        res["p_ref"] = self.par_scalar("p_ref", "Pa")

        self.declare_vector_keys("mu")


class LinearHeatCapacity(ThermoContribution):
    r"""This contribution implements a simple heat capacity, being linear
    in temperature. It normally builds on a reference state and yields the
    standard state model. It requires the following parameters:

    ========= ====================================== ===============
    Parameter Description                            Symbol
    ========= ====================================== ===============
    ``cp_a``  Molar heat capacity coefficient vector :math:`c_{p,a}`
    ``cp_b``  Molar heat capacity coefficient vector :math:`c_{p,b}`
    ========= ====================================== ===============

    Heat capacity is here defined as

    .. math::

        \Delta T := T-T_\mathrm{ref}\\
        c_{p,i} = c_{p,a,i} + \Delta T\,c_{p,b,i}

    The contribution does not define new properties, but builds on entropy
    ``S`` (:math:`S`) and chemical potential ``mu`` (:math:`\mu`). With the
    following definitions:

    .. math::

        \Delta h_i &= \left (c_{p,a,i} + \frac12\Delta T\,c_{p,b,i}\right )\,
                      \Delta T\\
        \Delta s_i &= (c_{p,a,i} - c_{p,b,i}\,T_\mathrm{ref})\,
                      \ln \frac{T}{T_\mathrm{ref}} +
                      c_{p,b,i}\,\Delta T\\

    The properties are updated as

    .. math::

        S &\leftarrow S + \sum_i \Delta s_i\,n_i\\
        \mu_i &\leftarrow \mu_i + \Delta h_i - T\,\Delta s_i

    .. todo:
        - use independent T_0 for cp = A  + B * (T - T_0), so it is
          invariant to T_ref. Then
          H = A * (T - T_ref) + B / 2 * (T^2 - T_ref^2 - T_0 * T + T_0 * T_ref)

    The contribution model domain is defined for positive temperatures only,
    due to a logarithmic contribution to entropy.
    """

    def define(self, res):
        T, n = res["T"], res["n"]
        T_ref = res["T_ref"]
        d_T, f_T = T - T_ref, T / T_ref
        cp_a = self.par_vector("cp_a", self.species, "J/(mol*K)")
        cp_b = self.par_vector("cp_b", self.species, "J/(mol*K**2)")

        d_h = (cp_a + 0.5 * d_T * cp_b) * d_T
        d_s = (cp_a - cp_b * T_ref) * log(f_T) + cp_b * d_T
        res["S"] += d_s.T @ n
        res["mu"] += d_h - T * d_s

        self.add_bound("T", T)  # logarithm taken


class StandardState(ThermoContribution):
    """This contribution assumes the current state to be the standard state
    and tags the following properties:

    ========== =========
    Property   Origin
    ========== =========
    ``S_std``  ``S``
    ``p_std``  ``p_ref``
    ``mu_std`` ``mu``
    ========== =========

    The contribution does not support any options or parameters.
    """

    provides = ["S_std", "p_std", "mu_std"]

    def define(self, res):
        # tag current chemical potential and entropy as standard state
        # create copy, such that tagged objects will not be mutated.
        res["S_std"] = copy(res["S"])
        res["p_std"] = copy(res["p_ref"])
        res["mu_std"] = copy(res["mu"])

        self.declare_vector_keys("mu_std")  # self.species can be default


class IdealMix(ThermoContribution):
    r"""This contribution supplies the ideal mix entropy contribution,
    applicable for both liquid and gas phases. It does not require any
    parameters nor supports options.


    The contribution does not define new properties, but builds on entropy
    ``S`` (:math:`S`) and chemical potential ``mu`` (:math:`\mu`). With
    :math:`\Delta s_i := -R\,\ln n_i/N`, we update

    .. math::

        S &\leftarrow S + \sum_i \Delta s_i\,n_i\\
        \mu_i &\leftarrow \mu_i - T\,\Delta s_i

    The contribution model domain is limited to positive quantities in order to
    prevent negative arguments to the logarithmic functions. One could enable
    also the third quadrant (**all** mole numbers negative), but doing so and
    thereby allowing the solver to jump between these two nearly disconnected
    domains has proven to be challenging in terms of solver robustness.
    """

    def define(self, res):
        T, n = res["T"], res["n"]
        N = qsum(n)
        x = n / N
        gtn = R_GAS * log(x)
        res["S"] -= n.T @ gtn
        res["mu"] += T * gtn

        self.add_bound("n", n, self.species)


class GibbsIdealGas(ThermoContribution):
    r"""This contribution supplements the ideal gas entropy contribution and
    defines the volume property in Gibbs coordinates, i.e. with volume as a
    function of pressure. The main use of this class is to actually represent
    an ideal gas model, while most non-ideal gas phase models are equations
    of state and formulated in Helmholtz coordinates.

    Based on the previously provided reference pressure ``p_ref`` and
    :math:`\Delta s := -R\,\ln p/p_\mathrm{ref}`, the update is

    .. math::

        S &\leftarrow S + \Delta s\,N\\
        \mu_i &\leftarrow \mu_i - T\,\Delta s

    Volume is defined as

    .. math::

        V = \frac{N\,R\,T}{p}

    Due to the logarithmic term in entropy, the contribution model domain is
    limited to positive pressures.
    """

    provides = ["V"]

    def define(self, res):
        T, p, n, p_ref = res["T"], res["p"], res["n"], res["p_ref"]
        N = qsum(n)
        gtn = R_GAS * log(p / p_ref)

        res["S"] -= N * gtn
        res["V"] = N * R_GAS * T / p
        res["mu"] += T * gtn

        self.add_bound("p", p)


class HelmholtzIdealGas(ThermoContribution):
    r"""This contribution supplements the ideal gas entropy contribution and
    defines the pressure property in Helmholtz coordinates, i.e. with pressure
    as function of volume. This is the common base contribution for most
    equations of state. Based on the previously provided reference pressure
    ``p_ref`` and

    .. math::
        \Delta s := -R\,\ln \frac{N\,R\,T}{V\,p_\mathrm{ref}}

    the update is

    .. math::

        S &\leftarrow S + \Delta s\,N\\
        \mu_i &\leftarrow \mu_i - T\,\Delta s

    Pressure is defined as

    .. math::

        p = \frac{N\,R\,T}{V}

    Due to the logarithmic term in entropy, the contribution model domain is
    limited to positive volumes.
    """

    def define(self, res):
        T, V, n, p_ref = res["T"], res["V"], res["n"], res["p_ref"]
        N = qsum(n)
        p = N * R_GAS * T / V
        gtn = R_GAS * log(p / p_ref)

        res["S"] -= N * gtn
        res["p"] = p
        res["mu"] += T * gtn

        self.add_bound("V", V)

    def initial_state(self, state, properties):
        volume = qsum(state.mol_vector) * R_GAS * \
                 state.temperature / state.pressure
        return ([base_magnitude(state.temperature),
                 base_magnitude(volume)] +
                list(base_magnitude(state.mol_vector)))


class ConstantGibbsVolume(ThermoContribution):
    r"""This contribution defines a mixture with constant molar volumes and
    hence zero compressibility and thermal expansion. This is an easy and
    sufficient way to describe a condensed phase at moderate pressures in
    limited temperature ranges. Mixing volume is not considered either, but can
    be added via additional contributions.

    The sole parameter is the molar volume for each species:

    ========= =================== ===========
    Parameter Description         Symbol
    ========= =================== ===========
    ``v_n``   Molar volume vector :math:`v_n`
    ========= =================== ===========

    The system volume is then calculated as

    .. math::

        V = \sum_i v_{n,i}\,n_i

    For the chemical potential, this yields, given pressure ``p`` and reference
    pressure ``p_ref``:

    .. math::

        \mu_i \leftarrow \mu_i + v_{n,i}\,(p - p_\mathrm{ref})

    .. note::

        The volume definition has a slight impact on the chemical potential,
        as shown in above equation. For a typical condensed phase however, say
        water, :math:`v_n = 1.8\cdot 10^{-5}\ \mathrm{m^3/mol}`, and if
        pressures are in the range of few bars,
        :math:`\Delta \mu < 10\ \mathrm{J/mol}`. This yields a 0.4 % change in
        the activity coefficient.
    """
    provides = ["V"]

    def define(self, res):
        n, p, p_ref = res["n"], res["p"], res["p_ref"]
        v_n = self.par_vector("v_n", self.species, "m**3/mol")
        res["mu"] += v_n * (p - p_ref)
        res["V"] = v_n.T @ n


class MolecularWeight(ThermoContribution):
    """This contribution defines the molecular weights as a species vector,
    based on the underlying :class:`simu.SpeciesDefinition` definitions."""
    provides = ["mw"]

    def define(self, res):
        mw = qvertcat(*[s.molecular_weight
                        for s in self.species_definitions.values()])
        res["mw"] = mw
        self.declare_vector_keys("mw")


class ChargeBalance(ThermoContribution):
    r"""This contribution defines a charge balance residual, if there are
    charged species. A ValueError is raised if only either positive or negative
    charged species exist.

    Given the presents of opposite charged species, the residual defined, based
    on the molar quantities :math:`n_i` is

    .. math:: \sum_i n_i \cdot c_i = 0

    Here, the charges :math:`c_i` are extracted from the species
    definition.

    By default, the unit for the residual is ``e`` for states and ``e / s`` for
    flows. This can be altered by the following options:

      - The option ``state_unit`` defines the unit used for states. If
        ``flow_unit`` is not given, it will be defined as ``state_unit / s``.
      - The option ``flow_unit`` defines the unit used for flows. It does not
        impact the ``state_unit`` value.

    """

    def define(self, res):
        charges = [s.charge for s in self.species_definitions.values()]
        pos = len([c for c in charges if c > 0])
        neg = len([c for c in charges if c < 0])
        if pos == 0 and neg == 0:
            return
        if pos == 0 or neg == 0:
            raise ValueError("Charges of only one sign cannot be balanced")
        charges = qvertcat(*charges)
        residual = res["n"].T  @ charges

        state_unit = self.options.get("state_unit", "e")
        flow_unit = self.options.get("flow_unit", f"{state_unit} / s")
        try:
            self.add_residual("balance", residual, flow_unit)
        except DimensionalityError:
            self.add_residual("balance", residual, state_unit)
