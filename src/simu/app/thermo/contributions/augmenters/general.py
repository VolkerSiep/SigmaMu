from casadi import DM

from simu import (
    ThermoContribution, qsum, registered_contribution, Quantity,
    SpeciesDefinition, qvertcat)


@registered_contribution
class GenericProperties(ThermoContribution):
    r"""Provide basic derived thermodynamic properties:

    ========= ================================= ==============
    Property  Description                       Symbol
    ========= ================================= ==============
    ``G``     Total Gibbs free energy [J]/[W]   :math:`G`
    ``H``     Total enthalpy [J]/[W]            :math:`H`
    ``A``     Total Helmholtz energy [J]/[W]    :math:`A`
    ``U``     Total inner energy [J]/[W]        :math:`U`
    ``N``     Total moles [mol]/[mol/s]         :math:`N`
    ``M``     Total mass [kg]/[kg/s]            :math:`M`
    ``m``     Partial mass (flows) [kg]/[kg/s]  :math:`m_i`
    ``x``     Mole fractions [-]                :math:`x_i`
    ``w``     Mass fractions [-]                :math:`w_i`
    ``Mw``    Average molecular weight [kg/mol] :math:`\bar M`
    ========= ================================= ==============

    With entropy :math:`S`, chemical potential :math:`\mu_i` and molecular
    weights :math:`M_i`, it is

    .. math::
      :no-wrap:

      \begin{alignat*}{4}
        G &= \sum_i \mu_i\,n_i &\qquad H &= G + T\,S &\qquad
        A &= G - p\,V &\qquad U &= A + T\,S \\
        N &= \sum_i n_i & m_i &= n_i\,M_i &
        M &= \sum_i m_i & \bar M &= \frac{M}{N} \\
        x_i &= \frac{n_i}{N} & w_i &= \frac{m_i}{M}
      \end{alignat*}
    """

    provides = ["G", "H", "A", "U", "N", "m", "M", "x", "w", "Mw"]

    def define(self, res):
        res["G"] = res["n"].T @ res["mu"]
        st = res["T"] * res["S"]
        pv = res["p"] * res["V"]
        res["H"] = res["G"] + st
        res["A"] = res["G"] - pv
        res["U"] = res["A"] + st

        res["N"] = qsum(res["n"])
        res["m"] = res["n"] * res["mw"]
        res["M"] = qsum(res["m"])
        res["Mw"] = res["M"] / res["N"]

        res["x"] = res["n"] / res["N"]
        res["w"] = res["m"] / res["M"]

        for name in ("m", "x", "w"):
            self.declare_vector_keys(name)


@registered_contribution
class Elemental(ThermoContribution):
    r"""Provide flows or quantities per chemical element as follows:

    ========= ================================= ================
    Property  Description                       Symbol
    ========= ================================= ================
    ``n_e``   Elemental moles                   :math:`n_{e,j}`
    ``x_e``   Elemental mole fractions          :math:`x_{e,j}`
    ``N_e``   Total elemental moles             :math:`N_{e}`
    ``m_e``   Elemental masses                  :math:`m_{e,j}`
    ``w_e``   Elemental mass fractions          :math:`w_{e,j}`
    ========= ================================= ================

    Based on the parsed chemical formulae of each species, the contribution
    first determines the super-set of occurring elements
    :math:`\mathbb E_\cup = \bigcup_{i\in\mathbb S} \mathbb E_i`.
    The set is handled as a sorted list for reproducibility. Stoichiometric
    coefficients :math:`\nu_{ij}` describe the occurrence of element :math:`j`
    in species :math:`i`. Then, and with atomic weights :math:`M_j`:

    .. math::
      :nowrap:

      \begin{alignat}{3}
        n_{e,j} &= \sum_{i\in\mathbb S} \nu_{ij}\cdot n_i &\qquad
        \dot N_e & = \sum_{j\in\mathbb E_\cup} \dot n_{e,j} &\qquad
        x_e & = n_{e,j} / N_e\\
        m_{e,j} &= M_j \, n_{e,j} &&&
        w_{e,j} & =  m_{e,j} / \sum_{k\in\mathbb E_\cup} m_{e, k}
      \end{alignat}
    """
    provides = ["n_e", "x_e", "N_e", "m_e", "w_e"]

    def define(self, res):
        n = res["n"]

        # find all elements
        all_elements = set()
        for definition in self.species_definitions.values():
            all_elements |= definition.elements.keys()
        all_elements = sorted(all_elements)
        aw = SpeciesDefinition.formula_parser.atomic_weights
        w_e = qvertcat(*[aw[e] for e in all_elements])

        # create stoichiometry matrix
        stoichiometry = Quantity(DM(
            [[definition.elements.get(e, 0) for e in all_elements]
            for definition in self.species_definitions.values()]))

        res["n_e"] = stoichiometry.T @ n
        res["N_e"] = qsum(res["n_e"])
        res["x_e"] = res["n_e"] / res["N_e"]
        res["m_e"] = res["n_e"] * w_e
        res["w_e"] = res["m_e"] / qsum(res["m_e"])

        for name in ["n_e", "x_e", "m_e", "w_e"]:
            self.declare_vector_keys(name, all_elements)
