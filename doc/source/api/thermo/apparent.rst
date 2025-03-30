Thermodynamic models with apparent species
==========================================

Introduction
------------
The motivation to look closer into this topic has a slightly broader scope: Treating apparent species consistently is a use-case that creates a number of demands to the software design. In particular,

  - creating implicit equations, thus residuals, on the thermodynamic layer. These equations need to be solved with the process model's equations.
  - requiring partial second order derivatives that might be specific to the canonical variables of the underlying state function.

One further question is to which degree the concept of apparent species can be implemented generically, for instance just by defining the projection matrix.

When we introduce apparent species, we imply constraints, such as equilibrium reactions and/or charge balance being applied to the system at all times. Consequently, a projection of species is sufficient to describe the entire state vector of the thermodynamic model.

Good examples are any dissociating solutions, such as aqueous chloric acid. The true species are |H2O|, |HCl|, |H3O+|, and |Cl-|, but the equilibrium

.. math:: \mathrm{HCl} + \mathrm{H_2O} \rightleftarrows \mathrm{H_3O^+} + \mathrm{Cl^-}

is always considered established. It is convenient for most purposes to only refer to |H2O| and |HCl| as *apparent* species, while underlying properties, such as the pH value, electrical conductivity, and degree of dissociation can still be calculated based on the original species set.

Hypothetically, one could try to model water for instance as a mixture of |H2O|, |(H2O)2| and |(H2O)6| with individual standard states and consequent equilibrium constants between the entities. The client process model would not like to deal with all these species in the material balances, but only consider |H2O|, and we can even give it a new name, for instance ``water``. Such thermodynamic model is then interchangeable with a model that only considers |H2O|, but introduces interaction terms to describe the same physics.

The first part of defining quantities :math:`\hat n_i` of apparent species is the mapping. For the |HCl| example:

.. math::
   :nowrap:

    \begin{align*}
    \hat n_{\mathrm{H2O}} &:= n_{\mathrm{H_2O}} + n_{\mathrm{H_3O^+}}\\
    \hat n_{\mathrm{HCl}} &:= n_{\mathrm{HCl}} + n_{\mathrm{Cl^-}}
    \end{align*}

For that two species can replace four, we need two constraints, in this case the chemical equilibrium of dissociation and the charge balance:

.. math::
   :nowrap:

    \begin{align*}
      \mu_{\mathrm{H2O}} + \mu_{\mathrm{HCl}} &= \mu_{\mathrm{H3O^+}} + \mu_{\mathrm{Cl^-}}\\
      n_{\mathrm{H_3O^+}} &= n_{\mathrm{Cl^-}}
    \end{align*}

There is a lazy and a consistent approach to thermodynamic models with apparent species. For the lazy approach, we are basically done, only that all derived species-specific properties, such as mass flows, mole fractions and mass fractions are mapped as well, and the chemical potentials are what they are, as they are equal anyhow. The consistent approach however deserves its own section.

The consistent approach
-----------------------
As said, the client of the thermodynamic model shall not see a fundamental difference between an electrolyte model that considers |H3O+| and |Cl-| and an equally good model based on a molecular description and interaction.

Apparent heat capacity
^^^^^^^^^^^^^^^^^^^^^^

But what do these models calculate as heat capacity, defined as

.. math::  c_p := \pade{H}{T}{p, n}

The inconsistency here is what it means to keep :math:`n` constant. The molecular model clearly keeps |H2O| and |HCl| constant, while the electrolyte model has the choice to keep all internal species constant or only the apparent ones.

For the sake of consistency and physical interpretation, we must allow for the true speciation to change: A change in temperature yields a change in dissociation, which in turn, due to the heat of reaction, yields to a contribution of the observed -- or apparent -- heat capacity -- not that any process simulation software known to mankind would care about it. Now, heat capacity is only of secondary importance in process modelling, as energy balances shall always be formulated in an integrated matter as a difference between two defined states. Secondarily however, heat capacity enters for instance via the Prandtl number empirical relationships with direct impact on equipment design, and the plant manager would not ask whether you have used true or apparent species based heat capacity for the design of the new heat exchanger or reactor. Well, back to our equations ...

The above introduced mapping can more generally be written as

.. math:: \hat n = A \cdot n

Here, :math:`A` is row-deficient, and additional constraints are required to specify the part of :math:`n` that is not observable in `\hat n`. If ions are present, one of these constraints is the charge balance. We might even already define the charge (electrons) as an additional apparent species and by that enabling modelling of electrochemical processes. The remaining nullspace of :math:`A` is the reactivity matrix :math:`E` of the system:

.. math:: A \cdot E^\mathrm{T} = 0

Then, the missing constraints are:

.. math:: E \cdot \mu = 0


The total differential of enthalpy in Gibbs coordinates :math:`H(T, p, n)` is

.. math:: \mathrm{d}H = \pade{H}{T}{p, n}\, \mathrm{d}T + \pade{H}{n}{T, p} \, \mathrm{d}n

Further, the total differential of :math:`E \cdot \mu` is

.. math::

  \mathrm{d}(E \cdot \mu) =
    E \cdot \pade{\mu}{T}{p, n} \, \mathrm{d}T + E \cdot \pade{\mu}{n}{T, p} \, \mathrm{d}n = 0

Combined with the demand of constant apparent species :math:`\mathrm{d}(A\cdot n) = 0`, we have

.. math::

    E \pade{\mu}{T}{p, n} \, \mathrm{d}T +
    \begin{pmatrix} E \cdot \pade{\mu}{n}{T, p}\\ A \end{pmatrix} \cdot \mathrm{d}n = 0

This constraint between :math:`\mathrm{d}T` and :math:`\mathrm{d}n` can be substituted into above total differential of enthalpy to determine the apparent and therefore correct heat capacity:

.. math::
  :nowrap:

   \begin{align*}
   \hat c_p = \pade{H}{T}{p, \hat n}
       &=  \pade{H}{T}{p, n} - \pade{H}{n}{T, p} \,
      \begin{pmatrix} E \cdot \pade{\mu}{n}{T, p} \\ A \end{pmatrix}^{-1}\cdot
      E \cdot  \pade{\mu}{T}{p, n}\\
   &= c_p + \bar h\, \begin{pmatrix} E \cdot  \mu_n\\A \end{pmatrix}^{-1}\cdot E\cdot \bar s
   \end{align*}

In the last row, we solely replaced derivatives by their physical interpretation, using partial molar enthalpies :math:`\bar h` and entropies :math:`\bar s`. Basing the calculations on these can make the implementation less dependent on the canonical coordinate system, as specific deviations happen upstream while defining these partial properties.

Thanks to `CasADi`_, above expression can be defined in such manner that even :math:`\hat c_p` can still be derived with respect to the thermodynamic state for solving the model with a second order method.

In the same manner, all required second order properties have to be evaluated, including the thermal expansion coefficient, compressibility, and subsequent speed of sound. Already at this point, we can see a pattern, as the term

.. math::

    \pade{n}{\mu}{T, p, \hat n} :=
    \begin{pmatrix}E \cdot \pade{\mu}{n}{T, p}\\ A \end{pmatrix}^{-1}\cdot E

represents the change of speciation due to shift in chemical equilibrium, and above result takes a more generic form:

.. math::

   \hat c_p = \pade{H}{T}{p, n} - \pade{H}{n}{T, p} \cdot \pade{n}{\mu}{T, p, \hat n} \cdot \pade{\mu}{T}{p, n}

.. warning::

    This approach must fail when applied to systems in which the Gibbs phase rule prohibits coexistence of both species in a reaction, such as two solid phases of the same substance or a pure gas of a dimerizing species. In this case, :math:`\mu_n` is rank-deficient and the matrix not invertible.

    Physically, for instance the heat capacity had to approach infinity at the transition temperature.

Apparent thermal expansion coefficient
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
For the fun of it, let us calculate another property, the apparent thermal expansion coefficient:

.. math::

    \hat \varepsilon_T := \frac{1}{V}\pade{V}{T}{p, \hat n}

With the insight from above section, we can derive this easily as follows:

.. math::

    \hat \varepsilon_T = \varepsilon_T +\frac{1}{V}\,\bar v\cdot
        \pade{n}{\mu}{T, p, \hat n}\cdot \bar s

Conclusion
----------
Once the partial molar properties, including the matrix :math:`\mu_n` are calculated specifically based on the canonical variables of the state function, the apparent second order properties can most likely be expressed by means of these partial molar properties and thus be generic.

The general approach would be to treat the charge balance (in the presents of ions) as a normal balance equation, represented by a row in the projection matrix :math:`A`, and construct the equilibrium constraint matrix :math:`E` as the nullspace of :math:`A`.





.. |HCl| replace:: :math:`\mathrm{HCl}`
.. |H2O| replace:: :math:`\mathrm{H_2O}`
.. |H3O+| replace:: :math:`\mathrm{H_3O^+}`
.. |Cl-| replace:: :math:`\mathrm{Cl^-}`
.. |(H2O)2| replace:: :math:`\mathrm{(H_2O)_2}`
.. |(H2O)6| replace:: :math:`\mathrm{(H_2O)_6}`


