========
Tutorial
========

Creating a thermodynamic model
==============================

First some theory to set the scene
----------------------------------

.. note::

    By reading this part, you might learn about and / or gain a more structured view about what thermodynamic models are. Most text books, specially those of older dates, like to confuse at this point with fugacities, activities, and their coefficients. If you want to learn more about those, read an old thermodynamics book.

The basis for all process simulation is a physical description of the material(s) at hand. This is the task of thermodynamic models, and while these can be quite elaborate, we start with a simpler model -- primarily for the sake of this tutorial, but also because simple models are often sufficiently accurate for a given task.

Hearing about an *ideal gas*, the famous equation

.. math::
    p\,V = N\,R\,T

will flash in bright colors within the brain of most readers, accompanied by a slight scent of vanilla -- or is that only me? Indeed, above equation is the *ideal gas law*, and by that the most primitive equation that one might call an *equation of state* (EOS). However, it does not contain any standard state information as such, for instance heat capacity.

To start with, we need a **reference state**, describing the enthalpy and entropy of the species :math:`i` at a fixed reference temperature :math:`T_\mathrm{ref}` and pressure :math:`T_\mathrm{ref}`. This is expressed as

.. math::
    \mu_{i, \mathrm{ref}} = (T_\mathrm{ref}, p_\mathrm{ref}) = \Delta_f h_i - T\,s_{0,i}

Here, :math:`\Delta_f h_i` is the enthalpy of formation, and :math:`s_{0,i}` the standard entropy of species :math:`i`. They are *thermodynamic parameters* for which we later need to find values. Above equation yields the reference state chemical potential. Typically, we like to define :math:`T_\mathrm{ref} = 25 ^\circ`\ C and :math:`p_\mathrm{ref} = 1` bar, but this is really arbitrary, as long as we stay consistent.

Next, we describe by the so-called **standard state** how the substance reacts to temperature changes at constant pressure (still :math:`p_\mathrm{ref}`) and *standard state conditions*, which can be for instance as a pure substance. This is done by defining a standard state **heat capacity** :math:`c_{p, i}(T, p_\mathrm{ref})`, being the enthalpy increase of the system experiencing a 1 K temperature increase at standard state conditions. We can describe this as

.. math::
    \mu_i^\standard(T, p_\mathrm{ref}) = \mu_{i, \mathrm{ref}} (T_\mathrm{ref}, p_\mathrm{ref})
      + \int\limits_{T_\mathrm{ref}}^T c_{p, i}\,\mathrm{d}\tau
      - T\,\int\limits_{T_\mathrm{ref}}^T \frac{c_{p, i}}{\tau}\,\mathrm{d}\tau

By now, we already have a sufficient complete description of for instance *pure solids* and *pure liquids* under moderate pressures. Parameters can be found for instance in :cite:p:`Wagman_1982`. By adding two more contributions, we arrive at the ideal gas. These are the **ideal mix** (im) and the **ideal gas** (ig) contributions:

.. math::
    \begin{align*}
      \mu_i^\mathrm{im}(T, p_\mathrm{ref}) &= \mu_i^\standard(T, p_\mathrm{ref}) + R\,T\,\ln x_i\\
      \mu_i^\mathrm{ig}(T, p) &= \mu_i^\mathrm{im}(T, p_\mathrm{ref}) + R\,T\,\ln \frac{p}{p_\mathrm{ref}}
    \end{align*}

For a pure substance, we can even leave out the ideal mix term, as :math:`\ln x_i` is zero because of the mole fraction :math:`x_i` being unity. As a last step, we Euler-integrate the chemical potential to the Gibbs free energy via :math:`G = \sum_i n_i\,\mu_i`:

.. math::
    G = N\,\mu_i^\standard(T, p_\mathrm{ref}) + N\,R\,T\,\ln \frac{p}{p_\mathrm{ref}}

This last step was merely to show what happens when we derive the Gibbs energy with respect to pressure to obtain the volume of the system:

.. math::
    \left .\frac{\partial G}{\partial p}\right |_{T, N} = V = \frac{N\,R\,T}{p}

Surprise! We are back at the ideal gas law from above. This is simply because we silently integrated the ideal volume over pressure when throwing in the ideal gas contribution.

How to do this in SiMu?
-----------------------
As indicated above, thermodynamic models can easily be chunked into contributions, and these contributions are often additive, in this case exclusively. All above introduced contributions are already defined in ``SiMu``. They take care of defining the required parameters and adding the contributions up in order to form a state function like the Gibbs free energy. As a first step, we create a :class:`simu.ThermoFactory` object, register the Gibbs state and all thermodynamic model contributions that are available:

.. literalinclude:: ../examples/ideal_gas.py
   :language: python
   :linenos:
   :lines: 1-6

Next, a :class:`simu.ThermoFrame` can be constructed, stacking the contributions as discussed in the previous section, based on a Gibbs thermodynamic state, that is using the state function :math:`G(T, p, n)`:

.. literalinclude:: ../examples/ideal_gas.py
   :language: python
   :linenos:
   :lineno-start: 8
   :lines: 8-19

Here we decided to create an ideal gas model for methane, and used the following contributions:

 - :class:`simu.app.thermo.contributions.H0S0ReferenceState`
 - :class:`simu.app.thermo.contributions.LinearHeatCapacity`
 - :class:`simu.app.thermo.contributions.StandardState`
 - :class:`simu.app.thermo.contributions.IdealMix`
 - :class:`simu.app.thermo.contributions.GibbsIdealGas`

Note that we could have left out the :class:`StandardState <simu.app.thermo.contributions.StandardState>` class, as it only preserves the standard state properties and does no own calculations. Only if we later would like to introduce the concepts of fugacity :math:`f_i` and activity :math:`a_i`, we would require to access these properties. However, this is not interesting for an ideal gas, where :math:`f_i = p_i` and :math:`a_i = x_i`.

The :class:`simu.SpeciesDefinition` defines in this simple case some basic generic properties based on the chemical formula, that is molecular weight, elementary composition and electrical charge. For fun, you can analyse some formulae, for instance

>>> from simu import SpeciesDefinition
>>> print(SpeciesDefinition("Ca(CH3-CH2-COO)2").elements)
{'Ca': 1, 'C': 6, 'H': 10, 'O': 4}
>>> print(f"{SpeciesDefinition('KMnO4').molecular_weight:~.3f}")
158.032 g / mol

However, now we have our frame object and can ask it for the required thermodynamic parameters:

.. literalinclude:: ../examples/ideal_gas.py
   :language: python
   :linenos:
   :lineno-start: 21
   :lines: 21-22

We receive the following structure, indicating the required physical dimension as a place-holder for each parameter to be provided::

    {'H0S0ReferenceState': {'T_ref': 'K',
                            'dh_form': {'Methane': 'kg * m ** 2 / mol / s ** 2'},
                            'p_ref': 'kg / m / s ** 2',
                            's_0': {'Methane': 'kg * m ** 2 / K / mol / s ** 2'}},
     'LinearHeatCapacity': {'cp_a': {'Methane': 'kg * m ** 2 / K / mol / s ** 2'},
                            'cp_b': {'Methane': 'kg * m ** 2 / K ** 2 / mol / s ** 2'}}}

It is to be noted that the actual parameters can be given in different units, as long as they are compatible -- we come back to this.

Likewise, the frame can be queried for the properties that are calculated:

.. literalinclude:: ../examples/ideal_gas.py
   :language: python
   :linenos:
   :lineno-start: 24
   :lines: 24

This yields::

    {'S': 'kg * m ** 2 / K / s ** 2',
     'S_std': 'kg * m ** 2 / K / s ** 2',
     'T': 'K',
     'T_ref': 'K',
     'V': 'J * m * s ** 2 / kg',
     '_state': 'dimless',
     'mu': 'kg * m ** 2 / mol / s ** 2',
     'mu_std': 'kg * m ** 2 / mol / s ** 2',
     'mw': 'g / mol',
     'n': 'mol',
     'p': 'kg / m / s ** 2',
     'p_ref': 'kg / m / s ** 2',
     'p_std': 'kg / m / s ** 2'}


.. todo::
  - put in some parameters for methane and calculate properties
  - can I add calculation of density? How did I plan this again!?!??