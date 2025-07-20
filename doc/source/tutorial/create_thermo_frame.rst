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

will flash in bright colors within the brain of most readers, accompanied by a slight scent of vanilla -- or is that only me? Indeed, above equation is the *ideal gas law*, and by that the most primitive equation that one might call an *equation of state* (EOS). However, it does not contain any standard state information as such, that is enthalpy, entropy, and heat capacity.

To start with, we need a **reference state**, describing the enthalpy and entropy of the species :math:`i` at a fixed reference temperature :math:`T_\mathrm{ref}` and pressure :math:`p_\mathrm{ref}`. This is expressed as

.. math::
    \mu_{i, \mathrm{ref}} = (T_\mathrm{ref}, p_\mathrm{ref}) = \Delta_f h_i - T\,s_{0,i}

Here, :math:`\Delta_f h_i` is the enthalpy of formation, and :math:`s_{0,i}` the standard entropy of species :math:`i`. These are *thermodynamic parameters* for which we later need to find values. Above equation yields the reference state chemical potential. Typically, we like to define :math:`T_\mathrm{ref} = 25 ^\circ`\ C and :math:`p_\mathrm{ref} = 1` bar, but this is really arbitrary, as long as we stay consistent.

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

For a pure substance, we can even leave out the ideal mix term, because :math:`\ln x_i = 0`, given :math:`x_i = 1`. As a last step, we Euler-integrate the chemical potential to the Gibbs free energy via :math:`G = \sum_i n_i\,\mu_i`:

.. math::
    G = N\,\mu_i^\standard(T, p_\mathrm{ref}) + N\,R\,T\,\ln \frac{p}{p_\mathrm{ref}}

This last step was merely to show what happens when we derive the Gibbs energy with respect to pressure to obtain the volume of the system:

.. math::
    \left .\frac{\partial G}{\partial p}\right |_{T, N} = V = \frac{N\,R\,T}{p}

Surprise! We are back at the ideal gas law from above. This is simply because we silently integrated the ideal volume over pressure when throwing in the ideal gas contribution.

How to do this in SiMu?
-----------------------
As indicated above, thermodynamic models can easily be chunked into contributions, and these contributions are often additive, in this case exclusively. All above introduced contributions are already defined in ``SigmaMu``. They take care of defining the required parameters and adding the contributions up in order to form a state function like the Gibbs free energy. As a first step, we do some necessary imports for this session and create a :class:`simu.app.RegThermoFactory` object, which has already registered all available states and contributions:

.. exampleinclude:: ideal_gas.py
   :language: python
   :linenos:
   :lines: 1-3

Note that generic entities are imported directly from the ``simu`` root module, while specific objects, such as :class:`simu.app.RegThermoFactory` are provided by the ``simu.app`` module.

Next, a :class:`simu.ThermoFrame` can be constructed, stacking the contributions as discussed in the previous section, based on a Gibbs thermodynamic state, that is using the state function :math:`G(T, p, n)`:

.. exampleinclude:: ideal_gas.py
   :language: python
   :linenos:
   :lineno-start: 4
   :lines: 5-17

Here we decided to create an ideal gas model for methane, and used the following contributions (don't be afraid to follow the links to see their parameters, calculated properties, and mathematical formulations:

 - :class:`~simu.app.thermo.contributions.basic.H0S0ReferenceState`
 - :class:`~simu.app.thermo.contributions.basic.LinearHeatCapacity`
 - :class:`~simu.app.thermo.contributions.basic.IdealMix`
 - :class:`~simu.app.thermo.contributions.basic.GibbsIdealGas`

The :class:`simu.SpeciesDefinition` defines in this simple case some basic generic properties based on the chemical formula, that is molecular weight, elementary composition and electrical charge. For fun, you can analyse some formulae, for instance

>>> from simu import SpeciesDefinition
>>> print(SpeciesDefinition("Ca(CH3-CH2-COO)2").elements)
{'Ca': 1, 'C': 6, 'H': 10, 'O': 4}
>>> print(f"{SpeciesDefinition('KMnO4').molecular_weight:~.3f}")
158.032 g / mol

However, now we have our frame object and can ask it for the required thermodynamic parameters:

.. testsetup::

    >>> from pprint import pprint
    >>> from simu.examples.ideal_gas import frame, parameters

>>> pprint(frame.parameter_structure, width=90)
{'H0S0ReferenceState': {'T_ref': 'K',
                        'dh_form': {'Methane': 'J / mol'},
                        'p_ref': 'Pa',
                        's_0': {'Methane': 'J / K / mol'}},
 'LinearHeatCapacity': {'cp_a': {'Methane': 'J / K / mol'},
                        'cp_b': {'Methane': 'J / K ** 2 / mol'}}}

We receive the structure above, indicating the required physical dimension as a place-holder for each parameter to be provided. It is to be noted that the actual parameters can be given in different units, as long as they are compatible -- we come back to this.

..  note::

    To start with, the listed units will be based on the SI units (more precise: internal `Pint`_  units), but we apply an algorithm that introduces derived SI units, such as ``Pa``, ``J`` and ``W`` if this yields simpler expressions.

Likewise, the frame can be queried for the properties that will be calculated:

>>> pprint(frame.property_structure)
{'bounds': {'GibbsIdealGas': {'p': 'Pa'},
            'IdealMix': {'n': {'Methane': 'mol'}},
            'LinearHeatCapacity': {'T': 'K'}},
 'props': {'S': 'J / K',
           'T': 'K',
           'T_ref': 'K',
           'V': 'm ** 3',
           '_state': '',
           'mu': 'J / mol',
           'n': 'mol',
           'p': 'Pa',
           'p_ref': 'Pa'}}

Here, the bounds are variables that the thermodynamic model requires to be positive.

Getting back to the parameter structure, let us fill in some values and convert the structure into a dictionary of quantities:

.. exampleinclude:: ideal_gas.py
   :language: python
   :linenos:
   :lineno-start: 17
   :lines: 19-30

Here, we use the function :func:`simu.parse_quantities_in_struct` to convert the dictionary values into quantities.

We can call the model at this point, just be aware of an important design feature in ``SigmaMu``:

.. important::

    Thermodynamic states (for instance :math:`(T, p, \vec n)` or :math:`(T, V, \vec n)` are considered as plain arrays of their values in SI units. As an example, ``[373.15, 1e5, 1.0]`` represents a thermodynamic state for our model, describing a system of 1 mol methane at 1 bar pressure and 100 |degC|. This is because ``SigmaMu`` uses thermodynamic states as independent variables during numerical solving. Only the thermodynamic models themselves are to interpret the array elements.

We show in a moment how to create such states from physical quantities, but for now, the following code computes our ideal gas model:

>>> result = frame([400, 1e5, 1.0], parameters)["props"]
>>> for key, value in result.items():
...     print(f"{key}: {value:.5g~}")
S: 199.86 J / K
T: 400 K
T_ref: 298.15 K
V: 0.033258 m ** 3
_state: [400 1e+05 1]
mu: -1.5092e+05 J / mol
n: 1 mol
p: 1e+05 Pa
p_ref: 1e+05 Pa

As this is a pure species ideal gas, not much exciting is going on, but we have a chemical potential, entropy and volume. *Wait, where is enthalpy, heat capacity, density, compressibility, and all this interesting stuff?* -- We will come back to that. First let's look at a better way to define the initial state:

.. exampleinclude:: ideal_gas.py
   :language: python
   :linenos:
   :lineno-start: 36
   :lines: 36-38

Here, we use the :class:`simu.InitialState` helper class to generate a tuple of temperature, pressure, and molar quantities ``tpn``. Subsequently, the model (:class:`simu.ThermoFrame` object) takes this definition and turns it into a valid state for itself. In this case it trivially returns the SI values of the given state, but a proper initialization is and must be performed by the relevant thermodynamic contributions, if we define a model Helmholtz coordinates (an equation of state) or more exotic models, that could use entropy or enthalpy as free variables.

Summary / Outlook
-----------------
- We defined an ideal gas thermodynamic model based on already defined thermodynamic contributions.
- Parameters and calculated properties are represented by ``pint`` quantities in dictionaries.
- Thermodynamic states are pure numerical objects (arrays), only to be interpreted by the thermodynamic models themselves.

Parts of the above might seem cumbersome for this *application*, but we are off for bigger things than calculating single species ideal gas properties. Later, we will among other things do the following:

- Collect thermodynamic parameters in :class:`simu.ThermoParameterStore` objects that can be shared among multiple models of different types.
- Add the calculation of more (derived) physical properties in the thermodynamic models.
- Create :class:`simu.MaterialDefinition` class objects as a glue towards the actual purpose: process modelling -- finally.

.. important::
    The entire setup of thermodynamic models including material definitions is common for all projects dealing with the same type of process. This can be developed and maintained even in a separate repository, available and encapsulated within the authoring organization. This approach assures

     - protection of intellectual property by limited access rights
     - separation of thermodynamic modelling from process modelling, likely involving different teams and experts.
     - an efficient way to provide updates to derived work, as the thermodynamic models might be improved over time.

