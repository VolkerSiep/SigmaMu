=====================================
The first thermodynamic parameter fit
=====================================

What are thermodynamic parameters?
==================================
``SigmaMu`` distinguishes thermodynamic parameters from process parameters. A thermodynamic parameter describes the behaviour of a material as a function of only its thermodynamic state. It is independent of the context, such as equipment geometry. In ``SigmaMu`` thermodynamic parameters are defined in :class:`simu.ThermoContribution` instances. Typical
examples for thermodynamic parameters are binary interaction parameters, coefficients of empirical transport property expressions, and kinetic constants to determine reaction rates.

In contrast, process model parameters are defined in :class:`simu.Model` implementations and contain equipment-dependent parameters, such as heat exchanger surfaces or efficiencies. Specified temperatures, pressures and flows are still process model parameters, as they are determining the state of the material rather than describing material properties.

Introducing the example
=======================
The detailed and the mathematical description of thermodynamic data fits is described here: :ref:`fit-of-thermodynamic-parameters`.
The short version: Let's find values of thermodynamic parameters that let our models describe the available data best.


The following is a simplest possible example to demonstrate thermodynamic parameter fitting. The demonstrated code will seem to be -- and actually be -- overkill for this case. However, this demonstrates the same approach as used for fitting interaction parameters in a modern thermodynamic model for a multi-component mixture.

Background
----------
.. |alpha-tin| replace:: :math:`\alpha`-tin
.. |beta-tin| replace:: :math:`\beta`-tin

One of the likely simplest examples is the transition between grey |alpha-tin| and white |beta-tin|.
|alpha-tin| is a brittle non-metalic form, stable below 13.2 |degC|, while |beta-tin| is metallic and has a body-centered tetragonal crystal structure :cite:p:`Wikipedia_tinn_2026`. From :cite:`Wagman_1982`, standard state data is available:

============ ====================== ============= ===============
Modification :math:`\Delta_f h_i^0` :math:`s_i^0` :math:`c_{p,i}`
             [kJ/mol]               [J/(mol K)]   [J/(mol K)]
============ ====================== ============= ===============
|alpha-tin|  -2.090                 44.14         25.77
|beta-tin|   0                      51.55         26.99
============ ====================== ============= ===============

Considering pure solid phases and approximately constant heat capacities, this is a complete description of the equilibrium properties, excluding volumetric properties. We can use the :class:`~simu.app.thermo.contributions.basic.H0S0ReferenceState` and the :class:`~simu.app.thermo.contributions.basic.LinearHeatCapacity` contributions for a sufficient description.

For the aspiring tin specialist, there is of course more literature :cite:p:`Khvan_2019` with more detailed insight.

Simulation of transition temperature
------------------------------------
Note again that the following simulation is an overkill version of solving for the temperature at which the chemical potentials of both tin forms are equal:

.. math::

    \mu_\alpha = \mu_\beta\quad\text{with}\quad
    \mu_i = \Delta_f h_i^0 - T\,s_i^0 + c_{p,i}\,\left (
        T - T^\mathrm{ref} - T^\mathrm{ref}\,\ln \frac{T}{T^\mathrm{ref}}
      \right )


The first step is to simply simulate the transition temperature based on above parameters (``wagman_tin.yml``):

.. exampleinclude:: tin_parameter_fit/wagman_tin.yml
   :language: yaml
   :linenos:

Note that ``cp_b`` is left zero for both forms, as no data is readily available, and only temperatures close to reference temperature are relevant. Next, the material definition is constructed from the model configuration, the species definitions, and the parameters (``thermo.py``):

.. exampleinclude:: tin_parameter_fit/thermo.py
   :language: python
   :linenos:

There are a couple of aspects worth being mentioned:

  - Normally, we would not care about the identifier of the thermodynamic source once it is added to the store. In this case however, we are interested in changing parameters in that source, so we need ``SOURCE_ID`` for later.
  - The configuration (line 12-15) contains no calculation of volumetric properties. Volume itself is simply not defined. We could easily add for instance the :class:`~simu.app.thermo.contributions.basic.ConstantGibbsVolume` contribution and assign molar volumes ``v_n`` to each form, but we do not need to for this case.
  - The species definition (line 16-17) is only concerned about the atomic composition and does not distinguish |alpha-tin| from |beta-tin|. This case is however a good example for the fact that we still can define such species multiple times and assign different thermodynamic parameters to it.

Based on the material definition, we can create a simple model and simulate the transition temperature:

.. exampleinclude:: tin_parameter_fit/simulation.py
   :language: python
   :linenos:

.. note::

    The result of 8.83 |degC| is very arguably different from the cited 13.2 |degC|. This is due to a fundamental thermodynamic principle: "*You can't always please everybody!*" The tabulated data is likely based on calorimetric measurements independently for both |alpha-tin| and |beta-tin|, and in its compilation not being constrained to reproduce the transition temperature.

Some explanation of above code:

  - Dealing with pure solid phases without mixing effects, the absolut quantities of each species are irrelevant. The example just specifies 1 mol of each form to be present in the system (lines 7, 15-16).
  - The model is pressure-independent, but pressure is still a state variable that needs to be constrained. The easiest is to simply specify the pressure value (lines 8, 17). If we had defined a volume term, any departure from reference pressure would give a contribution to the chemical equilibrium and hence stability.
  - Instead of specifying temperature, the model constrains :math:`\mu_\alpha = \mu_\beta` (line 18) and hence back-calculates temperature.
  - For the next step, the parameter ``T_measured`` is defined (line 9), as it enters the calculation of the penalty term ``dT_norm`` (lines 10, 29).


Parameter fit
=============
A parameter fit is configured via a data structure (``thermo_fig_definition.yml``):

.. exampleinclude:: tin_parameter_fit/thermo_fit_definition.yml
   :language: yaml
   :linenos:

There is only one data point in this case, stating that the transition temperature is at 12.3 |degC| (line 7).
To not blow this project up more than it already is, we reuse the previously developed model, which we will soon register under the given name ``transition_model`` in the :class:`~simu.ThermoFitSolver` class (line 12).

.. note::

    In more realistic cases with many data sets and samples per data set, it is beneficial for robustness and performance to create a new model that specifies the measured temperature and returns a normalized chemical potential difference, such as :math:`\Delta \mu_{\alpha,\beta} / (R\,T)`. The resulting equation system is then linear and robust.

    This approach can be used whenever equilibria are to be evaluated, as an alterative to solving the equilibrium constraints at hand.

The ``data_to_model`` section says: "Use the value of the ``T_trans`` column in the data set as the parameter ``T_measured`` in the model" (line 14). The data fit will pick up ``DT_norm`` (line 16) as part of the objective to minimize the square of all penalty contributions.

*Well, normally we would need more data samples to generate a real minimization problem. This case degrades to a square system with a resulting zero penalty, reproducing the measured temperature exactly. But then again: This is a minimal example, not made to impress anybody.*

Finally, the standard entropy of |alpha-tin| is to be parameterized (line 19-20).

.. exampleinclude:: tin_parameter_fit/parameter_fit.py
   :language: python
   :linenos: