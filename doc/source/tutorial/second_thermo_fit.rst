====================================================
A slightly more advanced thermodynamic parameter fit
====================================================
The previous example was a master-piece example of using a sledgehammer to crack a nut.

This example is slightly more complex and realistic, namely to model the vapour pressure of water over an aqueous solution of ``NaCl``. For this, we extract some values from :cite:`Washburn_1928`.

To fit this task into bite-size, let us do the following:

- The pure water (ideal) liquid model, consisting of the standard state and molar volume parameterization is developed as a separate example ``h2o_fit`` in the ``examples`` folder of this repository, using calculated data of the IAPWS EoS as reference. Here, the volume is modelled solely as a function of temperature.
- The SRK EoS, used for steam, is as well parameterized in ``h2o_fit`` based on the IAPWS model. Here, low pressure data is used to fit the standard state parameters, and elevated pressure data to fit the polar parameter of the Boston-Mathias alpha function.
- No chemical reactions, calorimetric, nor volumetric properties regarding ``Na+`` and ``Cl-`` are relevant for this example. Hence the standard state and molar volume parameters can be left to zero.

Pre-done work: Fit of pure water model for electrolyte applications
===================================================================

The VLE fit to the IAPWS model is shown in the figure below. Deviations of condensate properties are shown on the left, while deviations for steam properties are shown on the right.

.. image:: figures/h2o_vle_fit.png
    :align: center

The dashed lines represent the model with a fit standard state, but without volume parameters (for condensate) and the SRK polar parameter (for steam).

- The first row describes the deviation from equilibrium, meaning that condensate concentrations are reproduced with 0.15 % accuracy before the fit, and about 0.01 % after the fit. For steam, the model was improved from 1 % to 0.4 % accuracy.
- The second row describes the deviation in enthalpy, normalized with :math:`R\,T`. The deviation of condensate is improved from 0.01 to 0.003, that of steam from 0.06 to about 0.005.
- The third row describes the deviation in volume. The initial parameterization :math:`v_{\rm H_2O} = 18` cm3/mol results in up to 10 % deviation, while the second degree polynomial fit reduces this deviation to about 0.2 %. The steam densities are improved from 1 % to 0.4 % deviation.

The final parameter set is the following (stored in ``examples/h2o_fit/parameters_final.yml``):

.. exampleinclude:: h2o_fit/parameters_final.yml
   :language: yaml
   :linenos:

Admittingly, the polar parameter for in line 45 (:math:`\eta = -10.2726...`) is somewhat high, and the predictions are probably completely off at elevated temperatures far above 150 |degC|. However, the purpose of this model is to serve as a basis for electrolyte systems with a common temperature range between 0 and 150 |degC|.

Thermodynamic model
===================
The thermodynamic structure is defined in file ``examples/nacl_parameter_fit/thermo_config.yml``, starting with the species definition and the declaration of phases

.. exampleinclude:: nacl_parameter_fit/thermo_config.yml
   :language: yaml
   :lines: 1-5
   :linenos:

The ``phases`` definition is used in a moment to associate the actual model structure. For the **liquid phase**:

.. exampleinclude:: nacl_parameter_fit/thermo_config.yml
   :language: yaml
   :lines: 7-27
   :lineno-start: 7
   :linenos:

On top of the standard state and ideal mix contributions and the polynomial volume, the following contributions define the Pitzer model:

:class:`~simu.app.thermo.contributions.electrolytes.pitzer.ElectrolyteBasics`
  Defines some basic electrolyte properties, such as ionic strength and molality.
:class:`~simu.app.thermo.contributions.basic.ChargeBalance`
  Defines a constraint for each state, enforcing electro-neutrality.
:class:`~simu.app.thermo.contributions.electrolytes.pitzer.PitzerDebyeHueckel`
  Defines the long range interaction according to the extended Debye-Hückel model
:class:`~simu.app.thermo.contributions.electrolytes.pitzer.PitzerBinaryInteraction`
  Defines the short range binary interaction, specifically just linear in temperature and without the dependency of ionic strengths.

Note that none of these contributions alters the calculated properties of pure water.

The **gas phase** uses the Soave-Redlich-Kwong EoS with the Boston-Mathias :math:`\alpha`-function:

.. exampleinclude:: nacl_parameter_fit/thermo_config.yml
   :language: yaml
   :lines: 28-49
   :lineno-start: 28
   :linenos:

Given that only pure water is considered for the gas phase, the mixing rules are not having any effect.

The **default parameters** for Na+ and Cl- for the liquid phase are defined in the same file:

.. exampleinclude:: nacl_parameter_fit/thermo_config.yml
   :language: yaml
   :lines: 51-
   :lineno-start: 51
   :linenos:

For the purpose of this example, none of the standard state parameters require to be adjusted, as they do not impact the saturation pressure of water. The Debye-Hückel parameters are fixed for water as a solvent. The only two parameters to be fit are the interaction parameters (line 81 and line 84).

The module ``examples/nacl_parameter_fit/thermo.py`` creates the **materials** for the gas and liquid phase based on the above definitions.
It reads the definitions from ``examples/nacl_parameter_fit/thermo_config.yml`` and augments the parameters with those defined in ``examples/h2o_fit/parameters_final.yml``.

Setup for fitting aqueous NaCl vapour pressure data
===================================================
In :cite:`Washburn_1928`, page 370, the vapour pressure of aqueous NaCl solution is tabulated for a range of temperatures and NaCl weight fractions. The first step is to convert a subset of the data into our ``yaml`` format:

.. exampleinclude:: nacl_parameter_fit/Washburn_1928_all.yml
   :language: yaml
   :linenos:
   :lines: 1-4, 16-26, 49-59, 82-92, 117-128
