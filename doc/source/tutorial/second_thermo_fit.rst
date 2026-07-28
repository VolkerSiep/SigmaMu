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

Admittingly, the polar parameter for in line 45 (:math:`\eta = -10.2726...`) is somewhat high. The model probably fails completely at elevated temperatures far above 150 |degC|. However, the purpose of this model is to serve as a basis for electrolyte systems with a common temperature range between 0 and 150 |degC|.

.. todo::

  - Find data for e.g. aqueous NaCl (freezing temperature and vapour pressure as function of temperature).
  - find standard state data in Wagman (don't need to use data that challenges that).
  - Use Pitzer model.
  - Fit at least binary parameters to match data.
  - Show evaluation of original (no interaction) and fitted model
