Release notes
=============

V1.0b1
------
Changes up to this point are not documented as they are part of the initial development with no prior release.
However, the most recent changes made are

- Added Pitzer Gibbs excess model, including Pitzer-Debye-Hückel contribution, binary and ternary parameters.
- Added an extended Barin heat capacity contribution
- Added IAPWS steam/water model
- Improved performance of flowsheet solver, as it used dominant time to simplify units of measurements, which was completely unnecessary.

.. note::

    As of this point, |SigmaMu| has been used successfully in an actual project with a system size of over 1000 variables.
    The model solves in 0.3 seconds, while the preparation of the model for solving can still take 2 about seconds.
    The bottle-neck here is in the preparation of the `CasADi`_ function to calculate ten-thousands of individual properties.
    We will address this matter soon - easiest by allowing to exclude many properties that are not of interest in the
    context of process modelling results (for instance chemical potentials, standard state properties, etc.)



