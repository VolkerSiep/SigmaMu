Standard property names of thermodynamic model results
======================================================

Introduction
------------
:class:`~mushell.thermo.ThermoContribution` objects generate calculated properties that are eventually via the :class:`~mushell.thermo.ThermoFrame` objects exported to the client code. The following section puts some standard on the names of these properties, in particular for generic properties that can be expected to be evaluated by any thermodynamic model.

Obligatory properties
---------------------

============================== ===========================================================================
Variable name                  Description
============================== ===========================================================================
``x``                          This property is not calculated, but initially provided to the first
                               contribution. It simply contains the state as a vector with:math:`n+2`
                               elements.
``T``, (``p`` or ``V``), ``n`` Based on the state above, the first contribution defines the interpretation
                               of the state, i.e. whether it is for instance a Gibbs or a Helmholtz model.
                               (``n`` is a vector) 
``S``                          Entropy is first defined by the reference state and supplemented throughout
                               the contribution chain.
``mu``                         Chemical potential is first defined by the reference state and supplemented
                               throughout the contribution chain. (vector)
``V`` or ``p``                 The conjugated variable to the state, either volume or pressure.
============================== ===========================================================================

Optional common properties
--------------------------

Some of the above properties that are updated from one contribution to another (e.g. ``S``, ``mu``),
can be frozen in various states with below suffixes.
  
============= =========================
Variable name Description
============= =========================
``*_std``     Standard state properties
``*_im``      Ideal mix properties
``*_ig``      Ideal gas properties
``*_res``     Residual property
``*_ex``      Excess property
============= =========================

Optionally, a model can export second order properties, for which the physical interpretation is dependent on the coordinate system:

============= ======================================================
Variable name Description
============= ======================================================
``dSdT``      Temperature derivative of entropy at constant pressure
``dVdT``      Isobaric expansion
``dVdp``      Isothermal compression
``s_bar``     Partial molar entropies (vector)
``v_bar``     Partial molar volumes (vector)
============= ======================================================

.. note::

	We encourage to not publish derivable properties on this layer, as these are very numerous and not
	required for each model evaluation. These include for instance enthalpy, inner ennergy, compositional
	fractions, heat capacities, adiabatic exponent, compressibility, thermal expansitivity, Joule-Thomson
	coefficient, partial molar enthalpy and speed of sound. 
