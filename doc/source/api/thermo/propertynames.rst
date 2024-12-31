.. _standard property names:

Standard property names of thermodynamic model results
======================================================

Introduction
------------
:class:`~simu.ThermoContribution` objects generate calculated properties that are eventually via the :class:`~simu.ThermoFrame` objects exported to the client code. The following section puts some standard on the names of these properties, in particular for generic properties that can be expected to be evaluated by any thermodynamic model.

Generally, properties starting with an underscore ``_`` are not exposed in the material objects and hence the process modelling context.

Standardized properties
-----------------------
The following properties shall be published under the listed names, such that downstream contributions and materials can rely on them:

============== ===========================================================================
Variable name  Description
============== ===========================================================================
``_state``     The canonical state of the thermodynamic model, as it may be required for
               derivative calculations in constrained systems.
``T``          Temperature
``p``          Pressure
``n``          Vector of molar quantities
``S``          Entropy (as the extensive quantity)
``mu``         Chemical potential vector
============== ===========================================================================

Standard property suffixes
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
