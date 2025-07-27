Implementing classes
====================
This part of the documentation describes implementations of the core thermodynamic base classes, such as :class:`simu.ThermoContribution` and :class:`simu.StateDefinition`.


ThermoFrame Factory
-------------------

.. autoclass:: simu.app.RegThermoFactory
  :show-inheritance:
  :exclude-members: __init__, __new__

.. autodecorator:: simu.registered_contribution
.. autodecorator:: simu.registered_state

Pre-defined data
----------------
.. rubric:: all_states (simu.core)

This is a list of all currently imported :class:`~simu.StateDefinition` class definitions that are decorated with :func:`~simu.registered_state`. It normally contains at least the :class:`~simu.app.thermo.state.GibbsState` and the :class:`~simu.app.thermo.state.HelmholtzState` class.

.. rubric:: all_contributions (simu.core)

This is a list of all currently imported :class:`~simu.ThermoContribution` class definitions that are decorated with :func:`~simu.registered_contribution`.

.. rubric:: predefined_parameters (simu.app)

This is a :class:`~simu.ThermoParameterStore` instance with some pre-defined thermodynamic parameters.

State definitions
-----------------
.. toctree::
   :maxdepth: 1

   statedefinitions

Thermodynamic contributions
---------------------------
.. toctree::
   :maxdepth: 2

   ideal
   cubic
   iapws
   special
   augmenters
