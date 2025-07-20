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
