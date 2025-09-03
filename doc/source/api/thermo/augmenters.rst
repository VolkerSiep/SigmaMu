Augmenters
==========

Augmenters are generic :class:`simu.ThermoContribution` objects that add derived physical properties to a thermodynamic model, but without impacting its definition. These can be purely informative properties, such as *average molecular weight*, but also add new information, such as transport properties.

.. currentmodule:: simu.app.thermo.contributions.augmenters.general

General Properties
------------------
.. autoclass:: GenericProperties
  :show-inheritance:
  :exclude-members: __init__, __new__

Elemental flows and fractions
-----------------------------
.. autoclass:: Elemental
  :show-inheritance:
  :exclude-members: __init__, __new__
