Cubic equations of state
========================

RedlichKwongEOS
---------------
The base-class is specialised into two sub-classes to adress in particular liquid and gas phases.
As such, the base-class is disabled from direct instantiation.

.. autoclass:: mushell.thermo.cubic.RedlichKwongEOS

.. autoclass:: mushell.thermo.cubic.RedlichKwongEOSLiquid
  :members: relax

.. autoclass:: mushell.thermo.cubic.RedlichKwongEOSGas
  :members: relax

Linear Peneloux volume shift
----------------------------
 .. autoclass:: mushell.thermo.cubic.LinearPenelouxVolumeShift


Non-symmetric mixing rule
-------------------------
 .. autoclass:: mushell.thermo.cubic.NonSymmmetricMixingRule
