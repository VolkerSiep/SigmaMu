Cubic equations of state
========================

Redlich Kwong EOS
-----------------
The base-class is specialised into two sub-classes to address in particular liquid and gas phases. As such, the base-class is disabled from direct instantiation.

Base class
..........
.. autoclass:: simu.app.thermo.contributions.cubic.rk.RedlichKwongEOS
  :show-inheritance:

Liquid root
...........
.. autoclass:: simu.app.thermo.contributions.cubic.RedlichKwongEOSLiquid
  :show-inheritance:
  :members: relax

Gas root
........
.. autoclass:: simu.app.thermo.contributions.cubic.RedlichKwongEOSGas
  :show-inheritance:
  :members: relax

Critical parameters
-------------------
.. autoclass:: simu.app.thermo.contributions.cubic.CriticalParameters
  :show-inheritance:

Non-symmetric mixing rule
-------------------------
 .. autoclass:: simu.app.thermo.contributions.cubic.NonSymmetricMixingRule
  :show-inheritance:

Linear mixing rule
-------------------------
 .. autoclass:: simu.app.thermo.contributions.cubic.LinearMixingRule
  :show-inheritance:

Redich Kwong A-function
-----------------------
 .. autoclass:: simu.app.thermo.contributions.cubic.RedlichKwongAFunction
  :show-inheritance:

Redich Kwong B-function
-----------------------
 .. autoclass:: simu.app.thermo.contributions.cubic.RedlichKwongBFunction
  :show-inheritance:

Mathias Boston-Mathias alpha-function
-------------------------------------
 .. autoclass:: simu.app.thermo.contributions.cubic.BostonMathiasAlphaFunction
  :show-inheritance:

Redlich Kwong m-factor
----------------------
 .. autoclass:: simu.app.thermo.contributions.cubic.RedlichKwongMFactor
  :show-inheritance:

Appendix
--------
.. toctree::
   :maxdepth: 1
   
   cubic/alpha_extensions
