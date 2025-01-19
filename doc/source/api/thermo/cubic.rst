Cubic equations of state
========================

Redlich Kwong EOS
-----------------
The base-class is specialised into two sub-classes to address in particular liquid and gas phases. As such, the base-class is disabled from direct instantiation.

Base class
..........
.. autoclass:: simu.app.thermo.contributions.cubic.rk.RedlichKwongEOS
  :show-inheritance:
  :exclude-members: __init__, __new__

Liquid root
...........
.. autoclass:: simu.app.thermo.contributions.cubic.RedlichKwongEOSLiquid
  :show-inheritance:
  :exclude-members: __init__, __new__

Gas root
........
.. autoclass:: simu.app.thermo.contributions.cubic.RedlichKwongEOSGas
  :show-inheritance:
  :exclude-members: __init__, __new__

Critical parameters
-------------------
.. autoclass:: simu.app.thermo.contributions.cubic.CriticalParameters
  :show-inheritance:
  :exclude-members: __init__, __new__

Non-symmetric mixing rule
-------------------------
 .. autoclass:: simu.app.thermo.contributions.cubic.NonSymmetricMixingRule
  :show-inheritance:
  :exclude-members: __init__, __new__

Linear mixing rule
-------------------------
 .. autoclass:: simu.app.thermo.contributions.cubic.LinearMixingRule
  :show-inheritance:
  :exclude-members: __init__, __new__

Redich Kwong A-function
-----------------------
 .. autoclass:: simu.app.thermo.contributions.cubic.RedlichKwongAFunction
  :show-inheritance:
  :exclude-members: __init__, __new__

Redich Kwong B-function
-----------------------
 .. autoclass:: simu.app.thermo.contributions.cubic.RedlichKwongBFunction
  :show-inheritance:
  :exclude-members: __init__, __new__

Mathias Boston-Mathias alpha-function
-------------------------------------
 .. autoclass:: simu.app.thermo.contributions.cubic.BostonMathiasAlphaFunction
  :show-inheritance:
  :exclude-members: __init__, __new__

Redlich Kwong m-factor
----------------------
 .. autoclass:: simu.app.thermo.contributions.cubic.RedlichKwongMFactor
  :show-inheritance:
  :exclude-members: __init__, __new__

Appendix
--------
.. toctree::
   :maxdepth: 1
   
   cubic/alpha_extensions
