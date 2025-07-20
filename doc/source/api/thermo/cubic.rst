Cubic equations of state
========================

Redlich Kwong EOS
-----------------
The base-class is specialised into two sub-classes to address in particular liquid and gas phases. As such, the base-class is disabled from direct instantiation.


Base class
..........
.. currentmodule:: simu.app.thermo.contributions.cubic.rk
.. autoclass:: RedlichKwongEOS
  :show-inheritance:
  :exclude-members: __init__, __new__

Liquid root
...........
.. autoclass:: RedlichKwongEOSLiquid
  :show-inheritance:
  :exclude-members: __init__, __new__

Gas root
........
.. autoclass:: RedlichKwongEOSGas
  :show-inheritance:
  :exclude-members: __init__, __new__

A-function
..........
 .. autoclass:: RedlichKwongAFunction
  :show-inheritance:
  :exclude-members: __init__, __new__

B-function
..........
 .. autoclass:: RedlichKwongBFunction
  :show-inheritance:
  :exclude-members: __init__, __new__

m-factor
........
 .. autoclass:: RedlichKwongMFactor
  :show-inheritance:
  :exclude-members: __init__, __new__


Critical parameters
-------------------
.. currentmodule:: simu.app.thermo.contributions.cubic.core
.. autoclass:: CriticalParameters
  :show-inheritance:
  :exclude-members: __init__, __new__

Non-symmetric mixing rule
-------------------------
 .. autoclass:: NonSymmetricMixingRule
  :show-inheritance:
  :exclude-members: __init__, __new__

Linear mixing rule
-------------------------
 .. autoclass:: LinearMixingRule
  :show-inheritance:
  :exclude-members: __init__, __new__

Mathias Boston-Mathias alpha-function
-------------------------------------
 .. autoclass:: BostonMathiasAlphaFunction
  :show-inheritance:
  :exclude-members: __init__, __new__

Appendix
--------
.. toctree::
   :maxdepth: 1
   
   cubic/alpha_extensions
