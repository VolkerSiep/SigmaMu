Cubic equations of state
========================

Redlich Kwong EOS
-----------------
The base-class is specialised into two sub-classes to adress in particular liquid and gas phases.
As such, the base-class is disabled from direct instantiation.

Base class
..........
.. autoclass:: simu.thermo.cubic.rk.RedlichKwongEOS

Liquid root
...........
.. autoclass:: simu.thermo.RedlichKwongEOSLiquid
  :members: relax

Gas root
........
.. autoclass:: simu.thermo.RedlichKwongEOSGas
  :members: relax

Critical parameters
-------------------
.. autoclass:: simu.thermo.CriticalParameters

Non-symmetric mixing rule
-------------------------
 .. autoclass:: simu.thermo.NonSymmetricMixingRule

Linear mixing rule
-------------------------
 .. autoclass:: simu.thermo.LinearMixingRule

Redich Kwong A-function
-----------------------
 .. autoclass:: simu.thermo.RedlichKwongAFunction

Redich Kwong B-function
-----------------------
 .. autoclass:: simu.thermo.RedlichKwongBFunction

Mathias Boston-Mathias alpha-function
-------------------------------------
 .. autoclass:: simu.thermo.BostonMathiasAlphaFunction

Redlich Kwong m-factor
----------------------
 .. autoclass:: simu.thermo.RedlichKwongMFactor


Appendix
--------

.. toctree::
   :maxdepth: 1
   
   cubic/alpha_extensions
