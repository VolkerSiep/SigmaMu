Cubic equations of state
========================

Redlich Kwong EOS
-----------------
The base-class is specialised into two sub-classes to adress in particular liquid and gas phases.
As such, the base-class is disabled from direct instantiation.

Base class
..........
.. autoclass:: simu.thermo.cubic.rk.RedlichKwongEOS
  :members: name, category, requires
  :undoc-members: 

Liquid root
...........
.. autoclass:: simu.thermo.cubic.rk.RedlichKwongEOSLiquid
  :members: relax 

Gas root
........
.. autoclass:: simu.thermo.cubic.rk.RedlichKwongEOSGas
  :members: relax

Linear Peneloux volume shift
----------------------------
 .. autoclass:: simu.thermo.cubic.LinearPenelouxVolumeShift
  :members: name, category, requires
  :undoc-members: 

Non-symmetric mixing rule
-------------------------
 .. autoclass:: simu.thermo.cubic.NonSymmmetricMixingRule
  :members: name, category, requires
  :undoc-members:  

Redich Kwong A-function
-----------------------
 .. autoclass:: simu.thermo.cubic.rk.RedlichKwongAFunction
  :members: name, category, requires
  :undoc-members: 

Redich Kwong B-function
-----------------------
 .. autoclass:: simu.thermo.cubic.rk.RedlichKwongBFunction
  :members: name, category, requires
  :undoc-members: 

Mathias Boston-Mathias alpha-function
-------------------------------------
 .. autoclass:: simu.thermo.cubic.BostonMathiasAlphaFunction
  :members: name, category, requires

Redlich Kwong m-factor
----------------------
 .. autoclass:: simu.thermo.cubic.rk.RedlichKwongMFactor
  :members: name, category, requires



Appendix
--------

.. toctree::
   :maxdepth: 1
   
   cubic/alpha_extensions
