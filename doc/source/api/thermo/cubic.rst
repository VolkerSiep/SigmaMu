Cubic equations of state
========================

Redlich Kwong EOS
-----------------
The base-class is specialised into two sub-classes to adress in particular liquid and gas phases.
As such, the base-class is disabled from direct instantiation.

Base class
..........
.. autoclass:: mushell.thermo.cubic.rk.RedlichKwongEOS
  :members: name, category, requires
  :undoc-members: 

Liquid root
...........
.. autoclass:: mushell.thermo.cubic.rk.RedlichKwongEOSLiquid
  :members: relax 

Gas root
........
.. autoclass:: mushell.thermo.cubic.rk.RedlichKwongEOSGas
  :members: relax

Linear Peneloux volume shift
----------------------------
 .. autoclass:: mushell.thermo.cubic.LinearPenelouxVolumeShift
  :members: name, category, requires
  :undoc-members: 

Non-symmetric mixing rule
-------------------------
 .. autoclass:: mushell.thermo.cubic.NonSymmmetricMixingRule
  :members: name, category, requires
  :undoc-members:  

Redich Kwong A-function
-----------------------
 .. autoclass:: mushell.thermo.cubic.rk.RedlichKwongAFunction
  :members: name, category, requires
  :undoc-members: 

Redich Kwong B-function
-----------------------
 .. autoclass:: mushell.thermo.cubic.rk.RedlichKwongBFunction
  :members: name, category, requires
  :undoc-members: 

Mathias Boston-Mathias alpha-function
-------------------------------------
 .. autoclass:: mushell.thermo.cubic.BostonMathiasAlphaFunction
  :members: name, category, requires

Appendix
--------

.. toctree::
   :maxdepth: 2
   
   cubic/alpha_extensions
