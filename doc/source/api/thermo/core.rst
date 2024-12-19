.. currentmodule:: simu

Generic classes
===============
A thermodynamic model in ``SiMu`` is composed of :class:`ThermoContribution` objects stacked on top of a :class:`StateDefinition` instance. The latter one interprets the numerical state vector as physical quantities, such as temperature, pressure and mole flows. Each contribution then builds on the already defined quantities to add new ones, until the model is complete.

The :class:`ThermoFactory` administers the construction of this structure by allowing to first register the contribution classes, and then construct models based on a given structure. The model is then encapsulated into a :class:`ThermoFrame` object.

ThermoContributionDict
----------------------
.. autodata:: simu.core.thermo.frame.ThermoContributionDict
    :annotation:

ThermoFactory
-------------
.. autoclass:: ThermoFactory
   :special-members: __init__
   :members:

ThermoFrame
-----------
.. autoclass:: ThermoFrame
   :members:
   :special-members: __call__

ThermoContribution
------------------
.. autoclass:: ThermoContribution
   :members:
   :private-members: _tensor_structure, _vector

StateDefinition
---------------
.. autoclass:: StateDefinition
   :members:

SpeciesDefinition
-----------------
.. autoclass:: SpeciesDefinition
   :members:

SpeciesDB
---------
.. autoclass:: SpeciesDB
   :special-members: __init__
   :members:

InitialState
------------
.. autoclass:: InitialState
   :members:
