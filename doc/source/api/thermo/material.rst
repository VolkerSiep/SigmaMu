.. currentmodule:: simu

Material classes
================

This section describes the classes required to generate :class:`Material` objects. This is done by a factory class called :class:`MaterialDefinition`. The latter accumulates three pieces of information:

- The :class:`ThermoFrame` object, representing the thermodynamic model itself,
- the :class:`InitialState` object, defining the default initialisation of the object, and
- the :class:`ThermoParameterStore` object, being the source for thermodynamic parameters.

Once the :class:`MaterialDefinition` object is created, it can fabricate :class:`Material` objects according to its definition.

A :class:`ThermoParameterStore` instance administers any number of :class:`AbstractThermoSource` implementations, the actual sources of data that can be represented by parsed data files or adapters to databases. We provide some basic implementations by :class:`NestedDictThermoSource` and :class:`StringThermoSource`.

ThermoParameterStore
--------------------
.. autoclass:: ThermoParameterStore
   :members:

AbstractThermoSource
--------------------
.. autoclass:: AbstractThermoSource
   :special-members: __getitem__

NestedDictThermoSource
----------------------
.. autoclass:: NestedDictThermoSource
   :show-inheritance:

StringDictThermoSource
----------------------
.. autoclass:: StringDictThermoSource
   :show-inheritance:

MaterialDefinition
------------------
.. autoclass:: MaterialDefinition
   :members:

Material
--------
.. autoclass:: Material
   :show-inheritance:
   :special-members: __getitem__, __setitem__
   :members:

MaterialSpec
------------
.. autoclass:: MaterialSpec
   :special-members: __init__
   :members:
