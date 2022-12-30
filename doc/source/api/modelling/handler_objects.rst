Handler objects
===============
These handler classes are instantiated by the
:class:`~simu.model.base.ModelInstance` class and handle the data flow
to combine all user models into a consistent overall model. A model developer
normally doesn't need to relate to these objects.

ParameterHandler
----------------
.. autoclass:: simu.model.parameter.ParameterHandler
    :members:

PropertyHandler
---------------
.. autoclass:: simu.model.property.PropertyHandler
    :members:

MaterialHandler
---------------
.. todo::

    invent

ResidualHandler
---------------
.. todo::

    invent

HierarchyHandler
----------------
.. autoclass:: simu.model.hierarchy.HierarchyHandler
    :members:

NumericHandler
--------------
.. autoclass:: simu.model.numeric.NumericHandler
    :members:
