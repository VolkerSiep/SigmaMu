Quantity related functionality
==============================

.. automodule:: simu.core.utilities.quantity

Objects
-------

Quantity
^^^^^^^^
.. autoclass:: simu.core.utilities.Quantity
    :members:

.. auto

SymbolQuantity
^^^^^^^^^^^^^^
.. autoclass:: simu.core.utilities.SymbolQuantity
    :members:

QFunction
^^^^^^^^^
.. autoclass:: simu.core.utilities.QFunction
    :members:

Global Functions
----------------
The following functions redefine mathematical functions on the symbolic
quantities.

jacobian
^^^^^^^^
.. autofunction:: simu.core.utilities.jacobian

sum1
^^^^
.. autofunction:: simu.utilities.sum1

log
^^^
.. autofunction:: simu.utilities.log

exp
^^^
.. autofunction:: simu.utilities.exp

sqrt
^^^^
.. autofunction:: simu.utilities.sqrt

qpow
^^^^
.. autofunction:: simu.utilities.qpow

conditional
^^^^^^^^^^^
.. autofunction:: simu.utilities.conditional

qvertcat
^^^^^^^^
.. autofunction:: simu.utilities.qvertcat

base_unit
^^^^^^^^^
.. autofunction:: simu.utilities.base_unit

base_magnitude
^^^^^^^^^^^^^^
.. autofunction:: simu.utilities.base_magnitude

flatten_dictionary
^^^^^^^^^^^^^^^^^^
.. autofunction:: simu.utilities.flatten_dictionary

unflatten_dictionary
^^^^^^^^^^^^^^^^^^^^
.. autofunction:: simu.utilities.unflatten_dictionary

extract_units_dictionary
^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: simu.utilities.extract_units_dictionary
