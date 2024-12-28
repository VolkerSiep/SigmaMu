Quantity related functionality
==============================

.. automodule:: simu.core.utilities.quantity

Objects
-------

Quantity
^^^^^^^^
.. autoclass:: simu.Quantity
    :members:

.. auto

SymbolQuantity
^^^^^^^^^^^^^^
.. autoclass:: simu.SymbolQuantity
    :members:

QFunction
^^^^^^^^^
.. autoclass:: simu.QFunction
    :members:

Symbolic Functions
------------------
The following functions redefine mathematical functions on the symbolic
quantities.

jacobian
^^^^^^^^
.. autofunction:: simu.jacobian

qsum
^^^^
.. autofunction:: simu.qsum

log
^^^
.. autofunction:: simu.log

sqrt
^^^^
.. autofunction:: simu.sqrt

qpow
^^^^
.. autofunction:: simu.qpow

conditional
^^^^^^^^^^^
.. autofunction:: simu.conditional

Utility functions
-----------------

qvertcat
^^^^^^^^
.. autofunction:: simu.qvertcat

base_unit
^^^^^^^^^
.. autofunction:: simu.base_unit

base_magnitude
^^^^^^^^^^^^^^
.. autofunction:: simu.base_magnitude

flatten_dictionary
^^^^^^^^^^^^^^^^^^
.. autofunction:: simu.flatten_dictionary

unflatten_dictionary
^^^^^^^^^^^^^^^^^^^^
.. autofunction:: simu.unflatten_dictionary

extract_units_dictionary
^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: simu.extract_units_dictionary

simplify_quantity
^^^^^^^^^^^^^^^^^
.. autofunction:: simu.simplify_quantity

parse_quantities_in_struct
^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: simu.parse_quantities_in_struct

quantity_dict_to_strings
^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: simu.quantity_dict_to_strings

extract_sub_structure
^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: simu.core.utilities.extract_sub_structure
