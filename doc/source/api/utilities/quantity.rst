Quantity related functionality
==============================

.. automodule:: simu.core.utilities.quantity

Objects
-------
.. currentmodule:: simu

Quantity
^^^^^^^^
.. autoclass:: Quantity
    :members:

.. auto

SymbolQuantity
^^^^^^^^^^^^^^
.. autoclass:: SymbolQuantity
    :members:

QFunction
^^^^^^^^^
.. autoclass:: QFunction
    :members:

Symbolic Functions
------------------
The following functions redefine mathematical functions on the symbolic
quantities.

jacobian
^^^^^^^^
.. autofunction:: jacobian

qsum
^^^^
.. autofunction:: qsum

log
^^^
.. autofunction:: log

sqrt
^^^^
.. autofunction:: sqrt

qpow
^^^^
.. autofunction:: qpow

conditional
^^^^^^^^^^^
.. autofunction:: conditional

Utility functions
-----------------

qvertcat
^^^^^^^^
.. autofunction:: qvertcat

base_unit
^^^^^^^^^
.. autofunction:: base_unit

base_magnitude
^^^^^^^^^^^^^^
.. autofunction:: base_magnitude

flatten_dictionary
^^^^^^^^^^^^^^^^^^
.. autofunction:: flatten_dictionary

unflatten_dictionary
^^^^^^^^^^^^^^^^^^^^
.. autofunction:: unflatten_dictionary

extract_units_dictionary
^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: extract_units_dictionary

simplify_quantity
^^^^^^^^^^^^^^^^^
.. autofunction:: simplify_quantity

parse_quantities_in_struct
^^^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: parse_quantities_in_struct

quantity_dict_to_strings
^^^^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: quantity_dict_to_strings

extract_sub_structure
^^^^^^^^^^^^^^^^^^^^^
.. autofunction:: extract_sub_structure
