.. currentmodule:: simu.core

.. _pytest: https://docs.pytest.org/
.. _CasADi: https://web.casadi.org
.. _NumPy: https://numpy.org/
.. _SciPy: https://scipy.org/
.. _Pint: https://pint.readthedocs.io
.. _PyYAML: https://pyyaml.org/
.. _Sphinx: https://www.sphinx-doc.org
.. _sphinxcontrib-bibtex: https://github.com/mcmtroffaes/sphinxcontrib-bibtex
.. _Matplotlib: https://matplotlib.org/


Getting started
===============
Installation
------------
By intention, installation shall be as easy as

.. code-block::

    pip install SiMu

By the time this is read by anybody but me, the package is available on ``PyPi``, hence this shall work. The tests can then be run, if `pytest`_ is available (see below) via

.. code-block::

   pytest --pyargs simu

Also this is supposed to work without any failed tests.

``SiMu`` depends on the following packages:

========= =================================================================
Name      What for
========= =================================================================
`CasADi`_ All symbolic algebra in the background uses it-
`Pint`_   All modelling happens via ``pint``'s ``Quantity`` class
`NumPy`_  Number chrunching, in particular linear algebra
`SciPy`_  For advanced numeric methods, utilised for instance by the solver
`PyYAML`_ The chosen format for configuration files.
========= =================================================================

For development, we require additionally

======================= =================================================
Name                    What for
======================= =================================================
`pytest`_               For running the unit tests
`matplotlib`_           In examples, we like to plot results sometimes
`Sphinx`_               The documentation is built with it.
`sphinxcontrib-bibtex`_ Handling of bibliographics in documentation
======================= =================================================

Hello World
-----------
To call the following a *process model* is quite an insult to actual process models, but it is a start:

.. literalinclude:: examples/hello_world.py
   :language: python
   :linenos:
   :pyobject: Square

To create a model, we derive from the :class:`Model <model.base.Model>` class and need to implement the following two methods:

In the ``interface`` method, we define a model parameter called ``length`` and declare to calculate a property called ``area``.
In the ``define`` method, the property ``area`` is assigned to be the square of the ``length`` parameter.

Next, we want to do something with the model. To do that, we wrap it into an instance of :class:`NumericHandler <model.numeric.NumericHandler>` and create a :class:`QFunction <utilities.quantity.QFunction>` object:

.. code-block::

    numeric = NumericHandler(Square.top())
    func = numeric.function


Now we can print the argument structure of the function:

>>> print(func.arg_structure)
{'model_params': {'length': 'm'}, 'vectors': {'states': ''}}

The :class:`simu.core.utilities.quantity.QFunction` is a function object that works on nested dictionaries of `pint`_ quantities, whereas the magnitude is of type ``casadi.SX``. When querying the ``arg_structure``, the values of the dictionaries are the units of measurements as defined in the model's interface.

The ``length`` parameter is a model parameter and registered as such. We have not defined any materials, and hence there are no thermodynamic state variables.

Likewise, we can query what the model will return:

>>> print(func.result_structure)
{'model_props': {'area': 'm ** 2'}, 'vectors': {'residuals': ''}}

As there are no model constraints, the ``residuals`` section here is empty, but we find the calculated ``area`` as a model property.
