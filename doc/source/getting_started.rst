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

The model has also collected the default input values for us.

>>> args = numeric.arguments
>>> print(args)
{'vectors': {'states': <Quantity(0x1, 'dimensionless')>},
 'model_params': {'length': <Quantity(10, 'meter')>},
 'thermo_params': {}}

For a larger real-life problem, this would also include the initial set of independent variables (``state``) and all thermodynamic parameters, collected from the various data sources. Here we see only the ``length`` parameter as being 10 m.

We can overwrite that parameter by changing its value in the obtained structure

>>> args[NumericHandler.MODEL_PARAMS]["length"] = Quantity(20, "cm")
>>> result = func(args)
>>> print(f"{result[NumericHandler.MODEL_PROPS]['area']:.3fP~}")
0.040 m²

Note that the formatting of the physical quantities is utilising `Pint`_ functionality.


Normal project structure
------------------------
Above example is cute but no usecase for using ``SiMu``. Admittingly, the area of a square can instead be calculated in one line of code. A real project starts with the setup of **thermodynamic models** to calculate physical properties of the materials that are part of the model. In ``SiMu``, a thermodynamic model is represented by a :class:`ThermoFrame <thermo.frame.ThermoFrame>` object and consists of :class:`ThermoContribution <thermo.frame.contribution.ThermoContribution>` object, the latter of which can be combined and extended with high flexibility. Examples for such contributions are

  - Ideal gas heat capacity
  - Ideal mix
  - Mixing rules in various settings, for instance equations of states
  - :math:`\alpha`-functions in cubic equations of state
  - Poynting corrections in Gibbs excess models

Once the thermodynamic model structures are defined, data sources are organized to provide the **thermodynamic parameters**, such as standard state parameters or critical constants. Naturally, the set of required parameters depend on the model structure and the set of chemical species.

Now the thermodynamic models are in place, and we can define materials. Material definitions, on top of the thermodynamic model singeltons, define the utilised set of chemical species and a representative initial state. Material definitions are then used within the :class:`Model <model.base.Model>` class to define instances, defining flows of materials or stagnant states, such as phase interface conditions.

Each ``Material`` instance introduces independent variables (``states``) to the model, which uses the thermodynamic properties in combination with model parameters (mostly operational and design parameters) to evaluate both model properties of interest, but also ``Residual`` properties. In standard simulations, those residuals are bought to zero by solving over the state variables.

However, instead of plain solving, one can conduct

  - **Data reconciliation**: Minimising the deviation between measured data and calculated model properties over the state of the model;
  - **Parameter fit**: Fitting model parameters common over multiple data-sets in combination with individual model states, constrained by the model's residuals, to minimise deviation between measured data and calculated model properties;
  - **Thermodynamic parameter fit**: Like *Parameter fit*, but specifically with thermodynamic parameters, often utilising laboratory data concerning equilibrium or calorimetric data.
  - **Parameter optimisation**: Minimizing an objective function as function of model properties over the thermodynamic state and model parameters, constrained by the model's residuals.

The derivatives required to efficiently perform these disciplines can easily be obtained, based on `CasADi`_ functionality.

.. todo::

    Draw a diagram and link all the classes here.
