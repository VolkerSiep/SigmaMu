Getting started
===============

Installation
------------
By intention, installation shall be as easy as

.. code-block::

    pip install SiMu

By the time this is read by anybody but me, the package is available on ``PyPi``, hence this shall work. The tests can then be run, if ``pytest`` is available (see below) via

.. code-block::

   pytest --pyargs simu

Also this is supposed to work without any failed tests.

``SiMu`` depends on the following packages:

========== ============================================================
Name       Why?
========== ============================================================
``casadi`` All symbolic algebra in the background
``numpy``  Numerical processing, in particular linear algebra
``scipy``  The solver utilises advanced numerics, provided by ``scipy``
``pint``   All modelling happens via ``pint``'s ``Quantity`` class
``pyyaml`` The chosen format for configuration files.
========== ============================================================

For development, we require additionally

======================== =================================================
Name                     Why?
======================== =================================================
``Sphinx``               The documentation is built with it.
``sphinxcontrib-bibtex`` Handling of bibliographics in documentation
``pytest``               For running the unit tests
``matplotlib``           In examples, we like to plot results sometimes
======================== =================================================

Hello World
-----------
To call the following a *process model* is quite an insult to actual process models, but it is a start:

.. code-block::

    from simu import Model

    class Square(Model):
        """A model of a square"""

        def interface(self):
            self.parameters.define("length", 10, "m")
            self.properties.declare("area", "m^2")

        def define(self):
            self.properties["area"] = self.parameters["length"] ** 2

To create a model, we derive from the ``Model`` class and need to implement the following two methods:

In the ``interface`` method, we define a model parameter called ``length`` and declare to calculate a property called ``area``.
In the ``define`` method, the property ``area`` is assigned to be the square of the ``length`` parameter.

Well, while this is valid code, we provide a shorter syntax as follows:

.. code-block::

    from simu import Model

    class Square(Model):
        """A model of a square"""

        def interface(self):
            self.param("length", 10, "m")
            self.prop("area", "m^2")

        def define(self):
            self.y["area"] = self.p["length"] ** 2

.. todo::

    This is not implemented yet, but would be nice, wouldn't it?

Next, we want to do something with the model.
