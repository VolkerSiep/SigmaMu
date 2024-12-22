The first process model
=======================

.. _CasADi: https://web.casadi.org

Recap
-----

A brief description on how models are created is given in the :ref:`getting started hello world` paragraph of the getting started section. The example was

.. literalinclude:: ../examples/hello_world.py
   :language: python
   :linenos:

To recap, the model declares an interface telling other (parent) models that it has a parameter called ``length`` and calculates a property called ``area``. Then, the definition implements the relationship.

The attentive reader might at this point have realised that no thermodynamic model was involved. The entire example is somehow remote to classical process engineering unless one calculates the cross-section of a square duct. Let us get this rectified in this section, and bring in our ideal gas model from the previous section.

To do so, we first need to create a :class:`simu.MaterialDefinition` object, and by this follow the proper way to build up a simulation, and in practice, you might soon build up a repository of materials required for your field of application.

Creating a material definition
------------------------------

We can reuse much of the previous code, resulting into a :class:`simu.ThermoFrame` object for a pure methane ideal gas. Only now, we store configuration data in yaml files. Let us have one file with some chemical species and their formulae (``species_db.yml``):

.. literalinclude:: ../examples/species_db.yml
   :language: yaml
   :linenos:

Further, we store the thermodynamic model structures in another file (``thermo_model_structures.yml``):

.. literalinclude:: ../examples/thermo_model_structures.yml
   :language: yaml
   :linenos:

And finally, the actual thermodynamic parameters (``ideal_gas_param.yml``):

.. literalinclude:: ../examples/ideal_gas_param.yml
   :language: yaml
   :linenos:

Let's first import some classes and read in those files:

.. literalinclude:: ../examples/ideal_gas_material.py
   :language: python
   :lines: 1-19
   :linenos:

Now, as before, we can create the factory and from there the frame for our thermodynamic model:

.. literalinclude:: ../examples/ideal_gas_material.py
   :language: python
   :lines: 21-26
   :lineno-start: 21
   :linenos:

Next, we go for the material definition, requiring also the initial state and a :class:`simu.ThermoParameterStore` object:

.. literalinclude:: ../examples/ideal_gas_material.py
   :language: python
   :lines: 28-32
   :lineno-start: 28
   :linenos:

This parameter store can be shared among multiple -- normally all -- material definitions, and thus holds a global set of thermodynamic parameters. Multiple stores are only required if two materials containing the same chemical species need to receive distinct values for the same thermodynamic parameter.

So far, we did not provide the concrete parameters, but the material definition has already told the store which parameters are required. Line 32 prints the super-set of the parameter names and units of all missing parameters required::

    {'H0S0ReferenceState': {'T_ref': 'K ',
                            'dh_form': {'Methane': 'J / mol '},
                            'p_ref': 'Pa ',
                            's_0': {'Methane': 'J / K / mol '}},
     'LinearHeatCapacity': {'cp_a': {'Methane': 'J / K / mol '},
                            'cp_b': {'Methane': 'J / K ** 2 / mol '}}}

Finally, we provide the already read parameters to the store:

.. literalinclude:: ../examples/ideal_gas_material.py
   :language: python
   :lines: 33-34
   :lineno-start: 33
   :linenos:

This time, there are no more missing symbols, and the print statement prints an empty dictionary.

In real applications, storing the meta-data and parameters in ``yml`` files is not the most stupid idea, but you might connect to any other file format or database of your choice, as long as the source can provide the nested dictionary of properties as requested by the store.

.. note::

    As multiple parameter sources can be stacked in one store, we recommend to assign one source per bibliographic source of parameters. The models can then easily be queried for the names of the used sources and by that keep these sources traceable.

Using a material in a model
---------------------------
This is the big moment, as we now can use the material definition in an actual process model. The above created :class:`simu.MaterialDefinition` object can be global for the entire project along with all other material definitions that you might need.

The following model is a *hello world* example for using such material:

.. literalinclude:: ../examples/material_model.py
   :language: python
   :lines: 1-16
   :linenos:

Here we first define the three parameters ``T``, ``p`` and ``V`` that determine our system. These three parameters also constitute the interface of our model. The definition creates a methane flow :class:`simu.Material` object from our definition. Finally, lines 15-17 constrain the system to the direct specifications of the parameter variables.

Well, the above syntax is very verbose, but might get into the way with regards to coding efficiency and the ambitions to keep lines short and to the point. For this reason, we define a subclass to :class:`simu.Model`, namely :class:`simu.AModel` that does nothing but defining abbreviations. As such, we can reduce the above model to:

.. literalinclude:: ../examples/material_amodel.py
   :language: python
   :lines: 1-16
   :linenos:

This is as much as we can do without entirely drowning out the pythonic way of coding.

Either way, here we are with a complete process model. By creating a :class:`simu.NumericHandler`, we obtain a function object representing our model. The function argument and result is a nested structure of quantities:

.. literalinclude:: ../examples/material_amodel.py
   :language: python
   :lines: 21-24
   :lineno-start: 21
   :linenos:

This yields::

    {'model_params': {'T': <Quantity(25, 'degree_Celsius')>,
                      'V': <Quantity(10, 'meter ** 3 / hour')>,
                      'p': <Quantity(1, 'bar')>},
     'thermo_params': {'default': {'H0S0ReferenceState': {'T_ref': <Quantity(25, 'degree_Celsius')>,
                                                          'dh_form': {'Methane': <Quantity(-74.87, 'kilojoule / mole')>},
                                                          'p_ref': <Quantity(1, 'bar')>,
                                                          's_0': {'Methane': <Quantity(188.66, 'joule / kelvin / mole')>}},
                                   'LinearHeatCapacity': {'cp_a': {'Methane': <Quantity(35.69, 'joule / kelvin / mole')>},
                                                          'cp_b': {'Methane': <Quantity(50.0, 'millijoule / kelvin ** 2 / mole')>}}}},
     'vectors': {'states': <Quantity([400, 200000, 1], 'dimensionless')>}}

Firstly, we can recognize the model parameters, the thermodynamic parameters, and the thermodynamic state of our material. The latter is stored in a dimensionless vector for the purpose of numerical solving. Later-on, we show how this vector, and/or individual parameters can be substituted by `CasADi`_ symbols and thus become free variables in a calculation.

Further, we can query the result by calling the function with this argument:

.. literalinclude:: ../examples/material_amodel.py
   :language: python
   :lines: 27-30
   :lineno-start: 27
   :linenos:

This yields::

    {'residuals': {'T': <Quantity(-101.85, 'kelvin')>,
                   'V': <Quantity(-0.0138511475, 'meter ** 3 / second')>,
                   'p': <Quantity(-100000.0, 'pascal')>},
     'thermo_props': {'source': {'S': <Quantity(194.096662, 'watt / kelvin')>,
                                 'T': <Quantity(400.0, 'kelvin')>,
                                 'T_ref': <Quantity(298.15, 'kelvin')>,
                                 'V': <Quantity(0.01662892523630648, 'meter ** 3 / second')>,
                                 'mu': <Quantity(-148614.303, 'joule / mole')>,
                                 'mw': <Quantity(0.016042999999999998, 'kilogram / mole')>,
                                 'n': <Quantity(1.0, 'mole / second')>,
                                 'p': <Quantity(200000.0, 'pascal')>,
                                 'p_ref': <Quantity(100000.0, 'pascal')>}},
     'vectors': {'residuals': <Quantity([-1.01850000e+09 -4.98641309e+08 -1.00000000e+07], 'dimensionless')>}}


Here we see the residuals as physical quantities, but also converted to a dimensionless vector, representing the quotient of residuals and their tolerances. Thermodynamic properties are included, and model properties would, if there were any.

The volume is calculated to 59.86 m3/hr, but we specified 10 m3/hr, and also the pressure and temperature are not yet as desired. The specifications are only fulfilled once the residuals are brought down to values below their tolerances.