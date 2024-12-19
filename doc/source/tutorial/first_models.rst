Creating the first models
=========================

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
   :lines: 28-31
   :lineno-start: 28
   :linenos:

This parameter store can be shared among multiple -- normally all -- material definitions, and thus holds a global set of thermodynamic parameters. Multiple stores are only required if two materials containing the same chemical species need to receive distinct values for the same thermodynamic parameter.

So far, we did not provide the concrete parameters, but the material definition has already told the store which parameters are required.
Here we retrieve the super-set of the parameter names and units of all missing parameters required::

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

This time, there are no more missing symbols, and the print statement yields an empty dictionary.

In real applications, storing the meta-data and parameters in ``yml`` files is not the most stupid idea, but you might connect to any other file format or database of your choice, as long as the source can provide the nested dictionary of properties as requested by the store.

.. note::

    As multiple parameter sources can be stacked in one store, we recommend to assign one source per bibliographic source of parameters. The models can then easily be queried for the names of the used sources and by that keep these sources traceable.


Using a material in a model
---------------------------
This is the big moment, as we now can use the material definition in an actual process model.