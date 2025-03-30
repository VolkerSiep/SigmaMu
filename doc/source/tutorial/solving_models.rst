Process simulation
==================

While process models are pretty by themselves, their real purpose is to be numerically solved and their results to be analyzed.
In the previous section, we created a model with a pure methane flow and specified temperature, pressure and volume flow, as for each material with :math:`n` species, `math:`n+2` residuals are to be defined for the model to be *square*.

For ordinary process simulation, we need our models to be square, that is, balancing the number of residuals and independent variables. Furthermore, the system needs to be well posed in the sense that the residuals are not dependent to each other, and the free variables are thus determined by the residuals. Let us review our example:

.. exampleinclude:: material_model.py
   :language: python
   :lines: 1-16
   :linenos:

Here, we specify :math:`T`, :math:`p`, and :math:`V`, which is valid, because the volume spec determines the molar flow at given temperature and pressure. As a counter-example, specifying density instead of volume would fail. The density of the material is already fixed by specifying temperature and pressure, but the system size would then still be undetermined.

A trivial check: **Each material must be part of at least one process constraint in extensive variables, and there must be at least as many constraints in extensive variables as materials.** This condition is not specific to ``SigmaMu``, but a general condition for a system not to be singular.

Running the solver
------------------

.. testsetup::

    >>> from simu.examples.material_model import Source

To solve the example module, we first wrap the model into a numeric handler:

>>> from simu import NumericHandler, SimulationSolver
>>> numeric = NumericHandler(Source.top())

Next steps are to instantiate a :class:`~simu.SimulationSolver` on this numeric handler, and then to solve the system:

>>> solver = SimulationSolver(numeric)
>>> result = solver.solve()
Iter   LMET   Alpha   Time                           Limit on bound                             Max residual
----- ----- ------- ------ ---------------------------------------- ----------------------------------------
    0   9.0    0.83   ...                source/IdealMix/n/Methane                                        T
    1   8.2       1   ...                                                                                 T
    2   6.4       1   ...                                                                                 V
    3  -7.3       1   ...                                                                                 V

The second column ``LMET`` is the logarithmic maximum error to tolerance ratio. An initial value of 7 to 9 is typical, for instance if a temperature is 1 to 100 K away from its solution, while the tolerance is 1e-7 K. When the Newton-Raphson method grips, and full steps can be taken, the decrease of LMET ideally doubles in each iteration. Once being less than zero, all residual values are below their tolerances, and the model is solved. Linear or nearly linear models show above behaviour: The LMET jumps straight down to the level of numerical precision (here roughly 1e-14 of the value magnitudes).

The relaxation factor is called ``Alpha``. The first iteration attempts to reduce the molar flow in one step dangerously close to the domain limit (:math:`n > 0`). The solver hence stabilises the step by only admitting 83 %. The residual of maximal value to tolerance ratio is shown in the right-most column.

The solving process returns a :class:`~simu.core.solver.simulation.SimulationSolverReport` object, which contains the content of above printed table for further analysis, but also the model's state vector in the solution point and a property to evaluate the model's properties. In this example, the only interesting part are the thermodynamic properties of the source stream:

>>> from pprint import pprint
>>> pprint(result.properties["thermo_props"]["source"])
{'S': <Quantity(21.140163, 'watt / kelvin')>,
 'T': <Quantity(298.15, 'kelvin')>,
 'T_ref': <Quantity(298.15, 'kelvin')>,
 'V': <Quantity(0.002777777777777777, 'meter ** 3 / second')>,
 'mu': {'Methane': <Quantity(-131118.979, 'joule / mole')>},
 'n': {'Methane': <Quantity(0.112054293180843, 'mole / second')>},
 'p': <Quantity(100000.0, 'pascal')>,
 'p_ref': <Quantity(100000.0, 'pascal')>}

The least boring result of this simulation is the calculated molar flow of methane:

>>> print(result.properties["thermo_props"]["source"]["n"]["Methane"].to("kmol/day"))
9.68149... kilomole / day

Changing parameters
-------------------
The process model already defines values for each parameter. These values are however only meant to be default values, suitable for testing and to self-document what kind of values will make sense as input to the model. These values can be changed when a model becomes a sub-model in a hierarchical setting, or in the :class:`~simu.NumericHandler` interface.

For extra convenience, the solver object provides direct mutable access via :meth:`simu.SimulationSolver.model_parameters` to the model's parameters:

>>> pprint(solver.model_parameters)
    {'model_params': {'T': <Quantity(25, 'degree_Celsius')>,
                      'V': <Quantity(10, 'meter ** 3 / hour')>,
                      'p': <Quantity(1, 'bar')>},
     'thermo_params': {'default': {'H0S0ReferenceState': {'T_ref': <Quantity(25, 'degree_Celsius')>,
                                                          'dh_form': {'Methane': <Quantity(-74.87, 'kilojoule / mole')>},
                                                          'p_ref': <Quantity(1, 'bar')>,
                                                          's_0': {'Methane': <Quantity(188.66, 'joule / kelvin / mole')>}},
                                   'LinearHeatCapacity': {'cp_a': {'Methane': <Quantity(35.69, 'joule / kelvin / mole')>},
                                                          'cp_b': {'Methane': <Quantity(50.0, 'millijoule / kelvin ** 2 / mole')>}}}},
     'vectors': {}}

Let's modify the input:

>>> from simu import Quantity
>>> solver.model_parameters["model_params"]["T"] = Quantity(120, "degC")
>>> result = solver.solve(output="none")
>>> print(result.properties["thermo_props"]["source"]["n"]["Methane"].to("kmol/day"))
7.3420... kilomole / day

As a result of increasing the temperature, the molar flow at constant volume flow becomes less. We could also change thermodynamic parameters at this point -- only that due to the applied ideal gas law, none of them has impact on the calculated molar flow.

In above call, we also omitted the output by setting the ``output`` stream to ``None``.

Using the callback function
---------------------------
Sometimes, for instance for debugging, it is useful to assess the model's state during the solving process in each iteration, and possibly even decide to stop the iterations based on custom conditions. The :class:`~simu.SimulationSolver` object offers to install a callback function:

>>> def my_callback(iteration, iter_report, state, prop_func):
...     props = prop_func(state)
...     print(iteration, props["thermo_props"]["source"]["n"]["Methane"].to("kmol/day"))
...     return True

>>> solver.set_option("call_back_iter", my_callback)
>>> solver.model_parameters["model_params"]["T"] = Quantity(-20, "degC")
>>> result = solver.solve()
    0 9.9565... kilomole / day
    1 11.402... kilomole / day

Here we observe the calculated molar flow for each iteration. The callback function returns ``True`` to proceed with the iterations until convergence is obtained.

Handling starting values
------------------------
The model hosts its initial state, defined through the :class:`~simu.MaterialDefinition` objects. For our freshly instantiated example model, this is

>>> numeric = NumericHandler(Source.top())
>>> pprint(numeric.export_state())
{'non-canonical': {},
 'thermo': {'source': {'T': '400 K',
                       'n': {'Methane': '1 mol'},
                       'p': '200000 Pa'}}}

Once we solve the model, the solver will by default retain the solution state (see ``retain_solution``):

>>> solver = SimulationSolver(numeric)
>>> result = solver.solve(output="none")
>>> print(f"This took {len(result.iterations)} iteration(s).")
This took 4 iteration(s).

>>> pprint(numeric.export_state())
{'non-canonical': {},
 'thermo': {'source': {'T': '298.1... K',
                       'n': {'Methane': '0.11205... mol'},
                       'p': '100000 Pa'}}}

If we run the solver again without changing any input, we get:
>>> result = solver.solve(output="none")
>>> print(f"This took {len(result.iterations)} iteration(s).")
This took 1 iteration(s).

We can import the state back into the model, for instance if we would like to start from a previous solution. Typically, one would store the initial state structure in a json file. Here, for simplicity:

>>> numeric = NumericHandler(Source.top())  # Let's create a fresh model
>>> start = {
...     'non-canonical': {},
...     'thermo': {'source': {'T': '25 degC',
...                           'n': {'Methane': '0.11205429318084301 mol'},
...                           'p': '1 bar'}}
... }
>>> _ = numeric.import_state(start)

Now we only need 1 iteration to solve the model, as we picked up the prior solution as start values:

>>> solver = SimulationSolver(numeric)
>>> result = solver.solve(output="none")
>>> print(f"This took {len(result.iterations)} iteration(s).")
This took 1 iteration(s).

In fact, the performed iteration was not a full evaluation, but only a positive check for convergence.

>>> print(result.iterations[0].lmet)
-7.3...