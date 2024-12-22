Solving models
==============

While process models are pretty by themselves, their real purpose is to be numerically solved and their results to be analyzed.
In the previous section, we created a model with a pure methane flow and specified temperature, pressure and volume flow, as for each material with :math:`n` species, `math:`n+2` residuals are to be defined for the model to be *square*.

Process simulation
------------------
For ordinary process simulation, we need our models to be square, that is, balancing the number of residuals and independent variables.

