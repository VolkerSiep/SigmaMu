Process simulation solver
=========================
One of the most common tasks is to simply solve a process model for its independent variables at constant thermodynamic and process parameters. The model is solved, if all residuals evaluate to absolute values lower their tolerance.

The process model has the form

.. math:: r(x, p, c) = 0

Here, :math:`x` is the complete set of independent variables, :math:`p` and :math:`c` the process model and thermodynamic parameters respectively, and :math:`r` the set of residuals. For the model to be solvable, the system must be square (:math:`\mathrm{dim}\ r = \mathrm{dim}\ x`), and the system matrix :math:`J_{rx}` must not be singular.

Iteratively, we can determine a raw Newton-Raphson update via

.. math:: \Delta x^0\ \leftarrow\ J_{rx}\cdot \Delta x^0 = -r


Additionally, the model domain :math:`\mathbb D \ni x` is limited, and its bounds described by

.. math:: b(x, p, c) > 0

These are not to be confused with inequality constraints. In particular, the case :math:`b = 0` is generally not feasible, and, as we solve an equation system, defining such bounds as inequality constraints would make the system non-square and thus not solvable.

Instead, we use :math:`b` to limit the step size of the current iteration step, assuming and hoping that the solver will eventually find a solution within the domain. We define

.. math::

    A = \left \{-\frac{b_i}{\Delta b_i} \bigg | \Delta b_i < 0 \right \}\quad\text{with}\quad
    \Delta b = J_{bx}\cdot \Delta x^0

Based on a margin parameter :math:`\gamma \approx 0.9` that describes how much we allow to approach the bounds, we determine the relaxation factor :math:`\alpha` as

.. math:: \alpha = \min(\{1, \gamma\cdot \min(A)\})

The applied variable update is then

.. math:: \Delta x = x^\mathrm{(k+1)} - x^\mathrm{(k)} = \alpha\cdot \Delta x^0

**Example**: A pressure is 2 bar, and the raw update vector suggests a step of -3 bar. We receive an element in :math:`A` of value 2/3. If this is the most limiting variable, the allowed step-size considering :math:`\gamma` is 0.6. As this value is less than unity, it will be applied, and the pressure will be updated from 2 to 0.2 bar - we allow to approach the margin by 90 % of the possible step size.

.. note::

    If we would allow the pressure to become truly zero, we would face two issues: Firstly, the model might divide by pressure or utilise the logarithm of it, and zero is already outside the domain. Secondly, numerical noise could still create negative values and thus consequences, even if zero still was in the domain.


With above scheme, the model is solved iteratively until the residual values drop below their tolerance.

SimulationSolver
----------------
.. autoclass:: simu.SimulationSolver
   :show-inheritance:
   :members:

SimulationSolverIterationReport
-------------------------------
.. autoclass:: simu.core.solver.simulation.SimulationSolverIterationReport
  :exclude-members: __init__
  :members:

SimulationSolverReport
----------------------
.. autoclass:: simu.core.solver.simulation.SimulationSolverReport
  :exclude-members: __init__
  :members:

SimulationSolverCallback
------------------------
.. autodata:: simu.core.solver.simulation.SimulationSolverCallback
