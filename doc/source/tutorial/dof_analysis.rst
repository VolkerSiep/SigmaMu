==========================
Degree of freedom analysis
==========================

This example of the steam generation and downstream turbine demonstrated a decent amount of couplings that contradict the simple principle that each unit operation closes the degrees of freedom (DOF) in the model to calculate its output streams and process properties based on its inlet streams and process parameters. That principle often applies in *process flowsheet* modelling, while *custom* modelling, in particular beneath the level of unit operations, sometimes creates the challenge of keeping track of the degrees of freedom.

Moreover, it is crucial to pose the model well-defined. Even with equal number of variables and constraints, models can be singular, as linear dependent residuals have accidentally been included.

Simple DOF analysis
===================
The simple way of analysing the DOF-balance of smaller models (say <100 variables) is summarized as follows:

  - Use a spreadsheet program and list all independent variables in one column. These are typically temperature, pressure, and partial molar flows for each material.
  - Leave the neighbouring column free, but copy the names of all residuals into a third column.
  - Now assign the residuals to state variables by moving them into the corresponding line in the second column.
    The rule here is that the residual must be dependent on the chosen variable - ideally dominantly dependent.

If a redundant residual occurs, one will not find a state variable to assign this or a closely related residual to, and also a state variable remains unassigned. It is in close relation to this state variable that a residual is missing.

The best is to limit this method to smaller model parts in the hierarchical context. Once established, these model parts can then be cleared and included in a simplified form into the parent context.

Advanced DOF analysis
=====================
The graphical analysis is best performed graphically.

.. todo::

    Talk about materials being defined and connected in a context. Connections can take "ownership" or not.
    Based on this, the unit then has a net number of DOF (positive or negative). It should document these
    DOFs well.
