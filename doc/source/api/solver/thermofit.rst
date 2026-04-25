===============================
Fit of thermodynamic parameters
===============================

Terms & Conditions
==================
In a thermodynamic parameter fit, the following terms are introduced:

**Thermophysical parameters**
  Model parameters that are independent of geometry and equipment, but describe the material properties. Typically, these include calorimetric, volumetric and entropic parameters in its core, transport properties, and reaction kinetic correlations.

   Examples are binary interaction coefficients, volume translation coefficients, diffusion coefficients, and kinetic reaction constants.

**Dataset**
  A dataset is a table of data, of which each row describes a reproducible *thermodynamic* observation. Here, *thermodynamic* denotes a relationship whose mathematical description includes *thermophysical parameters*. Within a dataset, each row provides the same data, hence the table is full and rectangular. Each column has a header that represents the name of the sampled variable, and a unit of measurement that is applied to all values in that column.

  An example dataset is one describing phase equilibrium for a binary mixture at constant temperature, and it can have the following form:

  ===== ===== =====
  p     x     y
  ===== ===== =====
  *bar* *--*  *--*
  3.07  0.021 0.832
  6.99  0.060 0.920
  11.05 0.110 0.946
  14.97 0.183 0.957
  ...
  ===== ===== =====

  - Instead of creating many small data sets, it can be a valid choice to augment them, for instance by adding a temperature column and including multiple isotherms.
  - The model might support a weight parameter to define the impact of each row by normalization of the defined penalties. This can be handled like any other parameter and its data being added as a column to the dataset.

**Penalty model**
  To interpret a row in a *dataset*, its values are provided to a small model that evaluates its contribution to a penalty function. Via *contributions*, one penalty model can be applied for one or more datasets, and one dataset can be interpreted by one or more penalty models.

  A penalty model is a :class:`~simu.Model` implementation that calculates dimensionless penalty terms as function of model parameters. These parameters are defined and mapped from the linked data set.

  As an example, a penalty model can receive two-phase flash data and evaluate the distance from equilibrium as a penalty, for instance expressed as a normalized difference in chemical potentials :math:`\mu_i^\mathrm{gas} - \mu_i^\mathrm{liquid}`.

**Contribution**
  A contribution [to a thermodynamic parameter fit] links a data set to a penalty model, defines the mapping between dataset columns and model parameters, and specifies which model parameters are interpreted as penalty contributions.

**Evaluation model**
  For assessing a thermodynamic data fit, the raw penalties used by the *penalty model* are not always suitable for human interpretation. As mentioned above, a penalty model might penalize a difference in chemical potential, while a human being prefers assessing deviations in tangible properties, like pressures and compositions. Of course, in some cases, one might simply be able to reuse the penalty model as the evaluation model.

  As such, the *evaluation model* receives parameter values from a *dataset* and calculates either tangible deviations or calculated versions of further dataset columns that represent observed result.

  As for the two-phase flash example, the evaluation model can solve the flash at given temperature and liquid fraction to then report the calculated pressure and gas fraction for comparison with the tabulated data.

**Evaluation**
  Alike to a **contribution**, an evaluation maps a *dataset* to an *evaluation model*, maps input parameters and defines properties to be reported.

  There can be multiple evaluations (and models) applied to one dataset. This can be required to visualize deviations in multiple dimensions or to offer alternative interpretations of deviations.

Mathematical description
========================
Data fit
--------
The *penalty model* is represented by residuals :math:`r(x, p, \tau) = 0`, whereas :math:`x` is the internal state of the model, :math:`p` are the model parameters, and :math:`\tau` are the thermodynamic parameters to be estimated. We demand :math:`\dim r = \dim x` and a well-formed (non-singular) equation system with a solution within the domain of :math:`x`. Further, the model calculates dimensionless one or more penalty terms :math:`q(x, p)`.

For each row :math:`j` in the dataset :math:`i`, the model receives its data as parameters :math:`p_{ij}` and provides :math:`r_{ij}` as function of :math:`x_ {ij}`, further :math:`q_{ij}`, as well as the later required Jacobian matrices.
For practical convenience, we also introduce a weight factor :math:`w_i` to control the impact for entire data sets.

The data fit is the minimization

.. math::

    \min_{x_{ij},\tau} \frac12 \sum_{ij} w_i^2\,q_{ij}(x_{ij}, p_{ij}, \tau)^2
    \quad\text{s.t.}\quad
    r_{ij}(x_{ij}, p_{ij}, \tau) = 0\quad \forall_{ij}

The pairs :math:`(ij)` can be combined into :math:`k`, and once the parameters are set from the dataset, the :math:`p_{ij}` can be omitted:

.. math::

    \min_{x_k,\tau} \frac12 \sum_k w_k^2\,q_k(x_k, \tau)^2
    \quad\text{s.t.}\quad r_k(x_k, \tau) = 0\quad \forall_k

Notable, for a given :math:`\tau`, the systems :math:`r_k(x_k, \tau) = 0` are square and solvable, as they are small independent process models.
This allows to solve each of them independently (*hint: and in parallel*). At the solution point, :math:`r_k = 0`, the Jacobian is derived via the total differential:

.. math::
   :nowrap:

   \begin{align*}
   \mathrm{d}r &= \left . \frac{\partial r_k}{\partial x_k} \right |_{\tau}\,\mathrm{d}x_k +
                  \left . \frac{\partial r_k}{\partial \tau} \right |_{x_k}\,\mathrm{d}\tau = 0\\
   \Rightarrow J_{x_k, \tau} &:= \frac{\mathrm{d}x_k}{\mathrm{d}\tau} =
     -\left (\left . \frac{\partial r_k}{\partial x_k} \right |_{\tau} \right )^{-1}\cdot
     \left . \frac{\partial r_k}{\partial \tau} \right |_{x_k}
   \end{align*}

For the penalty properties :math:`q_k`, the sensitivity is

.. math::

   J_{q_k, \tau} = \left . \frac{\partial q_k}{\partial x_k} \right |_{\tau}\,J_{r_k, \tau}
      + \left . \frac{\partial q_k}{\partial \tau} \right |_{x_k}

The stationary condition then becomes

.. math:: \sum_k w_k^2\,q_k\,J_{q_k, \tau} = 0

Neglecting higher order derivatives (*normally not a problem here*), the linearized least squares residual formulation is

.. math::

   \left [\sum_k w_k^2\,q_k\,J_{q_k, \tau} \right ]+
     \left [ \sum_k w_k^2 \left ( J_{q_k, \tau} \right )^\mathrm{T}\,J_{q_k, \tau} \right ]\,\Delta\tau \approx 0

This system is used as a second order update scheme for :math:`\tau`.



Evaluation
----------
