.. _alpha_extensions:

Extensions of alpha-functions to super-critical regions
=======================================================

Introduction
------------
Alpha-functions in cubic equations of state are designed to describe the vapour-liquid equilibrium well. The end-point is defined by the critical point. Beyond this critical point, extrapolation is non-physical and yields large errors in prysical property calculations, such as density, enthalpy, entropy and chemical potentials.

Therefore several authors, e.g. :cite:p:`Boston1980` have suggested alternative alpha-functions to be used for reduced temperatures :math:`\tau':= T / T_{c}` above unity. To be suitable, the following properties must hold for the extension:

  1. It must evaluate to unity for :math:`\tau' = 1`.
  2. It must be monotonically decreasing and approach zero for large :math:`\tau`.
  3. It must be differentiable at :math:`\tau' = 1`. Otherwise, entropic and calorimetric properties will be calculated discontinuously at the critical temperature.
  4. It should be differentiable twice  at :math:`\tau' = 1`. Otherwise, calorimetric second order properties, such as heat capacity will be calculated discontinuously at the critical temperature, and the performance of numerical solvers can suffer.
  
Defining the interface to the sub-critical expression
-----------------------------------------------------
To simplify the problem further, we can transform the variables without compromising above conditions. Firstly, we consider the square of the alpha-function: :math:`\alpha' := \alpha^2`, and secondly, we use the square root of the reduced temperature as independent variable: :math:`\tau = \sqrt{\tau'}`.

Next, we calculate the properties at :math:`\tau = 1` for the sub-critical alpha-function - which of course is individual. For the Mathias alpha-function :cite:p:`Mathias1983`, we have

.. math::

    \alpha' = 1 + m(1-\tau) - \eta\,(1-\tau)\,(0.7 - \tau^2)
    
We can easily validate that :math:`\alpha' = 1` for :math:`\tau = 1`. Further

.. math::

    \frac{\mathrm{d}\alpha'}{\mathrm{d}\tau} 
      &= -m + \eta\,(0.7 - \tau^2) + 2\,\eta\,\tau\,(1-\tau) = -m + \eta\,(-3\tau^2 + 2\tau  + 0.7)\\
    \frac{\mathrm{d}^2\alpha'}{\mathrm{d}\tau^2} 
      &= \eta\,(-6\tau + 2)
      
Evaluated at the critical temperature, we have

.. math::

    \alpha'_\tau := \left .\frac{\mathrm{d}\alpha'}{\mathrm{d}\tau}\right |_{\tau=1} = -m - 0.3\eta
    \quad\text{and}\quad
    \alpha'_{\tau\tau} := \left .\frac{\mathrm{d}^2\alpha'}{\mathrm{d}\tau^2}\right |_{\tau=1} = -4\eta

As said above, the actual expressions for :math:`\alpha'_\tau` and :math:`\alpha'_{\tau\tau}` are specific for the used expression.

Design of the super-critical expression
---------------------------------------
The Boston-Mathias extrapolation :cite:p:`Boston1980` proposes the following shape of the super-critical term (due to our definition of :math:`\tau`, the parameters :math:`c` and :math:`d` are not the same as published, but transformed for later convenience by simpler terms):

.. math::

	\alpha' = \exp \left [\frac{c}{d}(1-\tau^{d})\right ]

We can already see that :math:`\alpha'(\tau = 1) = 1` independent of the values of :math:`c` and :math:`d`. The derivatives are:

.. math::

    \frac{\mathrm{d}\alpha'}{\mathrm{d}\tau} 
      &= -c\,\tau^{d-1}\alpha'\\
    \frac{\mathrm{d}^2\alpha'}{\mathrm{d}\tau^2} 
      &= c\,\left [c\,\tau^{2d-2} + (1 - d)\tau^{d-2} \right ] \alpha'

Above is valid for :math:`d\not\in \{0, 1\}`. Evaluated at the critical temperature, we have

.. math::

    \left .\frac{\mathrm{d}\alpha'}{\mathrm{d}\tau}\right |_{\tau=1} = -c
    \quad\text{and}\quad
    \left .\frac{\mathrm{d}^2\alpha'}{\mathrm{d}\tau^2}\right |_{\tau=1} = c\,(c + 1 - d)

From here, we can solve for the parameters:

.. math::

    c = - \alpha'_\tau \quad\text{and}\quad
    d = 1 + \frac{\alpha'_{\tau\tau}}{\alpha'_\tau} - \alpha'_\tau

For the example above:

.. math::

	c = m + 0.3\eta \quad\text{and}\quad
	d = 1 + \frac{4\,\eta}{m + 0.3\eta} + m + 0.3\eta
	
Results
-------
The figure below shows the (square root of the) alpha function and its first derivative as function of reduced temperature, applied to ammonia.

.. image:: boston_mathias_siepmann_extrapolation.*
   :width: 80%

Aleady the original function :cite:p:`Boston1980` assumes the correct function value (unity) and first derivative at :math:`\tau = 1`, but the second derivative is not smooth originally, but zero at :math:`\tau = 1`. The extrapolation by :cite:p:`Boston1980` was developed for the standard SRK model (:math:`\eta = 0`), for which the curvature is zero at :math:`\tau = 1`. The authors then correctly set the derivative to zero based on this fact. No further intension in doing so has been documented.

Three years later, when :cite:p:`Mathias1983` publishes the extended alpha-function for polar substances, the correction in the extrapolation was only done to refit the first derivative. The curvature however, now having the value :math:`-4\eta`, was ignored :cite:p:`AspenTech2001`.

.. todo::

    When having a complete model running (for NH3), check how a sub-critical VLE fit extrapolates to the super-critical region.
    Compare! I could just calculate :math:`p(T, V)` for a given volume, which is determined without further tuning parameters -
    except :math:`c`. I can fit :math:`c` to sub-critical data.

 
