Frequently asked questions
==========================

Well, this heading is a bit fake, to be honest, as I start(ed) writing this before publishing the software. Hence the frequency of at least the initial questions is rather low, if you count myself out. Nevertheless, these are questions somebody might ask.

Can I use `Pint`_ quantities?
-----------------------------
**Q:** Do I need to use ``from simu import Quantity``? Can I just use ``from pint import Quantity`` instead?

**A:** In most cases you probably can, as ``simu.Quantity`` derives from ``pint.Quantity``. However, besides using a dedicated ``pint.Registry``, we also add some unit operations, such as gauge pressure units, and functionality fpr custom ``json`` encoding. Not least, we define ``SymbolQuantity`` as a quantity object that holds a `CasADi`_ symbol as its magnitude. In conclusion: just use ``simu.Quantity``, please.

Can I use other (external) numerical solvers?
---------------------------------------------
**Q:** Why does |SigmaMu| define its own numerical solvers? Can't I better use well established high performant generic solvers?

**A:** In general yes, and there may be cases where this is a good idea. However, one inherent concept of |SigmaMu| is to use the dominantly thermodynamic nature of the models. This already reflects in the choice of thermodynamic state variables as independent variables. Further, the in-built solver respects the bounds that are defined for both thermodynamic models and process models to not leave the domain. This is extremely important when iterating on liquid volumes and chemical equilibria of trace components. If you decide to use a generic solver, try to handle the case of the solver stepping outside the model's domain gracefully.

As an example, a generic solver may step into a non-physical volume root for a given pressure. The model might converge, but give wrong results in terms of physical properties, in particular liquid density.

What about phase stability analysis?
------------------------------------
**Q:** In other simulation tools, flash calculations automatically determine the phase-region and calculate one-phase or two-phase stream properties depending on this stability analysis. In |SigmaMu|, I need to specify whether my fluid is liquid, solid or gaseous, even for equation of states. Why cannot it be more flexible?

**A:** For a sequential simulation tool, it is rather straight forward to perform phase stability analyses in a unit operation calculation, and simply forward the stream information to the downstream entities. In the equation-oriented approach, this is less straight forward, as the (dis-)appearance of phases causes a change in system size, and as such a change in degrees of freedom. As a second part of this answer, it can be stated that phase stability is seldom a practical problem when modelling steady-state systems. For a process to work, the number of phases must normally be as designed / desired. Now, there are exceptions, such as a discretised heat exchanger with on-setting condensation or evaporation. For most of such cases, there are tricks however to elegantly circumvent the need to dynamically change the number of phases in the model.

How to render degree-circles in files with property data
--------------------------------------------------------
**Q:** When I write ``yaml`` files with property data, the ``°C`` is escaped to ``\xB0C``. This is annoying, as I like to expose the files to my users and myself.

**A:** This is really just a detail on how to store the data using the `PyYaml`_ library. Use

.. code-block::

    dump(your_data, your_file, allow_unicode=True)

Reading of the same data can be done without setting the flag. However, if you use ``json``, be aware that ``json`` is designed for non-human information exchange, and the general recommendation is not to temper with it. Use ``yaml`` for files that users might need to work with.

When should I use inheritance over model hierarchy?
---------------------------------------------------
**Q:** As models are python classes, one can sometimes chose whether to instantiate another model as a child model within an outer model or to create a sub-class of that model instead. When shall I do what?

**A:** The short answer is to not use inheritance if there is no special reason for it, and even less, if the ``interface`` method had to be overwritten -- bad sign! If the definition of a specific model is in parts overly dependent on its configuration, inheritance can be considered, though the parametrisation might also just be used to include different child modules to bring in the individual aspects.

As an example, given a heat exchanger with a hot and cold fluid inlet and outlet, the actual heat transfer model could certainly be implemented by (generalisation) inheritance, but probably equally well by plugging in different sub-models.

Inheritance might rather be used if |SigmaMu| models are to be integrated into another framework, as for instance as unit operations into a graphical application. In this case, for instance, the configuration of material ports might be augmented with attributes that link in the geometric location of those ports.


