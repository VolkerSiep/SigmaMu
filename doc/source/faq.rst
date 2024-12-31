Frequently asked questions
==========================

Well, this heading is a bit fake, to be honest, as I start(ed) writing this before publishing the software. Hence the frequency of at least the initial questions is rather low, if you count myself out. Nevertheless, these are questions somebody might ask.

Can I use `Pint`_ quantities?
-----------------------------
**Q:** Do I need to use ``from simu import Quantity``? Can I just use ``from pint import Quantity`` instead?

**A:** In most cases you probably can, as ``simu.Quantity`` derives from ``pint.Quantity``. However, besides using a dedicated ``pint.Registry``, we also add some unit operations, such as gauge pressure units, and functionality fpr custom ``json`` encoding. Not least, we define ``SymbolQuantity`` as a quantity object that holds a `CasADi`_ symbol as its magnitude. In conclusion: just use ``simu.Quantity``, please.

Can I use other (external) numerical solvers?
---------------------------------------------
**Q:** Why does ``SiMu`` define its own numerical solvers? Can't I better use well established high performant generic solvers?

**A:** In general yes, and there may be cases where this is a good idea. However, one inherent concept of ``SiMu`` is to use the dominantly thermodynamic nature of the models. This already reflects in the choice of thermodynamic state variables as independent variables. Further, the in-built solver asks the (thermodynamic) model(s) in each iteration for the allowed step-size to not leave the domain of the model. This is extremely important when iterating on liquid volumes and chemical equilibria of trace components. If you decide to use a generic solver, try to handle the case of the solver stepping outside the model's domain gracefully.

What about phase stability analysis?
------------------------------------
**Q:** In other simulation tools, flash calculations automatically determine the phase-region and calculate one-phase or two-phase stream properties depending on this stability analysis. In `SiMu`, I need to specify whether my fluid is liquid, solid or gaseous, even for equation of states. Why cannot it be more flexible?

**A:** For a sequential simulation tool, it is rather straight forward to perform phase stability analyses in a unit operation calculation, and simply forward the stream information to the downstream entities. In the equation-oriented approach, this is less straight forward, as the (dis-)appearance of phases causes a change in system size, and as such a change in degrees of freedom. As a second part of this answer, it can be stated that phase stability is seldom a practical problem when modelling steady-state systems. For a process to work, the number of phases must normally be as designed / desired. Now, there are exceptions, such as a discretised heat exchanger with on-setting condensation or evaporation. For most of such cases, there are tricks however to elegantly circumvent the need to dynamically change the number of phases in the model.

How to render degree-circles in files with property data
--------------------------------------------------------
**Q:** When I write ``yaml`` files with property data, the ``°C`` is escaped to ``\xB0C``. This is annoying, as I like to expose the files to my users and myself.

**A:** This is really just a detail on how to store the data using the `PyYaml`_ library. Use

.. code-block::

    dump(your_data, your_file, allow_unicode=True)

Reading of the same data can be done without setting the flag. However, if you use ``json``, be aware that ``json`` is designed for non-human information exchange, and the general recommendation is not to temper with it. Use ``yaml`` for files that users might need to work with.