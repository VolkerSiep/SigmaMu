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

