# What is *SigmaMu*?

***SigmaMu* is a python library for first principle steady-state modelling of chemical processes.**

<img src="doc/source/figures/simu_logo.jpeg" alt="SigmaMu logo" style="display: block; margin-left: auto; margin-right: auto; width: 400pt">

![Logo](https://raw.githubusercontent.com/VolkerSiep/SigmaMu/doc/source/figures/simu_logo.jpeg)

Installation is easily done via `pip` from `pypi` via

```
   pip install SigmaMu
```

**Once public, the detailed documentation will be on [readthedocs](https://docs.readthedocs.io) and linked from here. Also, the package will be registered in [PyPi](https://pypi.org/).**

With focus on **rigorous thermodynamic models**, the first layer provides functionality to flexibly define, combine thermodynamic and parameterize model contributions into accurate descriptions of the processed materials. One can use pre-defined model structures, exchange single contributions, such as the &alpha;-function in equations of states, or implement entirely new contributions and model structures.

The second layer provides a framework for **steady-state process modelling** in a fully hierarchical context. Concepts such as material and energy balances can be defined and reused without code duplication. Models define a self-documenting interface of interacting materials, parameters and properties. 

As the toplevel model is used to generate a numeric handler, it is fed to the third layer of functionality, represented by the **numerical solvers**. These are tailor-made for process models built on thermodynamic models, characterized by restricted domains and solutions close to the domain limits, badly scaled variables, and strong non-linearity. 

**SigmaMu** is heavily based on [CasADi](https://web.casadi.org) for the mathematical representation and differentiation of the models, and [Pint](https://pint.readthedocs.io) for flexible and consistent handling of physical dimensions and units of measurements.

Being a python library built on well establish dependencies only, **SigmaMu** is platform independent and moderate regarding its own footprint. In summary, the focus points are

- Scalable process modelling environment, suitable for modelling industrial production processes with high accuracy and in detail, utilising available CPUs.
- Covering standard process simulation, but also for instance optimization, data reconciliation, data fit of process and/or thermodynamic parameters.
- Extensible environment, enabling encapsulation and protection of intellectual property on repository level
  - in-house development and efficient reuse of thermodynamic models and process model parts.
  - The *LGPL* license allows use of **SigmaMu** even in commercial settings and without forcing *GPL*-type licenses onto the derived software, as long as **SigmaMu** itself (the core) is kept as is.
- Compatibility to other operative systems, for instance to run models in cloud-deployed Linux docker containers.
- Enabling model development with elements of modern software development, including 
  - concurrent development by utilising version control systems
  - standard ways of code/model testing using standard frameworks, such as [pytest](https://docs.pytest.org/)
  - in-code documentation of models using [sphinx](https://www.sphinx-doc.org) with for instance `autodoc` and `mathjax`.
- Maximising model reliability and reducing modelling errors through
  - hierarchical modelling, enabling efficient independent testing of minor sub-model parts
  - forcing consistent modelling with respect to physical dimensions by using [Pint](https://pint.readthedocs.io) for all exchange of physical quantities.
- Extensive documentation, use of `doctest` and *unit tests*, and examples to make the software most accessible for new users.  

## ... and what not?
- **SigmaMu** makes no attempt to provide a **graphical user interface** to make its functionality easier accessible and increase usability. This can surely be done and might be a good idea for some use-cases, but is not the primary goal of this project.
- **SigmaMu** is not a general purpose modelling tool, such as for instance [Modelica](https://modelica.org/). While it is to some degree possible to model non-chemical systems, one is surely better off by using such generic tools.
- While the modelling approach of **SigmaMu** strongly promotes better understanding of thermodynamics and chemical process modelling, it is not removing the necessity to define the model on equation level. However, some models are readily available for reuse, and experts might develop further in-house model parts that need not be understood by downstream developers.
- Physical property data is less *public domain* than you might think, and while it is perfectly legal to use numerous online databases and textbooks as free sources in academic and cooperate work, extracting compiled data and republishing them open source is not generally admitted. Therefore, ``SigmaMu`` cannot contain significant physical property data by itself.
 
  The idea is that physical data for the particular application is gathered legally in derived projects under much less strict conditions. However, we might seek permission and be granted to provide data for a subset of common species, not least for demonstration purposes. 

# Example
Below code imports a material definition for methane as an ideal gas and then defines a simple process model that simply specifies a methane flow by temperature, pressure and volume flow.

```
from simu import Model
from .material_factory import ch4_ideal_gas

class Source(Model):
    """A model of a methane source"""

    def interface(self):
        self.parameters.define("T", 25, "degC")
        self.parameters.define("p", 1, "bar")
        self.parameters.define("V", 10, "m^3/hr")

    def define(self):
        src = self.materials.create_flow("source", ch4_ideal_gas)
        self.residuals.add( "T", self.parameters["T"] - src["T"], "K")
        self.residuals.add( "p", self.parameters["p"] - src["p"], "bar")
        self.residuals.add( "V", self.parameters["V"] - src["V"], "m^3/h")


```

This model can be solved as follows:

```
from simu import NumericHandler, SimulationSolver
from simu.examples.material_model import Source

numeric = NumericHandler(Source.top())
solver = SimulationSolver(numeric)
solver.solve()
```

The solver applies the *Newton-Raphson* method, whereas the system matrix (here pathetically 3x3) is obtained using [CasADi](https://web.casadi.org). The linearized systems are solved using the multicore *Intel oneAPI* solver provided by the [PyPardiso](https://pypi.org/project/pypardiso/) package:
```
Iter   LMET   Alpha   Time   Limit on bound Max residual
----- ----- ------- ------ ---------------- ------------
    0   9.0    0.83   0.07 source/n/Methane            T
    1   8.2       1   0.07                             T
    2   6.4       1   0.07                             V
    3  -7.3       1   0.07                             V
```
Here, the relaxation factor &alpha; is limited in the first iteration for the quantity of methane to remain positive. The `LMET` column shows the maximum logarithmic (base 10) residual error to tolerance ratio. The model is converged when `LMET < 0`. 
