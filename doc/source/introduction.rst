Introduction
============
Naming
------
The name ``SigmaMu`` [sɪɡ.mə mjuː] is obviously based on the Greek letters ``sigma`` and ``mu``:

.. math::

    \Huge{\Sigma\mu}

In short, you might read *simu*, which superficially can stand for *simulation*, and as such also gives the name to the python package (``simu``). The Greek letters can also be interpreted as the sum over chemical potentials. ``SigmaMu`` is built on the concept of canonical modeling, that is to represent the system's constraints in their natural form. The chemical potential, being the partial derivative of Gibbs energy with respect to molar quantities, is the natural way to describe chemical equilibrium, phase equilibrium, and driving forces.
The constraint is then expressed as a stoichiometric sum of chemical potentials, weighted by the stoichiometry of the process. One can write

.. math::

    \sum_{i\in\mathrm{reactants}}\hspace{-1em} \nu_i\mu_i = 0

This is a much easier and modern approach than dealing with Henry, fugacity and activity coefficients, and other constructs of thermodynamics that have been invented to facilitate simplified calculations by pencil, constrained to the two dimensions that old fashion paper could offer.

Logo and history
----------------
The logo pictures a goose riding a python for obvious reasons, that is to get from A to B efficiently. The real question is: Why a goose, and why is it not afraid of the snake?

.. image:: figures/simu_logo.jpeg
    :width: 400
    :align: center

The name of the goose became *Muriel* at a later point in time. Muriel was first alone, and I met it in Canada on a business trip. It nested on a site that was scheduled for excavation. Due to protection laws, the project got delayed by 2 weeks and everybody was angry with her. At that time, I was writing in-house code for a custom process simulation in another project. My wife drew the first logo, which I included in the report - a bit to tease the engineers on site.

.. image:: figures/goose.png
    :align: center

At first, Muriel was afraid of the snake and didn't trust python for being its modelling language. That in-house model was built on my phd work :cite:p:`Siepmann_2006` and written with help of python, but defined in text files with their own format. After the project, Muriel started to trust the snake more and I decided to convert the software into a pure python library and thus make python itself the modelling language. This is when Muriel started to ride the python, and my wife created above logo.

When you now look back at the letter combination :math:`\Sigma\mu`, with wrong glasses and some imagination you can picture :math:`\Sigma` being the python and :math:`\mu` being Muriel.
Her name in a way starts with this letter, and *Muriel* is a name with origins in the Irish and Gaelic language, describing the shining sea. In Norway, the expression *morild* (sea fire) stands for the effect of sea being illuminated by plankton that is excitated by movements in the water - a nice experience when taking a late-summer night swim in the fjords.

In other words: **Muriel brings light into the dark when you are up to the neck in deep water.**


So what is SiMu?
----------------
``SigmaMu`` is a python library for chemical process modeling. It supports and strongly encourages **hierarchical modeling** of reusable modules. A module can be for instance a phase interface, a set of balance equations, represent a tray in a column, a unit operation, a process section or the entire process.

By being a **pure python library** ``SigmaMu`` supports and strongly encourages professional software development techniques for quality assurance. This includes unit testing, use of ``git`` (or other distributed versioning and control systems), and in-code documentation with `sphinx`_ and ``mathjac``. As such, models can be more maintainable and reusable, and one can work in teams on the development - just like in real software development.

``SigmaMu`` uses `casadi`_ to establish a symbolic representation of the entire model, then to be solved in an **equation oriented** manner. A tailor-made solver interacts with the model to adapt the step-length not to leave the model domain. With accurate derivatives being available, process optimization, parameter fitting, and data reconciliation are strong sides of the software.

Equation-oriented process modeling software brings the challenge of **good fault analysis** if the model is not posed correctly. ``SigmaMu`` analyses singular system matrices for the set of likely process variables involved in the fault. Starting values can easily be stored in human-readable ``json`` or ``yaml`` files .

``SigmaMu`` is as such **by design not interactive** and has **no graphical user interface**. Python however enables to interface the model in uncountable ways, for instance towards a spread-sheet calculator, with an own custom GUI (using a GUI toolkit library), or as a web-service towards an IoT platform. ``SigmaMu`` is intended to run on MS Windows, Linux and MAC OS. Efficient deployment of models in Linux docker containers as web-servers is one attractive way to go.

``SigmaMu`` is a **custom modeling tool** and is as such not shipped with a wide library of generic unit operations. The concept is to develop the models from equation basis, but then reuse the modules to maximal degree.

``SigmaMu`` is developed for industrial use-cases. As such, **proprietary data is to be handled and encapsulated**. ``SigmaMu`` is designed to be slim by nature, but new python packages can be built on top to extend the basic functionality for instance with custom thermodynamic models - or just model parameters, reaction chemistry data, and in-house process models.

``SigmaMu`` is originally designed for lumped process modeling (flow-sheet modeling). Distributed systems must be discretized manually by the developer. However, 2D distributed systems have been successfully modeled already with the predecessor, such as reactive absorption of gas into a falling laminar film, based on a Gibbs excess model to describe the liquid phase properties.
