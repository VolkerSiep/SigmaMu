from collections.abc import Container, Iterator

from abc import ABC
from simu import AModel, QuantityDict, Quantity


class MultiNode(AModel, ABC):
    """This abstract model base defines a model with ``num_in`` inlets and
    ``num_out`` outlets.

    The inlet ports are called ``in_01``, ``in_02``, ``in_03``, etc.
    The outlet ports are called ``out_01``, ``out_02``, ``out_03``, etc.

    Subclasses still need to implement the :meth:`~simu.Model.define` method.

    Using this base-class is a good practice to maintain a standard naming of
    inlet and outlet ports for models with generic character and indistinct
    ports, such as in mixers or splitters. Please do not do this for Models
    where ports have individual roles. In that case, individual naming, such as
    ``cold_in`` or ``gas_out`` are preferable.

    Further, the ports defined in this class do not apply specific material
    specifications.

    :param num_in: The number of inlet streams to consider
    :param num_out: The number of outlet streams to consider
    """
    def __init__(self, num_in: int = 1, num_out: int = 1):
        self._num_in, self._num_out = num_in, num_out
        super().__init__()

    def interface(self):
        for name in self.in_names():
            self.md(name)
        for name in self.out_names():
            self.md(name)

    def in_names(self) -> Iterator[str]:
        """A generator to yield the names of all input ports"""
        for i in range(self._num_in):
            yield f"in_{i + 1:02d}"

    def out_names(self) -> Iterator[str]:
        """A generator to yield the names of all output ports"""
        for i in range(self._num_out):
            yield f"out_{i + 1:02d}"


class SpeciesBalance(MultiNode):
    r"""This model balances all species flows between the given inlet and outlet
    streams.

    The model handles streams with individual species sets. Naturally, each
    occurring species must then occur in at least one outlet and one inlet.
    Otherwise, its partial flow would be forced to zero, which is not accepted
    by most thermodynamic models (e.g. the ideal mix term :math:`\ln x_i` would
    diverge to negative infinity.

    :param num_in: The number of inlet streams to consider
    :param num_out: The number of outlet streams to consider
    :param tol_unit: The tolerance unit, representing a typical order of
      magnitude for the involved flows
    """
    def __init__(self, num_in: int = 1, num_out: int = 1,
                 tol_unit: str = "kmol/h", ignore: Container[str] = None):
        self._tol_unit = tol_unit
        self._ignore = [] if ignore is None else ignore
        super().__init__(num_in, num_out)

    def define(self):
        inlets = [self.m[name] for name in self.in_names()]
        outlets = [self.m[name] for name in self.out_names()]

        d_n = sum(m["n"] for m in inlets) -  sum(m["n"] for m in outlets)
        for species, dn_i in d_n.items():
            if not species in self._ignore:
                self.ra(species, dn_i, self._tol_unit)


class ElementBalance(MultiNode):
    """This model balances all element flows between the given inlet and outlet
    stream"""
    def __init__(self, num_in: int = 1, num_out: int = 1,
                 tol_unit: str = "kmol/h", ignore: Container[str] = None):
        self._tol_unit = tol_unit
        self._ignore = [] if ignore is None else ignore
        super().__init__(num_in, num_out)

    def define(self):
        inlets = [self.m[name] for name in self.in_names()]
        outlets = [self.m[name] for name in self.out_names()]

        d_e = sum(m["n_e"] for m in inlets) -  sum(m["n_e"] for m in outlets)
        for element, de_i in d_e.items():
            if not element in self._ignore:
                self.ra(element, de_i, self._tol_unit)



class EnthalpyBalance(MultiNode):
    """This model balances the enthalpy for an adiabatic multi-node.

    :param num_in: The number of inlet streams to consider
    :param num_out: The number of outlet streams to consider
    :param tol_unit: The tolerance unit, representing a typical order of
      magnitude for the involved enthalpy flows
    """
    def __init__(self, num_in: int = 1, num_out: int = 1, tol_unit = "MW"):
        self._tol_unit = tol_unit
        super().__init__(num_in, num_out)

    def interface(self):
        super().interface()
        self.pad("Duty", 0.0, "MW")

    def define(self):
        inlets = [self.m[name] for name in self.in_names()]
        outlets = [self.m[name] for name in self.out_names()]

        d_h = sum(m["H"] for m in inlets) -  sum(m["H"] for m in outlets)
        self.ra("h_balance", d_h + self.pa["Duty"], self._tol_unit)


class PhaseEquilibrium(AModel):
    r"""This model defines phase equilibrium between two streams by equalizing
    temperature, pressure, and chemical potential for all common species.

    If the species sets of both phases differ, no constraints are defined for
    the non-common ones.

    .. math::

        T_1 = T_2,\quad p_1 = p_2,\quad\text{and}\quad
        \mu_{1,i} = \mu_{2,i} \forall i\in \mathbb S_1 \cup \mathbb S_2

    If both materials belong to the same definition, the solution will be
    trivial and hence the equations linear dependent, so do not do that.

    The model does not prescribe the nature of the phases, as long as they are
    independent.

    """
    def interface(self):
        self.md("phase_1")
        self.md("phase_2")

    def define(self):
        p1, p2 = self.m["phase_1"], self.m["phase_2"]
        self.ra("T_eq", p1["T"] - p2["T"], "K")
        self.ra("p_eq", p1["p"] - p2["p"], "bar")

        for s in p1["mu"].keys() & p2["mu"].keys():
            self.ra(f"mu_eq_{s}", p1["mu"][s] - p2["mu"][s], "kJ/mol")


