from dataclasses import dataclass
from collections.abc import Iterator
from casadi import DM

from simu.core.utilities.quantity import Quantity
from simu.core.utilities.types import Map
from simu.core.utilities.errors import DimensionalityError


@dataclass
class Residual:
    """Class representing a single residual.

    The post-constructor will raise a ``DimensionalityError`` if the
    ``tolerance`` unit is not compatible with the ``value`` unit.

    A possible offset in the tolerance unit, such as in degrees Celsius or gauge
    pressures, will be eliminated automatically. That is:

    >>> from simu import SymbolQuantity as S, Quantity as Q
    >>> val = S("r", "barg")
    >>> r1 = Residual(val, Q(1e-7, "barg"))
    >>> r2 = Residual(val, Q(1e-7, "bar"))
    >>> f"{r1.tolerance - r2.tolerance:.15f~}"
    '0.000000000000000 bar'

    *(there is a silly rounding error, hence the trick with the formatting to
    make this doctest work)*
    """
    value: Quantity
    """The value as a symbolic quantity"""
    tolerance: Quantity
    """The absolute tolerance of the residual as a numeric quantity in the same
    physical dimensions as the value"""

    def __post_init__(self):
        tol_unit = self.tolerance.units
        try:
            self.value.to(tol_unit)
        except DimensionalityError:
            msg = f"Incompatible tolerance unit in residual"
            raise DimensionalityError(self.value.units, tol_unit, extra_msg=msg)

        # eliminate impact of offset in units like degC and barg
        self.tolerance -= Quantity(0.0, tol_unit)


ResidualProxy = Map[Residual]
"""This is just a dictionary, mapping the names (*str*) to the :class:`Residual`
objects"""

class ResidualHandler(ResidualProxy):
    """This class, being instantiated as the :attr:`simu.Model.residuals`
    attribute, allows to define residuals, i.e. process constraints."""

    def __init__(self):
        self.__residuals = {}

    def add(self, name: str, residual: Quantity,
            tol_unit: str, tol: float = 1e-7):
        """Define a residual, to approach zero withing its tolerance when the
        model is solved.

        :param name: A unique name within the local context
        :param residual: The residual to be zero.
        :param tol_unit: The unit of the tolerance, which must be compatible
          with the unit of the residual quantity.
        :param tol: The tolerance value. A normal recommendation is to give
          ``1e-7`` for a typical unit of measurement, which for extensive
          variables is largely problem size dependent. An industrial flow is
          typically characterized by ``kg/h``, while a lab scale experiment
          can operate with ``g/h``.
        """
        if name in self.__residuals:
            raise KeyError(f"Residual {name} already defined")
        self.__residuals[name] = Residual(residual, Quantity(DM(tol), tol_unit))

    def __getitem__(self, key: str) -> Residual:
        return self.__residuals[key]

    def __len__(self) -> int:
        return len(self.__residuals)

    def __iter__(self) -> Iterator[str]:
        return iter(self.__residuals)

    def create_proxy(self) -> ResidualProxy:
        """Create a proxy object"""
        return self
