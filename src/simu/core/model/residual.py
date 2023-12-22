from dataclasses import dataclass
from collections.abc import Iterator

from ..utilities.quantity import Quantity
from ..utilities.types import Map
from ..utilities.errors import DimensionalityError


@dataclass
class Residual:
    """Class representing a single residual"""
    value: Quantity
    tolerance: Quantity


ResidualProxy = Map[Residual]


class ResidualHandler(ResidualProxy):
    """This class, being instanciated as the :attr:`Model.residuals` attribute,
    allows to define residuals, i.e. process constraints."""

    def __init__(self):
        self.__residuals = {}

    def add(self, name: str, residual: Quantity,
            tol_unit: str, tol: float = 1e-7):
        """Define a residual"""
        if name in self.__residuals:
            raise KeyError(f"Residual {name} already defined")
        try:
            residual.to(tol_unit)
        except DimensionalityError:
            msg = f"Incompatible tolerance unit in residual {name}"
            raise DimensionalityError(residual.units, tol_unit, extra_msg=msg)

        # eliminate impact of offset in units like degC and barg
        tolerance = Quantity(tol, tol_unit) - Quantity(0.0, tol_unit)
        self.__residuals[name] = Residual(residual, tolerance)

    def __getitem__(self, key: str) -> Residual:
        return self.__residuals[key]

    def __len__(self) -> int:
        return len(self.__residuals)

    def __iter__(self) -> Iterator[str]:
        return iter(self.__residuals)

    def create_proxy(self) -> ResidualProxy:
        """Create a proxy object"""
        return self
