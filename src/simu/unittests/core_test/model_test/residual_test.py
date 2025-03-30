from pytest import raises

from simu.core.utilities import Quantity
from simu.core.utilities.errors import DimensionalityError
from simu.core.utilities.residual import ResidualHandler


def test_create_residuals():
    handler = ResidualHandler()
    res = Quantity("3 m")   # this could never work, but not important here
    handler.add("Res 1", res, "km")
    handler.add("Res 2", res, "cm", tol=1e-5)
    tol = handler["Res 2"].tolerance
    assert f"{tol:~}" == "1e-05 cm"


def test_create_residual_inconsistent_unit():
    handler = ResidualHandler()
    res = Quantity("3 m")  # this could never work, but not important here
    with raises(DimensionalityError):
        handler.add("My Residual", res, "kJ")


def test_create_residual_same_name():
    handler = ResidualHandler()
    res = Quantity("3 m")  # this could never work, but not important here
    handler.add("My Residual", res, "km")
    with raises(KeyError):
        handler.add("My Residual",  2 * res, "cm")