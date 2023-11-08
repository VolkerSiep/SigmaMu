from simu.core.model import NumericHandler
from .models import *


def test_parameters():
    proxy = SimpleParameterTestModel.top()
    numeric = NumericHandler(proxy)
    args = numeric.function.arg_structure
    assert args['parameters/length'] == 'm'


def test_properties():
    proxy = PropertyTestModel.top()
    numeric = NumericHandler(proxy)
    results = numeric.function.result_structure
    assert results['properties/area'] == 'm ** 2'


def test_residuals():
    proxy = ResidualTestModel.top()
    numeric = NumericHandler(proxy)
    results = numeric.function.result_structure
    assert results['residuals/area'] == "m ** 2"

