from collections.abc import Sequence
from pathlib import Path

from pytest import fixture
from yaml import safe_load

from simu import StringDictThermoSource, NumericHandler
from simu.core.solver.thermofit.config import ThermoFitValidationContext


@fixture
def thermo_fit_configuration():
    file = Path(__file__).parent / "thermo_fit_definition.yml"
    with file.open() as file:
        return safe_load(file)


@fixture(scope="session")
def linear_system():
    from numpy import array
    from scipy.sparse import csr_array
    matrix = csr_array([[3.0, 1.0], [1.0, 2.0]])
    rhs = array([9.0, 8.0])
    expected_x = array([2.0, 3.0])
    return matrix, rhs, expected_x


@fixture(scope="session")
def contribution_context_stub():
    class NHStub:
        def __init__(self):
            self._parameters = {"T": "K", "p": "bar", "x": "", "y": "", "w": ""}
            self._properties = {"dmu_norm": {"H2O": "", "CO2": ""},
                                "p": "bar", "y": ""}

        def parameter_unit(self, path: Sequence[str]) -> str:
            if len(path) > 1 or path[0] not in self._parameters:
                raise KeyError(f"'{'.'.join(path)}' not found")
            return self._parameters[path[1]]

        def property_unit(self, path: Sequence[str]) -> str:
            if not path:
                raise KeyError("Empty path")
            res = self._properties
            for p in path:
                res = res[p]
            if not isinstance(res, str):
                raise KeyError("Incomplete path")
            return  res

    return ThermoFitValidationContext(
        model_contexts= {n: NHStub() for n in ("vle_fit", "vle_eval_p")},
        thermo_source= StringDictThermoSource({
                "a": {"b": "300 K", "c": "400 K"}
        })
    )
