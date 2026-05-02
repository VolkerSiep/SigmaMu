from pathlib import Path

from pytest import fixture
from yaml import safe_load

from simu import StringDictThermoSource
from simu.core.solver.thermofit import ThermoFitValidationContext
from simu.core.utilities.types import Map


@fixture
def thermo_fit_configuration():
    file = Path(__file__).parent / "thermo_fit_definition.yml"
    with file.open() as file:
        return safe_load(file)


@fixture
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
        @property
        def parameters(self) -> Map[str]:
            result = {"T": "K", "p": "bar", "x": "", "y": "", "w": ""}
            return {f"process.{n}": u for n, u in result.items()}

        @property
        def properties(self) -> Map[str]:
            names = (
                [f"process.dmu_norm/{n}" for n in ("H2O", "CO2")] +
                [f"process.{n}" for n in ("p", "y")]
            )
            return {n: ("bar" if n == "process.p" else "") for n in names}

    return ThermoFitValidationContext(
        model_contexts= {n: NHStub() for n in ("vle_fit", "vle_eval_p")},
        thermo_source= StringDictThermoSource({
                "a": {"b": "300 K", "c": "400 K"}
        })
    )
