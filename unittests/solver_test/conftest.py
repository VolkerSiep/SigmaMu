from pathlib import Path
from pytest import fixture
from yaml import safe_load
from simu.core.utilities.types import Map

@fixture
def thermo_fit_configuration():
    file = Path(__file__).parent / "thermo_fit_definition.yml"
    with file.open() as file:
        return safe_load(file)


@fixture(scope="session")
def contribution_context_stub():
    class NumericHandlerStub:
        @property
        def parameters(self) -> Map[str]:
            result = {"T": "K", "p": "bar", "x": "-", "y": "-", "w": "-"}
            return {f"process.{n}": u for n, u in result.items()}

        @property
        def properties(self) -> Map[str]:
            names = (
                [f"process.dmu_norm/{n}" for n in ("H2O", "CO2")] +
                [f"process.{n}" for n in ("p", "y")]
            )
            return {n: ("bar" if n == "process.p" else "") for n in names}

    return {n: NumericHandlerStub() for n in ("vle_fit", "vle_eval_p")}