from collections.abc import Sequence
from pathlib import Path
from pytest import fixture
from yaml import safe_load

@fixture
def thermo_fit_configuration():
    file = Path(__file__).parent / "thermo_fit_definition.yml"
    with file.open() as file:
        return safe_load(file)


@fixture(scope="session")
def contribution_context_stub():
    class NumericHandlerStub:
        @property
        def parameter_names(self) -> Sequence[str]:
            return [f"process.{n}" for n in ("T", "p", "x", "y", "w")]

        @property
        def parameter_units(self) -> Sequence[str]:
            return ["K", "bar", "", "", ""]

        @property
        def property_names(self) -> Sequence[str]:
            return (
                [f"process.dmu_norm/{n}" for n in ("H2O", "CO2")] +
                [f"process.{n}" for n in ("p", "y")]
            )
        @property
        def property_units(self) -> Sequence[str]:
            return ["", "", "bar", ""]

    return {n: NumericHandlerStub() for n in ("vle_fit", "vle_eval_p")}