from pathlib import Path
from pytest import fixture
from yaml import safe_load

@fixture(scope="session")
def example_thermo_fit_configuration():
    file = Path(__file__).parent / "thermo_fit_definition.yml"
    with file.open() as file:
        return safe_load(file)
