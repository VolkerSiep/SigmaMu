from pydantic import BaseModel
from pytest import fixture

from simu.core.utilities.pydantic_types import PAmount, PPressure, PTemperature
from simu.core.utilities.types import Map

@fixture(scope="session")
def state_type():
    class State(BaseModel):
        t: PTemperature
        p: PPressure
        n: Map[PAmount]
    return State