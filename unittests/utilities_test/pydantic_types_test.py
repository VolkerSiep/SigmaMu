from pydantic import ValidationError
from pytest import raises, mark

from simu import Quantity


def test_quantity_annotations(state_type):
    state = state_type.model_validate({
        "t": "300 K",
        "p": "2 Pa",
        "n": {"A": "1 mol"},
    })
    assert isinstance(state.t, Quantity)
    assert isinstance(state.p, Quantity)
    assert isinstance(state.n["A"], Quantity)


@mark.parametrize("d", (
    {"t": "-1 K", "p": "2 Pa", "n": {"A": "1 mol"}},
    {"t": "300 K", "p": "0 Pa", "n": {"A": "1 mol"}},
    {"t": "300 K", "p": "2 Pa", "n": {"A": "-1 mol"}},
))
def test_quantity_annotations_failing(state_type, d):
        with raises(ValidationError):
            state_type.model_validate(d)
