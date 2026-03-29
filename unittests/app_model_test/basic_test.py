from simu import NumericHandler
from simu.core.utilities.testing import assert_reproduction

def test_create_species_balance(species_balance_example_model):
    num = NumericHandler(species_balance_example_model.top())

def test_species_balance(species_balance_example_model_stub):
    model = species_balance_example_model_stub.top()
    bal = model.hierarchy.handler["n_bal"]
    res = {name: expr.value.m for name, expr in bal.residuals.items()}
    assert_reproduction(res)

