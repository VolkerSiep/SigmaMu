from simu import NumericHandler
from simu.examples.material_model import Source
from simu.core.utilities import assert_reproduction, SymbolQuantity, jacobian


def test_query_bounds():
    numeric = NumericHandler(Source.top())
    res = numeric.vector_res_names(numeric.BOUND_VEC)
    assert_reproduction(res)


def test_bound_sensitivity():
    numeric = NumericHandler(Source.top())
    args = numeric.arguments
    names = numeric.vector_arg_names(numeric.STATE_VEC)
    state = SymbolQuantity("x", "", names)
    args[numeric.VECTORS][numeric.STATE_VEC] = state
    res = numeric.function(args, squeeze_results=False)
    res = res[numeric.VECTORS][numeric.BOUND_VEC]
    jac = jacobian(res, state).magnitude
    assert_reproduction(str(jac))

    # TODO: this can still be functionality of NumericHandler!?
    #  not really, why should the NumericHandler care what the independent
    #  variables are?
