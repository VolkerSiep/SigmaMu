from typing import cast
from pytest import fixture, raises

from simu import NumericHandler, SimulationSolver, Quantity
from simu.core.utilities.residual import ResidualHandler
from simu.core.utilities.errors import NonSquareSystem
from simu.core.utilities.qstructures import quantity_dict_to_strings
from simu.core.utilities import assert_reproduction
from simu.examples.material_model import Source


def test_instantiate():
    numeric = NumericHandler(Source.top())
    _ = SimulationSolver(numeric)


def test_solve(sim_result):
    assert len(sim_result.iterations) < 5


def test_solve_res_small(sim_result):
    res = sim_result.properties[NumericHandler.VECTORS][NumericHandler.RES_VEC]
    assert (abs(res.m_as("")) < 1).all()


def test_reproduce_result(sim_result):
    n = sim_result.properties["thermo_props"]["source"]["n"]["Methane"]
    assert abs(n.m_as("mol/s") - 0.112054293180843) < 1e-7


def test_non_square():
    """add another residual and thus make system non-square"""
    model = Source().create_proxy().finalise()
    residuals = cast(ResidualHandler, model.residuals) # need to add one more
    residuals.add("Q", model.residuals["T"].value, "K")
    numeric = NumericHandler(model)
    with raises(NonSquareSystem):
        _ = SimulationSolver(numeric)


def test_model_parameters():
    numeric = NumericHandler(Source.top())
    solver = SimulationSolver(numeric, output="None")
    param = quantity_dict_to_strings(solver.model_parameters["model_params"])
    assert_reproduction(param)


def test_invalid_option_constructor():
    numeric = NumericHandler(Source.top())
    with raises(ValueError):
        _ = SimulationSolver(numeric, max_iter=-3)


def test_valid_option_set_option():
    numeric = NumericHandler(Source.top())
    solver = SimulationSolver(numeric)
    solver.set_option("max_iter", 20)


def test_invalid_option_set_option():
    numeric = NumericHandler(Source.top())
    solver = SimulationSolver(numeric)
    with raises(ValueError):
        solver.set_option("max_iter", -3)


def test_change_parameters():
    numeric = NumericHandler(Source.top())
    solver = SimulationSolver(numeric, output="None")
    param = solver.model_parameters["model_params"]
    param["p"] = Quantity(2, "MPa")
    res = solver.solve(max_iter=5)
    p = res.properties["thermo_props"]["source"]["p"]
    assert abs(p - param["p"]) < Quantity(1e-7, "bar")


    # The following code can be start of providing a report unit system
    # from simu import flatten_dictionary
    # props = flatten_dictionary(res.result["thermo_props"])
    # default_units = {
    #     "[temperature]": "degC",
    #     "[length]": "m",
    #     "[volume]": "m^3",
    #     "[length] ** 3 / [time]": "m^3/hr",  # Volumetric flow
    #     "[time]": "s",
    #     "[mass] / [length] / [time] ** 2": "bar"
    # }
    # for k, q in props.items():
    #     try:
    #         unit = default_units[str(q.dimensionality)]
    #     except KeyError:
    #         pass  # print(k, q.dimensionality)
    #     else:
    #         q = q.to(unit)
    #     print(f"{k:20s} {q:.3g~}")


@fixture(scope="module")
def sim_result():
    numeric = NumericHandler(Source.top())
    solver = SimulationSolver(numeric, output="None")
    return solver.solve()
