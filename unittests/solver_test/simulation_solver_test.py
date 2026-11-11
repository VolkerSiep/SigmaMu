from typing import cast
from pytest import raises

from simu import (
    NumericHandler, NHKeys, SimulationSolver,
    Quantity, quantity_dict_to_strings)
from simu.core.utilities.residual import ResidualHandler
from simu.core.utilities.errors import NonSquareSystem
from simu.core.utilities.testing import assert_reproduction
from simu.examples.material_model import Source


def test_instantiate():
    numeric = NumericHandler(Source.top())
    _ = SimulationSolver(numeric)


def test_solve(sim_result):
    assert len(sim_result.iterations) < 5


def test_solve_res_small(sim_result):
    res = sim_result.properties[NHKeys.VECTORS][NHKeys.RESIDUALS]
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
    solver = SimulationSolver(numeric, output=None)
    param = solver.model_parameters.to_nested_map()[NHKeys.MODEL_PARAMS]
    assert_reproduction(quantity_dict_to_strings(param))


def test_change_parameters():
    numeric = NumericHandler(Source.top())
    solver = SimulationSolver(numeric, output=None)
    spec = Quantity(2, "MPa")
    solver.model_parameters.set_by_path([NHKeys.MODEL_PARAMS, "p"], spec)
    res = solver.solve(max_iter=5)
    p = res.properties["thermo_props"]["source"]["p"]
    assert abs(p - spec) < Quantity(1e-7, "bar")


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

