from simu import NumericHandler, SimulationSolver
from simu.examples.material_model import Source
from simu.core.utilities import assert_reproduction




def test_instantiate():
    numeric = NumericHandler(Source.top())
    _ = SimulationSolver(numeric)


def test_solve():
    numeric = NumericHandler(Source.top())
    solver = SimulationSolver(numeric, output=None)
    res = solver.solve()
    assert len(res.iterations) < 5
    n = res.result["thermo_props"]["source"]["n"]["Methane"]
    assert abs(n.to("mol/s").magnitude - 0.112054293180843) < 1e-7

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
    #         print(k, q.dimensionality)
    #     else:
    #         q = q.to(unit)
    #     print(f"{k:20s} {q:.3g~}")
