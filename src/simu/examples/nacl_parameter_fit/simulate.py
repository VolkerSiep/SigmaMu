from model import PSatModel
from simu import NumericHandler, SimulationSolver, NHKeys


def main():
    model = NumericHandler(PSatModel.top())
    solver = SimulationSolver(model)
    result = solver.solve()
    thermo_props = result.properties[NHKeys.THERMO_PROPS]
    print(f"{thermo_props["liquid"]["p"].to("bar"):.3f~}")
    print(f"{thermo_props["liquid"]["mu"]["H2O"].to("kJ/mol"):.2f~}")
    print(f"{thermo_props["gas"]["mu"]["H2O"].to("kJ/mol"):.2f~}")


if __name__ == '__main__':
    main()