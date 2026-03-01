from pprint import pprint
from simu import NumericHandler, SimulationSolver, Quantity, NHKeys
from process import SteamGeneration


def main():
    numeric = NumericHandler(SteamGeneration.top(), port_properties=False)
    solver = SimulationSolver(numeric)
    result = solver.solve()
    print_result(result)

    temp = Quantity(35, "degC")
    solver.model_parameters[NHKeys.MODEL_PARAMS]["condenser"]["temperature"] = temp
    result = solver.solve()
    print_result(result)


def print_result(result):
    props = result.properties
    print("Stream properties:")
    excluded = "A G Mw S U V m mu mw n w x".split()
    thermo_props = filter_results(props[NHKeys.THERMO_PROPS], excluded)
    pprint(thermo_props)
    print()
    print("Model properties:")
    pprint(filter_results(props[NHKeys.MODEL_PROPS]))


def filter_results(results, excluded_keys=None):
    excluded_keys = [] if excluded_keys is None else excluded_keys
    try:
        items = results.items()
    except AttributeError:
        return f"{results:.5g~}"
    return {key: filter_results(value, excluded_keys)
            for key, value in items if key not in excluded_keys}


if __name__ == '__main__':
    main()

