from pathlib import Path
from yaml import safe_load
from simu import ThermoFitEvaluator, parse_quantities_in_struct
from pandas import DataFrame
from matplotlib import pyplot

from common import load_definition, define_models


PARAM_FILE = Path(__file__).parent / "parameters_fit.yml"


def evaluate(evaluator, param=None):
    definition = load_definition()
    result = evaluator.solve(definition, param)
    data = result["all"].results
    df = DataFrame(data.data, columns=list(data.columns))

    baseline = df[df["w_nacl"] == 0.0][["T", "p_meas"]].rename(
        columns={"p_meas": "p_meas_pure"}
    )
    df = df.merge(baseline, on="T", how="left")
    df["p_meas_red"] = df["p_meas"] / df["p_meas_pure"]
    df["p_calc_red"] = df["p_calc"] / df["p_meas_pure"]
    return df


def main():
    models = define_models()
    evaluator = ThermoFitEvaluator(models)
    df = evaluate(evaluator)

    try:
        with PARAM_FILE.open() as file:
            param = parse_quantities_in_struct(safe_load(file))
    except FileNotFoundError:
        pass
    else:
        df_fit = evaluate(evaluator, param)
        df = df.merge(df_fit, on=["T", "w_nacl"], suffixes=("", "_fit"))

    for t, data in df.groupby("T"):
        l, = pyplot.plot(data["w_nacl"], data["p_meas_red"], ".")
        pyplot.plot(
            data["w_nacl"], data["p_calc_red_fit"], "-",
            color=l.get_color(), label=f"T = {t} degC"
        )
        pyplot.plot(
            data["w_nacl"], data["p_calc_red"], "--",
            color=l.get_color()
        )
    pyplot.grid()
    pyplot.legend(loc="best")
    pyplot.xlabel("w(NaCl) [%]")
    pyplot.ylabel("Pressure ratio [-]")
    pyplot.show()



if __name__ == '__main__':
    main()