from pathlib import Path
from yaml import safe_load
from simu import ThermoFitEvaluator, NumericHandler
from model import PSatModel
from pandas import DataFrame
from matplotlib import pyplot

THIS_DIRECTORY = Path(__file__).parent
THERMO_FIT_DEFINITION_FILE = THIS_DIRECTORY / "thermo_fit_definition.yml"
WASHBURN_H2O_FILE = THIS_DIRECTORY / "Washburn_1928_h2o.yml"
WASHBURN_ALL_FILE = THIS_DIRECTORY / "Washburn_1928_all.yml"


def load_definition():
    with THERMO_FIT_DEFINITION_FILE.open() as f:
        data = safe_load(f)
    with WASHBURN_H2O_FILE.open() as f:
        data["datasets"]["washburn_h2o"] = safe_load(f)
    with WASHBURN_ALL_FILE.open() as f:
        data["datasets"]["washburn_all"] = safe_load(f)
    return data

def evaluate(models, definition):
    evaluator = ThermoFitEvaluator(models)
    result = evaluator.solve(definition)
    data = result["all"].results
    df = DataFrame(data.data, columns=data.columns)

    baseline = df[df["w_nacl"] == 0.0][["T", "p_meas"]].rename(
        columns={"p_meas": "p_meas_pure"}
    )
    df = df.merge(baseline, on="T", how="left")
    df["p_meas_red"] = df["p_meas"] / df["p_meas_pure"]
    df["p_calc_red"] = df["p_calc"] / df["p_meas_pure"]

    for (t, data), c in zip(df.groupby("T"), "krbgmckrbgmc"):
        pyplot.plot(data["w_nacl"], data["p_meas_red"], f"{c}.")
        pyplot.plot(data["w_nacl"], data["p_calc_red"], f"{c}-", label=f"T = {t} degC")
    pyplot.grid()
    pyplot.legend(loc="best")
    pyplot.xlabel("w(NaCl) [%]")
    pyplot.ylabel("Pressure ratio [-]")
    pyplot.show()


def main():
    definition = load_definition()
    models = {"p_sat": NumericHandler(PSatModel.top())}
    evaluate(models, definition)


if __name__ == '__main__':
    main()