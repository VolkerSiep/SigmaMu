from pathlib import Path
from yaml import safe_load
from pandas import DataFrame
from simu import ThermoFitEvaluator, parse_quantities_in_struct
from simu.examples.plotting import pyplot
from simu.examples.h2o_fit.common import define_models, load_definition


CURRENT_DIR = Path(__file__).parent
PARAM_FILE = CURRENT_DIR / "parameters_fit.yml"

def plot(axes, result, legend, **style):
    for name, ax in zip(result, axes.T):
        ax[0].set_title(name)
        phase = result[name].results
        df = DataFrame(phase.data, columns=list(phase.columns))
        for p, data in df.groupby("p"):
            label=f"{p}" if legend else None
            l, = ax[0].plot(data["T"], data["f_mu"] - 1, label=label, **style)
            c = l.get_color()
            ax[1].plot(data["T"], data["q_h"], color=c, **style)
            ax[2].plot(data["T"], data["q_rho"], color=c, **style)


def read_parameters():
    try:
        with PARAM_FILE.open() as file:
            return parse_quantities_in_struct(safe_load(file))
    except FileNotFoundError:
        print("Data fit has not been run yet - skipping plot.")
        return None


def main(figure_file: Path | None = None):
    definition = load_definition()
    models = define_models()

    fig, ax = pyplot.subplots(
        ncols=2, nrows=3, figsize=(9, 7), sharex=True, layout="tight"
    )

    evaluator = ThermoFitEvaluator(models, gamma=0.7)
    result = evaluator.solve(definition)
    plot(ax, result, legend=False, linestyle="--")

    if (param := read_parameters()) is not None:
        result = evaluator.solve(definition, param)
        plot(ax, result, legend=True)

    ax[2][0].set_xlabel(r"Temperature [$^\circ$C]")
    ax[2][1].set_xlabel(r"Temperature [$^\circ$C]")
    ax[0][0].set_ylabel(r"Deviation factor $\exp \frac{\Delta \mu}{R\,T} - 1$ [-]")
    ax[1][0].set_ylabel(r"Deviation $\frac{\Delta H}{R\,T}$ [-]")
    ax[2][0].set_ylabel(r"Deviation $\frac{\Delta \varrho}{\varrho}$ [-]")
    ax[0][0].legend(loc="best", title="Pressure [bar]", ncol=2)

    if figure_file is None:
        pyplot.show()
    else:
        pyplot.savefig(figure_file)


if __name__ == '__main__':
    main()