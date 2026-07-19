from pathlib import Path
from yaml import safe_load
from common import define_models, load_definition
from simu import ThermoFitEvaluator, parse_quantities_in_struct
from pandas import DataFrame
from matplotlib import pyplot


PARAM_FILE = Path(__file__).parent / "parameters_fit.yml"

def plot(axes, result, style, alpha, legend):
    colors = "kbgrmcb"
    for name, ax in zip(result, axes.T):
        ax[0].set_title(name)
        phase = result[name].results
        df = DataFrame(phase.data, columns=list(phase.columns))
        for (p, data), c in zip(df.groupby("p"), colors):
            label=f"{p}" if legend else None
            ax[0].plot(data["T"], data["f_mu"], f"{c}{style}", alpha=alpha, label=label)
            ax[1].plot(data["T"], data["q_h"], f"{c}{style}", alpha=alpha)
            ax[2].plot(data["T"], data["q_rho"], f"{c}{style}", alpha=alpha)


def main():
    definition = load_definition()
    models = define_models()

    fig, axes = pyplot.subplots(
        ncols=2, nrows=3, figsize=(9, 7), sharex=True, layout="tight"
    )

    evaluator = ThermoFitEvaluator(models, gamma=0.7)
    result = evaluator.solve(definition)
    plot(axes, result, style="--", alpha=0.5, legend=False)

    try:
        with PARAM_FILE.open() as file:
            param = parse_quantities_in_struct(safe_load(file))
    except FileNotFoundError:
        pass
    else:
        result = evaluator.solve(definition, param)
        plot(axes, result, style="-", alpha=1.0, legend=True)

    for r in range(3):
        axes[r][0].grid()
        axes[r][1].grid()

    axes[2][0].set_xlabel(r"Temperature [$^\circ$C]")
    axes[2][1].set_xlabel(r"Temperature [$^\circ$C]")
    axes[0][0].set_ylabel(r"Deviation factor $\exp \frac{\Delta \mu}{R\,T}$ [-]")
    axes[1][0].set_ylabel(r"Deviation $\frac{\Delta H}{R\,T}$ [-]")
    axes[2][0].set_ylabel(r"Deviation $\frac{\Delta \varrho}{\varrho}$ [-]")
    axes[0][0].legend(loc="best", title="Pressure [bar]", ncol=2)

    # pyplot.savefig("srk_rho_poly.png")
    pyplot.show()


if __name__ == '__main__':
    main()