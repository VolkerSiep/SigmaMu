from pathlib import Path
from importlib import import_module

from simu.examples.h2o_fit.evaluation import main as h2o_fit_main

CURRENT_DIR = Path(__file__).parent

def run_example_plot(package: str, module_name: str, figure_name: str,
                     dependent: list[str] | None = None):
    module = import_module(f".{module_name}", package=package)
    source = Path(module.__file__)
    figure = CURRENT_DIR / figure_name


    if figure.exists() and figure.stat().st_mtime >= source.stat().st_mtime:
        print(f"Skipping '{module_name}', as '{figure_name}' is up to date.")
        return
    print(f"Running '{module_name}', as '{figure_name}' is out of date.")
    for d in dependent:
        print(f"  First running dependent script '{d}'")
        dep = import_module(f".{d}", package=package)
        dep.main()
    module.main(figure)


def main():
    run_example_plot("simu.examples.h2o_fit", "evaluation", "h2o_vle_fit.png",
                     dependent=["fit"])

if __name__ == '__main__':
    main()