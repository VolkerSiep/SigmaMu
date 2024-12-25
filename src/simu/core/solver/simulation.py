from typing import Callable

# todo: change into relative import
from simu.core.utilities.configurable import Configurable
from simu import NumericHandler

class SimulationSolver(Configurable):
    # noinspection PyUnusedLocal
    # Options are parsed via inspection
    def __init__(self, max_iter: int = 30, gamma: float = 0.9, *,
                 verbose: bool = False,
                 call_back_iter: Callable[[int], bool] = None):
        """

        """
        super().__init__()
        print(self.options["call_back_iter"](2))


    def solve(self, model: NumericHandler, **kwargs):
        pass

    @property
    def _arg_validations(self):
        return {
            "max_iter": {
                "f": lambda x: 1 <= x <= 10000,
                "msg": "must be between 1 and 10000",
            },
            "gamma": {
                "f": lambda x: 0.1 < x < 0.999,
                "msg": "must be between 0.1 and 0.999"
            },
            "call_back_iter": {
                "f": lambda x: x is None or callable(x),
                "msg": "must be callable",
                "replace_none": lambda iteration: True
            }
        }


s = SimulationSolver(max_iter=30, verbose=False)