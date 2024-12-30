from symtable import Function
from typing import Callable
from copy import deepcopy


from casadi import MX, jacobian, jtimes, Function
from numpy import array, min, argmax, abs, squeeze
from scipy.sparse import csc_array
from pypardiso import spsolve

from ..model.numeric import NumericHandler
from ..utilities import Quantity, QFunction
from ..utilities.types import NestedMutMap
from ..utilities.configurable import Configurable

_VEC, _STATE = NumericHandler.VECTORS, NumericHandler.STATE_VEC
_RES, _BOUND = NumericHandler.RES_VEC, NumericHandler.BOUND_VEC


class SimulationSolver(Configurable):
    r"""
    The simulation solver assumes both thermodynamic and model parameters to
    be constant
    """
    # noinspection PyUnusedLocal
    # Options are parsed via inspection
    def __init__(self, model: NumericHandler, *,
                 max_iter: int = 30,
                 gamma: float = 0.9,
                 verbose: bool = False,
                 call_back_iter: Callable[[int], bool] = None):
        """

        """
        super().__init__(exclude=["model"])
        self._model = model

        # store arguments (parameters) so the user can change them
        args = deepcopy(model.arguments)
        # store size of state
        self.__state_size = args[_VEC][_STATE].magnitude.size()[0]
        # user shall not think that putting a state here has any effect
        del args[_VEC][_STATE]
        self.__model_parameters : NestedMutMap[Quantity] = args

    def solve(self):
        gamma = self.options["gamma"]
        max_iter = self.options["max_iter"]
        funcs = self._prepare_functions()
        x = self.initial_state

        for iteration in range(max_iter):
            # evaluate system (matrix and rhs)
            r, dr_dx = funcs["f_r"](x)
            r = squeeze(array(r))
            dr_dx = csc_array(dr_dx)

            # assess error
            max_err_idx = argmax(abs(r))
            max_err = abs(r[max_err_idx])
            print(max_err)
            if max_err < 1:
                break

            # calculate full update
            dx = -spsolve(dr_dx, r)

            # find relaxation factor
            a = squeeze(array(funcs["f_b"](x, dx)))
            a = a[0 < a]
            alpha = min([1, min(a) * gamma]) if len(a) else 1

            # apply update
            x = x + alpha * dx
        else:
            msg = f"Model did not converge after {max_iter} iterations"
            raise ValueError(msg)

    def _prepare_functions(self):
        # prepare
        #  - a casadi MX function x -> (r, dr/dx)
        #  - a casadi MX function: (x, dx) -> (a_i = b_i / (db_i/dx_j) * dx_j)
        # prepare a QFunction x -> (y_m, y_t)
        param = deepcopy(self.__model_parameters)

        param[_VEC][_STATE] = (x := Quantity(MX.sym("x", self.__state_size)))
        res = self._model.function(param, squeeze_results=False)[_VEC]
        r, b = res[_RES], res[_BOUND]
        dx = Quantity(MX.sym("x", self.__state_size))
        return {
            "f_r": Function("f_r", [x], [r, jacobian(r, x)]),
            "f_b": Function("f_b", [x, dx], [-b / jtimes(b, x, dx)]),
            "f_y": QFunction({"x:": x}, res)
        }




    @property
    def initial_state(self):
        """Freshly extract the initial values from the model. These might have
        been changed after the solver class was instantiated"""
        args = self._model.arguments
        return args[NumericHandler.VECTORS][NumericHandler.STATE_VEC]

    @property
    def model_parameters(self) -> NestedMutMap[Quantity]:
        return self.__model_parameters


    @property
    def _arg_validations(self):
        between = Configurable._validate_between
        return {
            "max_iter": between(1, 10000),
            "gamma": between(0.1, 0.999),
            "call_back_iter": {
                "f": lambda x: x is None or callable(x),
                "msg": "must be callable",
                "replace_none": lambda iteration: True
            }
        }
