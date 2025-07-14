from typing import TYPE_CHECKING
from collections.abc import Sequence
from casadi import SX, vertcat, jacobian, Function
from ..utilities import flatten_dictionary

if TYPE_CHECKING:
    from ..thermo import ThermoFrame, InitialState
    from ..utilities import Quantity
    from ..utilities.structures import NestedMap


class VolumeSolver:
    def __init__(self, frame, parameters: NestedMap[Quantity],
                 initial_state: InitialState):
        x = SX.sym("x", len(frame.species) + 2)
        z = SX.sym("z", len(frame.species) + 2)  # initial state (T, p, n)
        result = frame(x, parameters, squeeze_results=False)
        props = result["props"]

        t_res = (props["T"] / initial_state.temperature).to("").m - 1
        p_res = (props["p"] / initial_state.pressure).to("").m - 1
        n_res = (props["n"] / initial_state.mol_vector).to("").m - 1

        def unpack(name):
            flat = flatten_dictionary(result.get("name", {}))
            return [v.m for v in flat.values()]

        bounds = vertcat(*unpack("bounds"))
        residuals = vertcat(*[t_res, p_res, n_res])

        db_dx = jacobian(bounds, x)
        dr_dx = jacobian(residuals, x)
        self.__function = Function("f", [x], [residuals, bounds, dr_dx, dr_dx],
                                   ["x"], ["r", "b", "dr_dx", db_dx])


    def __call__(self, estimate: Sequence[float]) -> Sequence[float]:
        return estimate
        # TODO: refine result by iteration here.
        #  1. make flat state x as SX, call function with paramters as above.
        #  2. concatenate y = [res["T"]/state.temperature, p/p_ini, n/n_ini]
        #  3. J = jacobain(y, x)
        #  4. I might also need jacobian of the bounds with respect to state.
        #  5. And I need the residuals actually!
        #  5. solve system, limiting steps to stay within bounds.
        return internal_estimate