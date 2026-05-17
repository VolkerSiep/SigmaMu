from typing import Any
from collections.abc import Sequence
from dataclasses import dataclass

from numpy import squeeze, array, vstack, concatenate, sum, sqrt
from numpy.typing import NDArray
from numpy.linalg import lstsq, norm
from scipy.sparse import csr_array
from casadi import SX, jacobian, jtimes, Function, vertcat
from simu import (
    NumericHandler, NHKeys, AbstractThermoSource, Quantity)
from simu.core.utilities.types import Map, MutMap, NestedMap, NestedMutMap
from simu.core.utilities.errors import NonSquareSystem

from .config import (
    ThermoFitDefinition, ThermoFitSolverConfig, ThermoFitValidationContext,
    ThermoFitContribution, ThermoFitParameter, DataSet, ThermoFitReport
)
from ..common import not_finite, assess_residuals, relax, check_model_square


@dataclass
class _ParameterSymbol:
    name: str
    default_value: Quantity


@dataclass
class _FunctionCollection:
    f_r: Function  # x, p, t -> r, r_x
    f_bx: Function  # x, p, t, dx -> b, a
    f_bt: Function  # x, p, t, dt -> b, a
    f_q: Function  # x, p, t -> q, q_x, q_t, r_t


@dataclass
class _ContributionResult:
    q: Sequence[NDArray]
    dq_dt: Sequence[NDArray]

@dataclass
class _DataPointResult:
    x: NDArray
    dr_dx: csr_array


class ModelContext:
    def __init__(self, model: NumericHandler):
        self._parameters = model.function.arg_structure.get(
            NHKeys.MODEL_PARAMS, {}
        )
        self._properties = model.function.result_structure.get(
            NHKeys.MODEL_PROPS, {}
        )

    def parameter_unit(self, path: Sequence[str]) -> str:
        return self._extract(path, self._parameters)

    def property_unit(self, path: Sequence[str]) -> str:
        return self._extract(path, self._properties)

    @staticmethod
    def _extract(path: Sequence[str], structure: NestedMap[str]) -> str:
        result = structure
        try:
            for p in path:
                result = result[p]
        except (KeyError, TypeError) as e:
            raise KeyError(f"Invalid path: '{'.'.join(path)}'") from e
        if not isinstance(result, str):
            raise KeyError(f"Invalid path: '{'.'.join(path)}'")
        return result

class DataRowConverter:
    def __init__(self, from_uom: Sequence[str], to_uom: Sequence[str]):
        self._from = from_uom
        self._to = to_uom

    def __call__(self, row: Sequence[float]) -> Sequence[float]:
        return [
            Quantity(r, f).to(t).magnitude
            for r, f, t in zip(row, self._from, self._to)
    ]


class ThermoFitContributionWrapper:
    def __init__(self,
                 model: NumericHandler,
                 contribution: ThermoFitContribution,
                 dataset: DataSet,
                 parameters: Map[ThermoFitParameter],
                 config: ThermoFitSolverConfig):
        self._model = model
        self._dataset = dataset
        self._contribution = contribution
        self._funcs = _prepare_functions(model, contribution, parameters)
        self._config = config

        param_uom = [d.uom for d in contribution.data_to_model.values()]
        self._row_converter = DataRowConverter(dataset.uom, param_uom)

        state = model.arguments[NHKeys.VECTORS][NHKeys.STATES].magnitude
        self._states = [state] * len(dataset.data)


    def solve(self, tau: NDArray) -> _ContributionResult:
        states, model, data = self._states, self._model, self._dataset.data
        config = self._config
        q : list[NDArray] = []
        dq_dt: list[NDArray] = []

        for r, row in enumerate(data):
            param = self._row_converter(row)
            try:
                result = self._solve_point(states[r], param, tau)
                states[r] = result.x
            except ValueError:
                continue  # ignore contribution
            q_i, q_x, q_t, r_t = self._funcs.f_q(states[r], param, tau)
            x_t = -config.linear_solver.solve(result.dr_dx, r_t)
            q.append(array(q_i).squeeze(axis=1))
            dq_dt.append(array(q_t + q_x @ x_t))

        return _ContributionResult(q, dq_dt)

    def relax(self, tau: NDArray, d_tau: NDArray):
        states, model, data = self._states, self._model, self._dataset.data
        min_alpha = 1.0
        for r, row in enumerate(data):
            param = self._row_converter(row)
            alpha = self._relax_point(states[r], param, tau, d_tau)
            if alpha < min_alpha:
                min_alpha = alpha
        return min_alpha

    def _solve_point(self,
                     state: NDArray,
                     param: Sequence[float],
                     tau: Sequence[float]) -> _DataPointResult:
        """Solve a point, return dr_dx. state is updated """
        model, config = self._model, self._config
        residual_names = model.vector_res_names(NHKeys.RESIDUALS)
        bound_names = model.vector_res_names(NHKeys.BOUNDS)
        config.linear_solver.reset()
        count_down, max_err = 2, 1.0
        dr_dx = None
        for iteration in range(config.max_iter_inner):
            # evaluate system
            r, dr_dx = self._funcs.f_r(state, param, tau)
            r = squeeze(array(r))
            dr_dx = csr_array(dr_dx)

            if not_finite(r, residual_names):
                raise ValueError("No convergence")

            max_err, max_res_name = assess_residuals(r, residual_names)
            if max_err < 1:
                count_down -= 1
            else:
                count_down = 2
            if not count_down:
                break

            # calculate update and find relaxation factor
            dx = -config.linear_solver.solve(dr_dx, r)
            b, a = self._funcs.f_bx(state, param, tau, dx)
            alpha, min_alpha_name = relax(b, a, bound_names, config.gamma)
            if alpha < config.wall:
                msg = f"Relaxation factor is below {config.wall}, " \
                      "no solution found"
                raise ValueError(msg)

            # apply update
            state = state + alpha * dx
        else:
            if max_err > 1:  # accept solution if it was in count-down
                msg = f"No convergence after {config.max_iter_inner} iterations"
                raise ValueError(msg)

        return _DataPointResult(
            x=state,
            dr_dx=dr_dx
        )

    def _relax_point(self,
                     state: NDArray,
                     param: Sequence[float],
                     tau: Sequence[float],
                     d_tau: Sequence[float]) -> float:
        bound_names = self._model.vector_res_names(NHKeys.BOUNDS)
        b, a = self._funcs.f_bt(state, param, tau, d_tau)
        return relax(b, a, bound_names, self._config.gamma)[0]


class ThermoFitSolver:
    def __init__(self, models: MutMap[NumericHandler],
                 thermo_source: AbstractThermoSource,
                 config: ThermoFitSolverConfig | None = None,
                 **options: Any):
        self._config = (config or ThermoFitSolverConfig()).update(**options)
        self._thermo_source = thermo_source
        self._models = models
        for n, model in models.items():
            try:
                check_model_square(model)
            except NonSquareSystem as err:
                raise NonSquareSystem(
                    variables=err.variables,
                    equations=err.equations,
                    name=f"matrix of model {n}"
                ) from err

    def set_options(self, config: ThermoFitSolverConfig | None = None,
                    **options: Any):
        """Overwrite configuration for subsequent solver runs

       :param config: Options for the solver as defined in
          :class:`~simu.core.solver.thermofit.config.ThermoFitSolverConfig`.
       :param options: overwriting individual configurations directly
       """
        self._config = (config or self._config).update(**options)

    def solve(self, thermo_fit_definition: Map[Any],
              config: ThermoFitSolverConfig | None = None,
              **options: Any) -> ThermoFitReport:
        config = (config or self._config).update(**options)
        definition = self._parse_definition(thermo_fit_definition)
        wrappers = {
            n: ThermoFitContributionWrapper(
                self._models[c.model_id], c,
                definition.datasets[c.dataset],
                definition.parameters,
                config
            ) for n, c in definition.contributions.items()
        }
        tau = _extract_default_values(definition.parameters)
        for iteration in range(config.max_iter_outer):
            sub_results = [w.solve(tau) for w in wrappers.values()]
            jac = vstack([j for s in sub_results for j in s.dq_dt])
            q = concatenate([q_i for s in sub_results for q_i in s.q])
            d_tau, *_ = lstsq(jac, -q)

            alpha = min(w.relax(tau, d_tau) for w in wrappers.values())
            if alpha < config.wall:
                raise ValueError(
                    f"Relaxation factor is below {config.wall} in outer loop; "
                    "no solution found"
                )
            tau += alpha * d_tau

            q_norm = norm(q)
            criterion = abs(q @ jac) / (norm(jac, axis=0) * q_norm + 1e-30)
            # print(iteration, tau, alpha, q_norm, criterion)

            if max(criterion) < config.epsilon or q_norm < config.epsilon_q:
                break
        else:
            raise ValueError("No convergence in outer loop after "
                             f"{config.max_iter_outer} iterations")

        parameters = _generate_parameter_struct(tau, definition.parameters)
        # collect tau
        return ThermoFitReport(
            parameters=parameters
        )


    def _parse_definition(self, definition: Map[Any]) -> ThermoFitDefinition:
        models = {n: ModelContext(m) for n, m in self._models.items()}
        context = ThermoFitValidationContext(models, self._thermo_source)
        return ThermoFitDefinition.model_validate(definition, context=context)

    def _define_thermo_parameter(self, name: str, parameter: ThermoFitParameter
                                 ) -> _ParameterSymbol:
        default_value = parameter.default
        if default_value is None:
            default_value = self._thermo_source[parameter.path]
        return _ParameterSymbol(
            name=name,
            default_value=default_value
        )

def _prepare_functions(model: NumericHandler, cont: ThermoFitContribution,
                       tau_def: Map[ThermoFitParameter]
                      ) -> _FunctionCollection:
    args = model.arguments
    num_states = args[NHKeys.VECTORS][NHKeys.STATES].shape[0]

    # define symbols for function arguments
    x = SX.sym("x", num_states)
    d_x = SX.sym("d_x", num_states)
    t = SX.sym("t", len(tau_def))
    d_t = SX.sym("d_tau", len(tau_def))
    p = SX.sym("p", len(cont.data_to_model))

    # replace state
    args[NHKeys.VECTORS][NHKeys.STATES] = Quantity(x)
    # replace thermo parameters from t in arg
    for t_i, def_i in zip(t.nonzeros(), tau_def.values()):
        path = [def_i.store_name, *def_i.path]
        symbol = Quantity(t_i, def_i.default.units)
        _replace_qty(args[NHKeys.THERMO_PARAMS], symbol, path)
    # replace model parameters from p in arg
    for p_i, def_i in zip(p.nonzeros(), cont.data_to_model.values()):
        symbol = Quantity(p_i, def_i.uom)
        _replace_qty(args[NHKeys.MODEL_PARAMS], symbol, def_i.path)

    # evaluate model symbolically
    res = model.function(args, squeeze_results=False)

    # extract r, b
    vectors = res[NHKeys.VECTORS]
    r, b = vectors[NHKeys.RESIDUALS].m, vectors[NHKeys.BOUNDS].m

    # extract q
    model_props = res[NHKeys.MODEL_PROPS]
    q = vertcat(*[_extract_qty(model_props, path).to("").m
                  for path in cont.penalties])
    # apply weight of entire contribution
    q *= cont.weight

    # create Jacobian matrices
    r_x, r_t = jacobian(r, x), jacobian(r, t)
    q_x, q_t = jacobian(q, x), jacobian(q, t)

    # create functions
    return _FunctionCollection(
        f_r=Function("f_r", [x, p, t], [r, r_x]),
        f_bx=Function("f_bx", [x, p, t, d_x], [b, -b / jtimes(b, x, d_x)]),
        f_bt=Function("f_bt", [x, p, t, d_t], [b, -b / jtimes(b, t, d_t)]),
        f_q=Function("f_q", [x, p, t], [q, q_x, q_t, r_t])
    )


def _extract_qty(results: NestedMap[Quantity], path: Sequence[str]) -> Quantity:
    for p in path:
        results = results[p]
    return results


def _replace_qty(arguments: NestedMutMap[Quantity],
                 item: Quantity, path: Sequence[str]):
    """Replace an item"""
    prev = None
    for p in path:
        if not p in arguments:
            return  # Thermo-parameter is not in model, skip
        prev, arguments = arguments, arguments[p]
    prev[path[-1]] = item


def _extract_default_values(parameters: Map[ThermoFitParameter]) -> NDArray:
    return array([p.default.magnitude for p in parameters.values()])


def _generate_parameter_struct(
        tau: Sequence[float],
        parameters: Map[ThermoFitParameter]) -> NestedMutMap[Quantity]:
    result = {}
    for parameter, tau_i in zip(parameters.values(), tau):
        res = result
        for p in parameter.path[:-1]:
            if p not in res:
                res[p] = {}
            res = res[p]
        res[parameter.path[-1]] = Quantity(tau_i, parameter.default.units)
    return result

