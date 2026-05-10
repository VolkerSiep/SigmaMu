from typing import Any
from collections.abc import Sequence
from dataclasses import dataclass
from casadi import SX, jacobian, jtimes, Function, vertcat
from simu import (
    NumericHandler, NHKeys, AbstractThermoSource, Quantity, SymbolQuantity)
from simu.core.utilities.types import Map, MutMap, NestedMap, NestedMutMap

from .config import (
    ThermoFitDefinition, ThermoFitSolverConfig, ThermoFitValidationContext,
    ThermoFitContribution, ThermoFitParameter
)

@dataclass
class _ParameterSymbol:
    name: str
    default_value: Quantity


@dataclass
class _FunctionCollection:
    f_r: Function  # x, p, t -> r, r_x
    f_bx: Function  # x, p, t, dx -> b, a
    f_bt: Function  # x, p, t, dt -> b, a
    f_q: Function  # x, p, t -> q, q_x, q_t, r_x, r_t


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


class ThermoFitSolver:
    def __init__(self, models: MutMap[NumericHandler],
                 thermo_source: AbstractThermoSource,
                 config: ThermoFitSolverConfig | None = None,
                 **options: Any):
        self._config = (config or ThermoFitSolverConfig()).update(**options)
        self._thermo_source = thermo_source
        self._models = models
        # TODO: check models to be square, store sizes,

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
              **options: Any):
        config = (config or self._config).update(**options)
        definition = self._parse_definition(thermo_fit_definition)
        funcs = {n: _prepare_functions(self._models[c.model_id],
                                       c, definition.parameters)
                 for n, c in definition.contributions.items()}


        # for each contribution, collect the model and create the required functions


        # need to identify thermodynamic parameters in arguments





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
        symbol = Quantity(t_i, def_i.default.units)
        _replace_qty(args[NHKeys.THERMO_PARAMS], symbol, def_i.path)
    # replace model parameters from p in arg
    for p_i, def_i in zip(p.nonzeros(), cont.data_to_model.values()):
        symbol = Quantity(p_i, def_i.uom)
        _replace_qty(args[NHKeys.THERMO_PARAMS], symbol, def_i.path)

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
        f_q=Function("f_q", [x, p, t], [q, q_x, q_t, r_x, r_t])
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