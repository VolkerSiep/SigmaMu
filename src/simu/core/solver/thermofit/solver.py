from typing import Any
from collections.abc import Sequence
from dataclasses import dataclass
from casadi import SX, jacobian, jtimes, Function
from simu import NumericHandler, NHKeys, AbstractThermoSource
from simu.core.utilities.types import Map, MutMap, NestedMap

from .config import (
    ThermoFitDefinition, ThermoFitSolverConfig, ThermoFitValidationContext,
    ThermoFitContribution
)

@dataclass
class _FunctionCollection:
    f_r: Function  # t, x, p -> r, r_x
    f_bx: Function  # t, x, dx, p -> a, b
    f_bt: Function  # t, x, dt, p -> a, b
    f_q: Function  # t, x, p -> q, q_x, q_t, r_x, r_t


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

    def solve(self, definition: Map[Any],
              config: ThermoFitSolverConfig | None = None,
              **options: Any):
        config = (config or self._config).update(**options)
        setup = self._parse_definition(definition)
        funcs = {n: self._prepare_functions(c)
                 for n, c in setup.contributions.items()}


        # for each contribution, collect the model and create the required functions


        # need to identify thermodynamic parameters in arguments

    def _prepare_functions(self,
                           cont: ThermoFitContribution) -> _FunctionCollection:
        model = self._models[cont.model_id]




    def _parse_definition(self, definition: Map[Any]) -> ThermoFitDefinition:
        models = {n: ModelContext(m) for n, m in self._models.items()}
        context = ThermoFitValidationContext(models, self._thermo_source)
        return ThermoFitDefinition.model_validate(definition, context=context)
