from typing import Any
from simu import NumericHandler, NHKeys
from simu.core.utilities.types import Map, MutMap

from .config import (
    ThermoFitDefinition, ThermoFitSolverConfig, ThermoFitValidationContext)

class _ModelContext:
    def __init__(self, model: NumericHandler):
        # TODO:
        #  - change protocol and context, so that parameters and properties
        #    are nested structures
        #  - re-understand what I tried to do with thermo-parameters.
        #    maybe do not fiddle around with ThermoSources, but allow any path
        #    of thermo parameters that is available in the models
        #  - Should I separate the structures for data fit and evaluation?
        #    They really do not need to be provided at the same time,
        #    and the evaluator can even be an entirely different object.

        self.parameters = model.arguments[NHKeys.MODEL_PARAMS]


class ThermoFitSolver:
    def __init__(self, models: MutMap[NumericHandler],
                 config: ThermoFitSolverConfig | None = None,
                 **options: Any):
        self._config = (config or ThermoFitSolverConfig()).update(**options)
        self._models = models

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
        context = self._create_validation_context()
        ThermoFitDefinition.model_validate(config, context=context)


        # for each data set, collect the model and create the required functions
        # need to identify thermodynamic parameters in arguments


    def _create_validation_context(self) -> ThermoFitValidationContext:
        ...
