from typing import Any

from simu import NumericHandler, SimulationSolver, Quantity, PropertyFilter
from simu.core.utilities.types import MutMap, Map, NestedMap
from simu.core.utilities.errors import NonSquareSystem

from .config import (
    DataSet, ThermoFitEvaluationConfig, ThermoFitValidationContext,
    ThermoFitEvaluation, ThermoFitEvaluationDefinition)
from .report import ThermoEvaluationReport
from ..common import ModelContext, check_model_square, DataRowConverter


class NoThermoPropFilter(PropertyFilter):
    """This filter removes all stream properties, as only model properties
    are addressable by the evaluation."""

    def keep_property(self, name: str, sub_key: str = None) -> bool:
        return False


class ThermoFitSingleEvaluator:
    def __init__(
            self,
            model: NumericHandler,
            dataset: DataSet,
            evaluation: ThermoFitEvaluation,
            config: ThermoFitEvaluationConfig
    ):
        model.set_property_filter(NoThermoPropFilter())
        self._solver = SimulationSolver(
            model,
            max_iter=config.max_iter,
            gamma=config.gamma,
            wall=config.wall,
            output=None,
            retain_solutioon=False
        )
        self._dataset = dataset

    def set_thermo_parameters(self, parameters: NestedMap[Quantity]):
        # first entry is store name - done!
        pass

    def solve(self, ) -> ThermoEvaluationReport:
        #
        pass


class ThermoFitEvaluator:
    def __init__(
            self, models: MutMap[NumericHandler],
            config: ThermoFitEvaluationConfig | None = None,
            **options: Any
    ):
        self._config = (config or ThermoFitEvaluationConfig()).update(**options)
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

    def set_options(
            self, config: ThermoFitEvaluationConfig | None = None,
            **options: Any
    ):
        """Update the solver configuration for subsequent runs.

        :param config: A new configuration object to replace the current one.
        :param options: Individual configuration parameters to override.
        """
        self._config = (config or self._config).update(**options)

    def solve(
            self, thermo_evaluation_definition: Map[Any],
            parameters: NestedMap[Quantity] | None = None,
            **options: Any
    ) -> Map[ThermoEvaluationReport]:
        config = self._config.update(**options)
        definition = self.parse_definition(thermo_evaluation_definition)
        models = self._models
        datasets = definition.datasets
        parameters = parameters or {}
        reports = {}
        for name, evaluation in definition.evaluations.items():
            evaluator = ThermoFitSingleEvaluator(
                models[evaluation.model_id],
                datasets[evaluation.dataset_id],
                evaluation, config
            )
            evaluator.set_thermo_parameters(parameters)
            reports[name] = evaluator.solve()
        return reports

    def parse_definition(
            self, definition: Map[Any]
    ) -> ThermoFitEvaluationDefinition:
        """Parse and validate a raw thermodynamic fit definition.

        :param definition: A dictionary or mapping representing the evaluation
          configuration.
        :return: A validated
          :class:`~simu.core.solver.thermofit.config.ThermoFitEvaluationDefinition`
          object.
        """
        context = ThermoFitValidationContext(
            model_contexts={n: ModelContext(m) for n, m in self._models.items()},
            thermo_source=None
        )
        return ThermoFitEvaluationDefinition.model_validate(
            definition, context=context
        )
