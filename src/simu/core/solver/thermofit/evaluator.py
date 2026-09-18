from typing import Any

from simu import (
    NumericHandler, SimulationSolver, Quantity, PropertyFilter, NHKeys)
from simu.core.utilities.types import MutMap, Map, NestedMap, NestedMutMap
from simu.core.utilities.errors import NonSquareSystem

from .config import (
    DataSet, ThermoFitEvaluationConfig, ThermoFitValidationContext,
    ThermoFitEvaluation, ThermoFitEvaluationDefinition)
from .report import ThermoEvaluationReport
from ..common import ModelContext, check_model_square, replace_qty, extract_qty


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
        param_uom = [d.uom for d in evaluation.data_to_model.values()]
        self._evaluation = evaluation

    def set_thermo_parameters(self, parameters: NestedMap[Quantity]):
        target = self._solver.model_parameters[NHKeys.THERMO_PARAMS]
        _overwrite_nodes(target, parameters)

    def solve(self, name: str) -> ThermoEvaluationReport:
        dataset = self._dataset
        evaluation = self._evaluation
        num_failed = 0
        model_parameters = self._solver.model_parameters[NHKeys.MODEL_PARAMS]
        parameter_paths = [p.path for p in evaluation.data_to_model.values()]

        # prepend quoted data
        quotes = evaluation.quote_data
        columns = list(quotes.keys())
        uom_quotes = [q.uom for q in quotes.values()]
        name_quotes = [q.name for q in quotes.values()]
        quote_idx = [dataset.columns.index(n) for n in name_quotes]
        uom_orig = [dataset.uom[i] for i in quote_idx]

        # append result data
        columns += evaluation.properties.keys()
        uom = uom_quotes + [prop.uom for prop in evaluation.properties.values()]


        num_properties = len(columns)
        results = []

        for r, row in enumerate(dataset.data):  # TODO: parallelize this loop
            # set parameters to model
            for magnitude, uom_i, path in zip(row, dataset.uom, parameter_paths):
                replace_qty(model_parameters, Quantity(magnitude, uom_i), path)
            try:
                result = self._solver.solve()
            except ValueError:
                results.append([float("nan")] * num_properties)
                num_failed += 1
                continue
            model_props = result.properties[NHKeys.MODEL_PROPS]
            result_row = [
                extract_qty(model_props, prop.path).to(prop.uom).magnitude
                for prop in evaluation.properties.values()
            ]
            quote = [
                Quantity(row[idx], u_o).to(u_q).magnitude
                for idx, u_q, u_o in zip(quote_idx, uom_quotes, uom_orig)
            ]
            results.append(quote + result_row)



        return ThermoEvaluationReport(
            results=DataSet(
                columns=columns,
                uom=uom,
                data=results,
                source=f"Evaluation '{name}'"
            ),
            num_failed=num_failed
        )


class ThermoFitEvaluator:
    """Orchestrates the evaluation of thermodynamic models against experimental datasets.

    This class manages a collection of models and executes evaluations based on
    provided configurations. It requires that all models are square systems and
    well-formed."""
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
        """Executes the evaluation defined by the provided configuration.

        :param thermo_evaluation_definition: A dictionary or mapping defining
          the datasets and evaluations to be performed.
        :param parameters: Optional thermodynamic parameters to be applied to
          the models before evaluation. These can be a result from a previous
          thermo fit run:
          :attr:`~simu.core.solver.thermofit.report.ThermoFitReport.final_parameters`.
        :param options: Additional configuration overrides for this specific
          run.
        :return: A mapping of evaluation names to their respective reports.
        """
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
            reports[name] = evaluator.solve(name)
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


def _overwrite_nodes(
        target: NestedMutMap[Quantity] | Quantity,
        parameters: NestedMap[Quantity] | Quantity
) -> Quantity | None:
    """Overwrite nodes from parameters in target structure, but do not create
    new ones. Throw error if structure is incompatible (i.e. if one of the
    structures exposes a leaf while the other holds a sub-structure under the
    same key. Also throw error if physical dimensions mismatch.

    >>> t = {"a": {"b": Quantity(3, "m"), "c": Quantity(4, "s")}}
    >>> p = {"a": {"b": Quantity(5, "cm"), "d": Quantity(4, "K")}}
    >>> _overwrite_nodes(t, p)
    >>> print(t)
    {'a': {'b': <Quantity(5, 'centimeter')>, 'c': <Quantity(4, 'second')>}}
    """
    if isinstance(target, Quantity) and isinstance(parameters, Quantity):
        if not target.check(parameters):
            msg = "Incompatible physical dimension while overwriting nodes"
            raise ValueError(msg)
        return parameters

    if isinstance(target, Quantity) or isinstance(parameters, Quantity):
        raise ValueError("Incompatible structure overwriting nodes")

    for key in set(target.keys()) & set(parameters.keys()):
        target_value = _overwrite_nodes(target[key], parameters[key])
        if target_value is not None:
            target[key] = target_value
    return None
