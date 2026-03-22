from __future__ import annotations
from typing import TYPE_CHECKING
from dataclasses import dataclass

if TYPE_CHECKING:
    from collections.abc import Sequence, Iterable
    from simu.core.utilities.types import Map, MutMap
    from simu import NumericHandler


@dataclass
class DataSet:  # all data is held in data sets
    columns: Sequence[str]
    data: Iterable[Sequence[float]]


@dataclass
class DataSeries:
    contribution_id: str
    model_property: bool
    property_id: str


Evaluation = Map[DataSeries]


@dataclass
class ThermoFitContribution:
    data: DataSet
    fit_model_id: str
    eval_model_id: str
    data_to_fit_model: Map[str]  # maps names in data set to model parameters
    data_to_eval_model: Map[str]
    penalties: Sequence[str]  # model properties representing penalties
    properties: Sequence[str]  # model properties interesting for evaluation
    weight: float

class ThermoParameterFit:
    def __init__(self):
        self._contributions: MutMap[ThermoFitContribution] = {}
        self._evaluations: MutMap[Evaluation] = {}
        self._models: MutMap[NumericHandler] = {}

    def add_contribution(self, name: str,
                         contribution: ThermoFitContribution):
        if name in self._contributions:
            raise ValueError(f"Contribution of name {name} already defined")
        self._contributions[name] = contribution

    def add_evaluation(self, name: str, evaluation: Evaluation):
        if name in self._evaluations:
            raise ValueError(f"Evaluation of name {name} already defined")
        self._evaluations[name] = evaluation

    def add_model(self, name: str, model: NumericHandler):
        if name in self._models:
            raise ValueError(f"Model of name {name} already defined")
        self._models[name] = model
