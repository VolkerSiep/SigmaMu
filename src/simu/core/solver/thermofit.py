from __future__ import annotations
from typing import TYPE_CHECKING
from dataclasses import dataclass

if TYPE_CHECKING:
    from collections.abc import Sequence, Iterable
    from simu.core.utilities.types import Map, MutMap
    from simu import NumericHandler


@dataclass
class DataSet:
    columns: Sequence[str]
    data: Iterable[Sequence[float]]


@dataclass
class ThermoFitContribution:
    data: DataSet
    model: NumericHandler
    data_model_map: Map[str]  # maps names in data set to model parameters
    penalties: Sequence[str]  # model properties representing penalties
    properties: Sequence[str]  # model properties interesting for evaluation


@dataclass
class ThermoFitEvaluation:
    # map properties from contribution, both data set and properties


class ThermoParameterFit:
    def __init__(self):
        self._contributions: MutMap[ThermoFitContribution] = {}
