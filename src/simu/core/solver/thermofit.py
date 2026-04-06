from typing import Self
from collections.abc import Sequence
from pydantic import (
    BaseModel, ConfigDict, Field, ValidationInfo,
    model_validator, field_validator)
from pint.registry import Quantity as QtyType

from simu.core.utilities.types import Map
from simu.core.utilities.quantity import UnitRegistry


_U = UnitRegistry.Unit
_Q = UnitRegistry.Quantity

class DataSet(BaseModel):
    columns: Sequence[str]
    uom: Sequence[str]
    data: Sequence[Sequence[float]]
    source: str = Field(default="Unknown")

    model_config = ConfigDict(extra='forbid')

    @field_validator("uom", mode="after")
    @classmethod
    def check_uom(cls, value: Sequence[str]) -> Sequence[str]:
        for uom in value:
            try:
                _U(uom)
            except Exception:
                raise ValueError(f"Invalid unit of measurement '{uom}'")
        return value

    @model_validator(mode="after")
    def check_dimensions(self) -> Self:
        print("checking dimensions")
        l_columns = len(self.columns)
        l_uom = len(self.uom)
        if l_uom != l_columns:
            raise ValueError(f"Number of units ({l_uom}) does not match "
                             f"number of columns ({l_columns})")
        for k, row in enumerate(self.data):
            l_row = len(row)
            if l_row != l_columns:
                raise ValueError(f"Number of values in row {k} ({l_row}) does "
                                 f"not match number of columns ({l_columns})")
        return self


class ThermoFitContribution(BaseModel):
    dataset: str  # to be registered dataset
    model_id: str  # to be registered model
    data_to_model: Map[str]  # keys to be column titles, values model parameters
    penalties: Sequence[str]  # to be properties of model
    weight: float = Field(default=1.0)

    model_config = ConfigDict(extra='forbid')

    @field_validator("model_id", mode="after")
    @classmethod
    def validate_model_id(cls, value: str, info: ValidationInfo) -> str:
        if value not in info.context:
            raise ValueError(f"Model '{value}' is not registered")
        return value

    @model_validator(mode="after")
    def validate_model_parameters(self, info: ValidationInfo) -> Self:
        model_id = self.model_id
        model = info.context[model_id]

        # validate parameter existence
        for param in self.data_to_model.values():
            if not param in model.parameters:
                msg = f"Parameter '{param}' not defined in model '{model_id}'"
                raise ValueError(msg)

        for penalty in self.penalties:
            # validate penalty existence
            try:
                unit = model.properties[penalty]
            except KeyError:
                msg = f"Property '{penalty}' not defined in model '{model_id}'"
                raise ValueError(msg)
            # validate whether penalties are dimensionless
            if not _U(unit).dimensionless:
                msg = (f"Penalty property '{penalty}` in model "
                       f"'{model_id}' is not dimensionless: '{unit}'")
                raise ValueError(msg)
        return self


class ThermoFitProperty(BaseModel):
    name: str  # must exist in model
    uom: str  # must be consistent with unit from model
    model_config = ConfigDict(extra='forbid')

    # todo: check existence of parameter and dimensionality in overall model

class ThermoFitParameter(BaseModel):
    name: str
    default: QtyType = Field(default=None)  # default value (to be parsed as qty
    lower: QtyType = Field(default=None)
    upper: QtyType = Field(default=None)

    model_config = ConfigDict(extra='forbid', arbitrary_types_allowed=True)

    @field_validator("default", "lower", "upper", mode="before")
    @classmethod
    def validate_default(cls, value: str, info: ValidationInfo):
        if value is None:
            return None
        try:
            return _Q(value)
        except Exception as e:
            raise ValueError(f"Invalid {info.field_name}: {value} - {e}")

    # TODO: check existence, sequence and dimensional compatibility in overall model

class ThermoFitEvaluation(BaseModel):
    dataset: str  # to be registered dataset
    model_id: str  # to be registered model
    data_to_model: Map[str]  # keys to be column titles, values model parameters
    properties: Map[ThermoFitProperty]  # to be properties of model

    model_config = ConfigDict(extra='forbid')

    # TODO: check constraints


class ThermoFitConfiguration(BaseModel):
    datasets: Map[DataSet]
    contributions: Map[ThermoFitContribution]
    parameters: Map[ThermoFitParameter]
    evaluations: Map[ThermoFitEvaluation]

    model_config = ConfigDict(extra='forbid')

    # TODO on integration level - in contributions and evaluations
    #  - evaluate existence of data set
    #  - evaluate existence of columns in data set (data_to_model) and unit consistency of parameters


# TODO: document everything (well!)
