from collections.abc import Sequence
from dataclasses import dataclass
from typing import Self, Protocol, Optional

from pint import DimensionalityError
from pint.registry import Quantity as QtyType
from pydantic import (
    BaseModel, ConfigDict, Field, ValidationInfo,
    model_validator, field_validator)

from simu import Quantity, AbstractThermoSource
from simu.core.utilities.quantity import UnitRegistry
from simu.core.utilities.types import Map

_Unit = UnitRegistry.Unit

def are_units_compatible(first: str, second: str) -> bool:
    try:
        d1 = _Unit(first).dimensionality
    except Exception as e:
        raise ValueError(f"Invalid unit '{first}'") from e
    try:
        d2 = _Unit(second).dimensionality
    except Exception as e:
        raise ValueError(f"Invalid unit '{second}'") from e
    return d1 == d2


class ThermoFitModelContext(Protocol):
    @property
    def parameters(self) -> Map[str]:
        ...
    @property
    def properties(self) -> Map[str]:
        ...


@dataclass
class ThermoFitValidationContext:
    model_contexts: Map[ThermoFitModelContext]
    thermo_source: AbstractThermoSource

def get_context(info: ValidationInfo) -> ThermoFitValidationContext:
    return info.context


class DataSet(BaseModel):
    columns: Sequence[str]
    uom: Sequence[str]
    data: Sequence[Sequence[float]]
    source: str = Field(default="Unknown")

    model_config = ConfigDict(extra='forbid')

    @field_validator("uom", mode="after")
    def check_uom(cls, value: Sequence[str]) -> Sequence[str]:
        for uom in value:
            try:
                _Unit(uom)
            except Exception as e:
                raise ValueError(f"Invalid unit of measurement '{uom}'") from e
        return value

    @model_validator(mode="after")
    def check_dimensions(self) -> Self:
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

class ThermoFitEntity(BaseModel):
    dataset: str  # to be registered dataset
    model_id: str  # to be registered model
    data_to_model: Map[str]  # keys to be column titles, values model parameters

    model_config = ConfigDict(extra='forbid')

    @field_validator("model_id", mode="after")
    def validate_model_id(cls, value: str, info: ValidationInfo) -> str:
        context = get_context(info)
        if value not in context.model_contexts:
            raise ValueError(f"Model '{value}' is not registered")
        return value

    @model_validator(mode="after")
    def validate_model_parameters(self, info: ValidationInfo) -> Self:
        model = self._model_context(info)

        # validate parameter existence
        for param in self.data_to_model.values():
            if not param in model.parameters:
                msg = (f"Parameter '{param}' not defined in "
                       f"model '{self.model_id}'")
                raise ValueError(msg)
        return self

    def _model_context(self, info: ValidationInfo) -> ThermoFitModelContext:
        context = get_context(info)
        return context.model_contexts[self.model_id]


class ThermoFitContribution(ThermoFitEntity):
    penalties: Sequence[str]  # to be properties of model
    weight: float = Field(default=1.0)

    @model_validator(mode="after")
    def validate_penalties(self, info: ValidationInfo) -> Self:
        model_id = self.model_id
        model = self._model_context(info)

        for penalty in self.penalties:
            # validate penalty existence
            try:
                unit = model.properties[penalty]
            except KeyError as e:
                msg = f"Property '{penalty}' not defined in model '{model_id}'"
                raise ValueError(msg) from e
            # validate whether penalties are dimensionless
            if not _Unit(unit).dimensionless:
                msg = (f"Penalty property '{penalty}` in model "
                       f"'{model_id}' is not dimensionless: '{unit}'")
                raise ValueError(msg)
        return self


class ThermoFitProperty(BaseModel):
    name: str  # must exist in model
    uom: str  # must be consistent with unit from model
    model_config = ConfigDict(extra='forbid')


class ThermoFitEvaluation(ThermoFitEntity):
    properties: Map[ThermoFitProperty]

    @model_validator(mode="after")
    def validate_properties(self, info: ValidationInfo) -> Self:
        properties = self._model_context(info).properties
        for prop in self.properties.values():
            name, uom = prop.name, prop.uom
            # Does property exist in model?
            try:
                model_unit = properties[name]
            except KeyError as e:
                msg = f"Property '{name}' not in model '{self.model_id}'"
                raise ValueError(msg) from e

            # Is the unit string a valid unit of measurement?
            # Are the units compatible?
            if not are_units_compatible(uom, model_unit):
                msg = (f"Property '{name}' has incompatible unit `{uom}`"
                       f"to mapped model property (`{model_unit}`)")
                raise ValueError(msg)
        return self


class ThermoFitParameter(BaseModel):
    path: Sequence[str]
    default: QtyType = Field(default=None)
    lower: QtyType = Field(default=None)
    upper: QtyType = Field(default=None)

    model_config = ConfigDict(extra='forbid', arbitrary_types_allowed=True)

    @field_validator("default", "lower", "upper", mode="before")
    def convert_values(cls, value: str,
                         info: ValidationInfo) -> QtyType | None:
        if value is None:
            return None
        try:
            return Quantity(value)
        except Exception as e:
            raise ValueError(f"Invalid {info.field_name}: {value} - {e}") from e

    @model_validator(mode="after")
    def validate_sequence(self) -> Self:
        def in_seq(a: Optional[QtyType], b: Optional[QtyType]) -> bool:
            try:
                return True if None in (a, b) else a < b
            except DimensionalityError as e:
                msg = f"Incompatible units: '{a:P~}' vs. '{b:P~}'"
                raise ValueError(msg) from e

        l, d, u = self.lower, self.default, self.upper
        err = ""
        if not in_seq(l, u):
            err = f"Upper bound '{u:~P}' less than lower bound '{l:~P}'"
        if not in_seq(l, d):
            err = f"Default value '{d:~P}' less than lower bound '{l:~P}'"
        if not in_seq(d, u):
            err = f"Default value '{d:~P}' more than upper bound '{u:~P}'"
        if err:
            raise ValueError(err)
        return self

    @model_validator(mode="after")
    def validate_vs_source(self, info: ValidationInfo) -> Self:
        # check existence and dimensional compatibility in thermo source
        source = get_context(info).thermo_source
        try:
            parameter = source[self.path]
        except KeyError as e:
            msg = f"'{'.'.join(self.path)}' not found in ThermoSource object"
            raise ValueError(msg) from e
        for value in (self.default, self.lower, self.upper):
            if value is None:
                continue
            if not value.check(parameter.units):
                msg = f"Incompatible unit: {value.units} vs. {parameter.units}"
                raise ValueError(msg)
        return self


class ThermoFitConfiguration(BaseModel):
    datasets: Map[DataSet]
    contributions: Map[ThermoFitContribution]
    evaluations: Map[ThermoFitEvaluation]
    parameters: Map[ThermoFitParameter]

    @model_validator(mode="after")
    def validate_configuration(self, info: ValidationInfo) -> Self:
        context = get_context(info).model_contexts
        self._validate_item(self.contributions, context)
        self._validate_item(self.evaluations, context)
        return self

    def _validate_item(self, mapping: Map[ThermoFitEntity],
                       context: Map[str, ThermoFitModelContext]):
        for entity in mapping.values():
            # does dataset exist
            try:
                dataset = self.datasets[entity.dataset]
            except KeyError as e:
                msg = f"Dataset '{entity.dataset}' not defined"
                raise ValueError(msg) from e

            # validate mapping
            parameters = context[entity.model_id].parameters
            for key, target in entity.data_to_model.items():
                # is column defined in dataset
                if key not in dataset.columns:
                    msg = f"Column '{key}' not defined in dataset"
                    raise ValueError(msg)
                uom = dataset.uom[dataset.columns.index(key)]

                # is target defined in model (checked before?)
                try:
                    uom_model = parameters[target]
                except KeyError as e:
                    msg = f"Target '{target}' not a model parameter"
                    raise ValueError(msg) from e

                # are units compatible?
                if not are_units_compatible(uom, uom_model):
                    msg = f"Incompatible units '{uom}' vs. '{uom_model}'"
                    raise ValueError(msg)
