from __future__ import annotations
from typing import Self, Protocol, Optional, Any
from io import Writer
from collections.abc import Sequence
from dataclasses import dataclass
import sys

from pint import DimensionalityError
from pint.registry import Quantity as QtyType
from pydantic import (
    BaseModel, ConfigDict, Field, ValidationInfo,
    model_validator, field_validator)

from simu import Quantity, AbstractThermoSource
from simu.core.utilities.quantity import UnitRegistry
from simu.core.utilities.types import Map, LinearSolver
from ..linear import NumpySolver


class ThermoFitSolverConfig(BaseModel):
    """Configuration settings for the ThermoFitSolver."""

    max_iter_inner: int = Field(default=30, ge=1)
    """The maximum number of iterations (default 30) for solving the
    sub-models for each data point.

    .. note::

      Normally, 30 iterations should be sufficient. In other words, if the
      model is not converged after 30 iterations, chances are quite low
      that it still will converge at all. The advice would be to try to
      improve the starting values and to investigate whether the model is
      properly posed.
    """

    max_iter_outer: int = Field(default=30, ge=1)
    """The maximum number of iterations (default 30) for solving the
    outer iteration on finding the optimal parameter values."""

    output: Writer[str] | None = Field(default_factory=lambda: sys.stdout)
    """The stream to direct the solver output to, by default ``sys.stdout``.
    ``None`` suppresses output. 

    The stream can be any object that supports a ``write`` method that consumes
    a string argument.
    """

    gamma: float = Field(default=0.9, gt=0.0, lt=1.0)
    r""":math:`\gamma` (default 0.9) is the
    fraction of the step-length applied by the solver before hitting the
    domain boundary. Normally, changing the value is not required.
    Generally, a lower value makes the model more robust against
    non-linear domain boundaries (and thus linearisation errors causing
    the state to exit the domain). A higher value yields slightly faster
    convergence, if the solution is in comparison with the initial values
    very close to the domain boundary.
    """

    wall: float = Field(default=1e-20, ge=0.0, lt=0.01)
    r"""Either if there is no solution within the domain of the
    model (for instance: The material balance forces some of the species
    flows in a stream to be negative), or if the solver for other reasons
    is forced to try to leave the model domain, the state will move closer
    and closer to the domain boundary and not revert. At some point,
    :math:`\gamma` becomes ridiculously small, and we need to give up.
    This threshold value is defined by ``wall`` (default ``1e-20``).
    """

    linear_solver: LinearSolver = \
        Field(default_factory=NumpySolver)
    r"""An option to provide any other linear solver for solving the Newton-type
    updates for the inner solving of the sub models for each data point.
    
    As the process models of this type are typically small (say, less than 100
    variables), and not in particular sparse, the default solver is the
    standard dense ``numpy.linalg.solve`` version.
    """

    epsilon: float = Field(default=1e-8, gt=0, lt=1)
    r"""The convergence criterion for the parameter fit, expressed as the
    orthogonal distance of the penalty vector to the sensitivity direction of
    each parameter. For the augmented linearized system
    
    .. math:: J^\mathrm{T}\cdot J\cdot \Delta \tau = -J^\mathrm{T}\,q
    
    the stationary condition :math:`J^\mathrm{T}\,q = 0` is interpreted as the
    potential to improve the solution per parameter :math:`\tau_\alpha` as the
    angle between :math:`\sum_i J_{\alpha i}\,e_i` and :math:`\sum_i q_i e_i`:
    
    .. math::
    
        \max_\alpha \frac{|J^\mathrm{T}\,q|}
          {||J^\mathrm{T}_\alpha||\,||q|| + \varepsilon_m} < \varepsilon
        
    Here, :math:`\varepsilon_m = 10^{-30}` is an artificial parameter to prevent
    division by zero, while :math:`\varepsilon` is the true tolerance parameter.  
    """

    epsilon_q: float = Field(default=0, ge=0, lt=1)
    r"""Normally, :attr:`epsilon` is sufficient to detect convergence, but in
    degenerate cases, such as when the number of data points is equal to the
    number of parameters to fit, :math:`q` approaches very small values and is
    impacted by the remaining residual of the inner model convergence. The
    orthogonality criterion is then less reliable than simply demanding a small
    value of the objective: :math:`||q|| < \varepsilon_q`.
     
     The default value, :math:`\varepsilon_q = 0`, disables this criterion
     unless the objective norm is truly zero.
     """

    model_config = ConfigDict(arbitrary_types_allowed=True, extra="forbid")

    def update(self, **options) -> ThermoFitSolverConfig:
        data = self.model_dump() | options
        if "output" not in options:
            # reverse undesired irreversible serialization of IO stream
            data["output"] = self.output
        return ThermoFitSolverConfig.model_validate(data)


class ThermoFitEvaluationConfig(BaseModel):
    max_iter: int = Field(default=30, ge=1)
    """The maximum number of iterations (default 30) for solving the
    sub-models for each data point.

    .. note::

      Normally, 30 iterations should be sufficient. In other words, if the
      model is not converged after 30 iterations, chances are quite low
      that it still will converge at all.
    """

    gamma: float = Field(default=0.9, gt=0.0, lt=1.0)
    r""":math:`\gamma` (default 0.9) is the
    fraction of the step-length applied by the solver before hitting the
    domain boundary. Normally, changing the value is not required.
    Generally, a lower value makes the model more robust against
    non-linear domain boundaries (and thus linearisation errors causing
    the state to exit the domain). A higher value yields slightly faster
    convergence, if the solution is in comparison with the initial values
    very close to the domain boundary.
    """

    wall: float = Field(default=1e-20, ge=0.0, lt=0.01)
    r"""Either if there is no solution within the domain of the
    model (for instance: The material balance forces some of the species
    flows in a stream to be negative), or if the solver for other reasons
    is forced to try to leave the model domain, the state will move closer
    and closer to the domain boundary and not revert. At some point,
    :math:`\gamma` becomes ridiculously small, and we need to give up.
    This threshold value is defined by ``wall`` (default ``1e-20``).
    """

    linear_solver: LinearSolver = \
        Field(default_factory=NumpySolver)
    r"""An option to provide any other linear solver for solving the Newton-type
    updates for the inner solving of the sub models for each data point.

    As the process models of this type are typically small (say, less than 100
    variables), and not in particular sparse, the default solver is the
    standard dense ``numpy.linalg.solve`` version.
    """

    def update(self, **options) -> ThermoFitEvaluationConfig:
        data = self.model_dump() | options
        return ThermoFitEvaluationConfig.model_validate(data)

    model_config = ConfigDict(arbitrary_types_allowed=True, extra="forbid")


class ThermoFitModelContext(Protocol):
    def parameter_unit(self, path: Sequence[str]) -> str:
        ...

    def property_unit(self, path: Sequence[str]) -> str:
        ...


@dataclass
class ThermoFitValidationContext:
    model_contexts: Map[ThermoFitModelContext]
    thermo_source: AbstractThermoSource | None

def get_context(info: ValidationInfo) -> ThermoFitValidationContext:
    """Help type-analyzer to know regarding the context type for data fit"""
    context = info.context
    assert isinstance(context, ThermoFitValidationContext)
    return context


class DataSet(BaseModel):
    """Represents a collection of experimental data points.

    Printing renders them in restructured text tables:

    >>> d = DataSet(
    ...     columns=["temperature", "pressure", "mole fraction"],
    ...     uom=["K", "bar", "%"],
    ...     data=[[300, 2, 3], [400, 3, 10], [500, 4, 80]],
    ...     source="Invented VLE data")
    >>> print(d)
    --------------- -------------- -----------------
    temperature [K] pressure [bar] mole fraction [%]
    --------------- -------------- -----------------
                300              2                 3
                400              3                10
                500              4                80
    --------------- -------------- -----------------
    """

    columns: Sequence[str]
    """The names of the columns in the dataset."""

    uom: Sequence[str]
    """The units of measurement for each column."""

    data: Sequence[Sequence[float]]
    """The experimental data values, organized as a sequence of rows."""

    source: str = Field(default="Unknown")
    """The source or origin of the dataset."""

    def __str__(self):
        if not self.data:
            return f"<empty dataset>"
        header = [f"{h} [{u}]" for h, u in zip(self.columns, self.uom)]
        lengths = [len(h) for h in header]
        sep = " ".join(["-" * l for l in lengths])
        header = " ".join(header)
        data = "\n".join([" ".join([f"{x:{l}.{l - 4}g}"
                                    for x, l in zip(row, lengths)])
                          for row in self.data])
        return "\n".join([sep, header, sep, data, sep])

    model_config = ConfigDict(extra='forbid')

    def __getitem__(self, column: str) -> Sequence[Quantity]:
        """Return the requested column as a sequence of Quantity objects"""
        idx = self.columns.index(column)
        uom = self.uom[idx]
        return [Quantity(row[idx], uom) for row in self.data]

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


class ModelParameter(BaseModel):
    """Identifies a parameter in the model."""

    path: Sequence[str]
    """The path to the parameter within the model structure."""

    uom: str = Field(default="UNSET")
    """The unit of measurement for the parameter. This field is filled on
    validation based on a query to the model."""

    model_config = ConfigDict(extra="forbid")


class ThermoFitEntity(BaseModel):
    """Base class for entities involved in a thermodynamic fit, associating a
    data set with a model"""

    dataset_id: str
    """The identifier of the dataset to be used."""

    model_id: str
    """The identifier of the model to be used."""

    data_to_model: Map[ModelParameter]
    """Mapping of dataset column titles to model parameters."""

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
            try:
                unit = model.parameter_unit(param.path)
            except KeyError as e:
                param_path = ".".join(param.path)
                msg = (f"Parameter '{param_path}' not defined in "
                       f"model '{self.model_id}'")
                raise ValueError(msg) from e
            param.uom = unit
        return self

    def _model_context(self, info: ValidationInfo) -> ThermoFitModelContext:
        context = get_context(info)
        return context.model_contexts[self.model_id]


class ThermoFitContribution(ThermoFitEntity):
    """Defines a contribution to the objective function for the fit."""

    penalties: Sequence[Sequence[str]]
    """The properties of the model to be used as penalties as a list of paths. 
    A path is a list of strings, identifying the penalty property in the
    hierarchical context of the model. Penalty properties must be dimensionless.
    """

    weight: float = Field(default=1.0, ge=0)
    """The weight applied to this contribution as a non-negative value."""

    @model_validator(mode="after")
    def validate_penalties(self, info: ValidationInfo) -> Self:
        model_id = self.model_id
        model = self._model_context(info)

        for penalty in self.penalties:
            # validate penalty existence
            name = ".".join(penalty)
            try:
                unit = model.property_unit(penalty)
            except KeyError as e:
                msg = f"Property '{name}' not defined in model '{model_id}'"
                raise ValueError(msg) from e
            # validate whether penalties are dimensionless
            if not _Unit(unit).dimensionless:
                msg = (f"Penalty property '{name}` in model "
                       f"'{model_id}' is not dimensionless: '{unit}'")
                raise ValueError(msg)
        return self


class ThermoFitParameter(BaseModel):
    """Defines a thermodynamic parameter to be fitted."""

    path: Sequence[str]
    """The path to the parameter within the thermo source."""

    default: QtyType | None = Field(default=None)
    """The default value for the parameter."""

    lower: QtyType | None = Field(default=None)
    """The lower bound for the parameter."""

    upper: QtyType | None = Field(default=None)
    """The upper bound for the parameter."""

    store_name: str = Field(default="default")
    """The name of the :class:`~simu.ThermoParameterStore` containing the
    parameter."""

    model_config = ConfigDict(extra='forbid', arbitrary_types_allowed=True)

    @field_validator("default", "lower", "upper", mode="before")
    def convert_values(cls, value: str | None,
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
                assert a is not None and b is not None
                msg = f"Incompatible units: '{a:P~}' vs. '{b:P~}'"
                raise ValueError(msg) from e

        l, d, u = self.lower, self.default, self.upper
        err = ""
        if not in_seq(l, u):
            err = f"Upper bound '{u:~P}' less than lower bound '{l:~P}'"
        if not in_seq(l, d):
            assert d is not None and l is not None
            err = f"Default value '{d:~P}' less than lower bound '{l:~P}'"
        if not in_seq(d, u):
            assert d is not None and u is not None
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
        if self.default is None:
            self.default = parameter
        return self


class DataSetQuote(BaseModel):
    """Defines a transfer of a data set property to the evaluation result"""

    name: str
    """The name of the column in the data set"""

    uom: str
    """The target unit of measurement"""


class ThermoFitProperty(BaseModel):
    """Defines a property to be evaluated."""

    path: Sequence[str]
    """The path to the property within the model structure.
    Note that only model properties can be addressed. The model therefore must
    export any relevant thermodynamic (material) properties additionally as a
    model property, e.g. ``self.pr["T_calc"] = my_material["T"]``.
    """

    uom: str
    """The unit of measurement for the property, in which its values will
    be reported in the results."""

    model_config = ConfigDict(extra='forbid')


class ThermoFitEvaluation(ThermoFitEntity):
    """Defines a specific evaluation task, mapping a dataset to a model and
    specifying which properties to evaluate."""

    properties: Map[ThermoFitProperty]
    """The properties to be evaluated. The keys of the provided mapping will
    become the column headers in the resulting data set."""

    quote_data: Map[DataSetQuote] = Field(default_factory=dict)
    """Columns of the data set to be quoted in the evaluation table, supporting
    unit conversion."""

    @model_validator(mode="after")
    def validate_properties(self, info: ValidationInfo) -> Self:
        model = self._model_context(info)
        for prop in self.properties.values():
            path, uom = prop.path, prop.uom
            name = ".".join(path)
            # Does property exist in model?
            try:
                model_unit = model.property_unit(path)
            except KeyError as e:
                msg = f"Property '{name}' not in model '{self.model_id}'"
                raise ValueError(msg) from e

            # Is the unit string a valid unit of measurement?
            # Are the units compatible?
            if not _are_units_compatible(uom, model_unit):
                msg = (f"Property '{name}' has incompatible unit `{uom}`"
                       f"to mapped model property (`{model_unit}`)")
                raise ValueError(msg)
        return self

    model_config = ConfigDict(extra='forbid')


class ThermoFitDefinition(BaseModel):
    """Defines the complete thermodynamic fit problem."""

    datasets: Map[DataSet]
    """The datasets used in the fit."""

    contributions: Map[ThermoFitContribution]
    """The contributions to the objective function."""

    parameters: Map[ThermoFitParameter]
    """The parameters to be fitted."""

    evaluations: Map[ThermoFitEvaluation] | None = Field(default=None)
    """Ignored in the context of data fit, but accepted, so the same data
    structure can be used for data fit and evaluation."""

    @field_validator("evaluations", mode="before")
    def ignore_evaluation_fields(cls, value: Any) -> None:
        return None

    @model_validator(mode="after")
    def validate_configuration(self, info: ValidationInfo) -> Self:
        context = get_context(info).model_contexts
        for contribution in self.contributions.values():
            _validate_entity(contribution, self.datasets, context)
        return self

    model_config = ConfigDict(extra='forbid')


class ThermoFitEvaluationDefinition(BaseModel):
    """Defines the configuration for a thermodynamic evaluation run.

    This structure mirrors the fit definition but focuses on evaluating model
    properties against datasets rather than fitting parameters. It validates
    that all referenced datasets and models are correctly defined and
    compatible.
    """
    datasets: Map[DataSet]
    """The datasets used in the evaluation."""

    evaluations: Map[ThermoFitEvaluation]
    """The evaluations to be performed."""

    contributions: Map[ThermoFitContribution] | None = Field(default=None)
    """Ignored in the context of evaluation, but accepted, so the same data
    structure can be used for data fit and evaluation."""

    parameters: Map[ThermoFitParameter] | None = Field(default=None)
    """Ignored in the context of evaluation, but accepted, so the same data
    structure can be used for data fit and evaluation."""

    @field_validator("contributions", "parameters", mode="before")
    def ignore_datafit_fields(cls, value: Any) -> None:
        return None


    @model_validator(mode="after")
    def validate_configuration(self, info: ValidationInfo) -> Self:
        context = get_context(info).model_contexts
        for contribution in self.evaluations.values():
            _validate_entity(contribution, self.datasets, context)
        return self

    @model_validator(mode="after")
    def validate_quotes(self, info: ValidationInfo) -> Self:
        for evaluation in self.evaluations.values():
            dataset = self.datasets[evaluation.dataset_id]
            for quote in evaluation.quote_data.values():
                idx = dataset.columns.index(quote.name)
                ds_unit = dataset.uom[idx]
                if not _are_units_compatible(quote.uom, ds_unit):
                    msg = f"Incompatible units '{quote.uom}' vs. '{ds_unit}'"
                    raise ValueError(msg)
        return self

    model_config = ConfigDict(extra='forbid')


def _validate_entity(entity: ThermoFitEntity,
                     datasets: Map[DataSet],
                     context: Map[ThermoFitModelContext]):
    try:
        dataset = datasets[entity.dataset_id]
    except KeyError as e:
        msg = f"Dataset '{entity.dataset_id}' not defined"
        raise ValueError(msg) from e

    # validate mapping
    model = context[entity.model_id]
    for key, target in entity.data_to_model.items():
        # is column defined in dataset
        if key not in dataset.columns:
            msg = f"Column '{key}' not defined in dataset"
            raise ValueError(msg)
        uom = dataset.uom[dataset.columns.index(key)]

        # is target defined in model (checked before?)
        try:
            uom_model = model.parameter_unit(target.path)
        except KeyError as e:
            name = ".".join(target.path)
            msg = f"Target '{name}' not a model parameter"
            raise ValueError(msg) from e

        # are units compatible?
        if not _are_units_compatible(uom, uom_model):
            msg = f"Incompatible units '{uom}' vs. '{uom_model}'"
            raise ValueError(msg)


def _are_units_compatible(first: str, second: str) -> bool:
    try:
        d1 = _Unit(first).dimensionality
    except Exception as e:
        raise ValueError(f"Invalid unit '{first}'") from e
    try:
        d2 = _Unit(second).dimensionality
    except Exception as e:
        raise ValueError(f"Invalid unit '{second}'") from e
    return d1 == d2

_Unit = UnitRegistry.Unit
