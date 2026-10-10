from typing import Any, ClassVar, Callable
from pint.registry import Quantity as QtyType
from pydantic_core import core_schema

from simu import Quantity

class PQuantity:
    NAME: ClassVar[str]
    UNIT: ClassVar[str]
    CHECK: ClassVar[Callable[[Any], bool]] = staticmethod(lambda qty: True)

    @classmethod
    def __get_pydantic_core_schema__(cls, source, handler):
        def validate(value: Any) -> QtyType:
            try:
                qty = Quantity(value)
            except Exception as e:
                msg = f"Error parsing {cls.NAME} of value '{value}': {e}"
                raise ValueError(msg)
            if not qty.check(cls.UNIT):
                msg = f"Incompatible unit for {cls.NAME} in value '{value}'"
                raise ValueError(msg)
            if not cls.CHECK(qty.m):
                msg = f"Constraint validation for {cls.NAME} in value '{value}'"
                raise ValueError(msg)
            return qty

        return core_schema.no_info_after_validator_function(
            validate,
            core_schema.any_schema(),
        )

class PTemperature(PQuantity):
    NAME = "Temperature"
    UNIT = "K"
    CHECK = staticmethod(lambda qty: qty > 0)


class PPressure(PQuantity):
    NAME = "Pressure"
    UNIT = "Pa"
    CHECK = staticmethod(lambda qty: qty > 0)


class PAmount(PQuantity):
    NAME = "Amount"
    UNIT = "mol"
    CHECK = staticmethod(lambda qty: qty > 0)
