from typing import Any, Callable
from pint.registry import Quantity as QtyType
from pydantic_core import core_schema

from simu import Quantity

class PQuantity:
    NAME: str
    UNIT: str
    CHECK = staticmethod(lambda qty: True)

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

def new_p_quantity(name: str, unit: str,
                   check: Callable[[Any], bool]) -> type[PQuantity]:
    class NewPQuantity(PQuantity):
        NAME = name
        UNIT = unit
        CHECK = staticmethod(check)
    return NewPQuantity


PTemperature = new_p_quantity("Temperature", "K", lambda q: q > 0)
PPressure = new_p_quantity("Pressure", "Pa", lambda q: q > 0)
PAmount = new_p_quantity("Amount", "mol", lambda q: q > 0)
