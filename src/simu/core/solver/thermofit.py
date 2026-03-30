from typing import Self, Any
from collections.abc import Sequence
from pydantic import BaseModel, ConfigDict, model_validator

class DataSet(BaseModel):
    source: str = "Unknown"
    columns: Sequence[str]
    uom: Sequence[str]
    data: Sequence[Sequence[float]]

    model_config = ConfigDict(extra='forbid')

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

