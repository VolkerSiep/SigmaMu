"""This module defines exception types"""

from pint.errors import DimensionalityError, UndefinedUnitError

class DataFlowError(RuntimeError):
    """An error used when the data flow of objects is not configured correctly.
    """
