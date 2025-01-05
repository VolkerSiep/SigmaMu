"""This module defines exception types"""

from pint.errors import DimensionalityError, UndefinedUnitError


class DataFlowError(RuntimeError):
    """An error used when the data flow of objects is not configured correctly.
    """

class IterativeProcessInterrupted(RuntimeError):
    """Exception raised when an iterative process, such as a numerical solver
     has been interrupted by user intervention, normally callback functions."""