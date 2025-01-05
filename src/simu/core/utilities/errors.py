"""This module defines exception types"""

from pint.errors import DimensionalityError, UndefinedUnitError


class DataFlowError(RuntimeError):
    """An error used when the data flow of objects is not configured correctly.
    """

class IterativeProcessFailed(RuntimeError):
    """Base-class for exceptions in case an iterative process failed."""

class IterativeProcessInterrupted(IterativeProcessFailed):
    """Exception raised when an iterative process, such as a numerical solver
     has been interrupted by user intervention, normally callback functions."""

class NonSquareSystem(ValueError):
    """Exception raised if a system that is to be square is not."""
    def __init__(self, variables: int, equations: int,
                 name: str = "system matrix"):
        excess = "equations" if equations > variables else "variables"
        delta = abs(variables - equations)
        msg = f"Non-square {name}: {variables} variables vs. " \
              f"{equations} equations; {delta} too many {excess}."
        super().__init__(msg)
