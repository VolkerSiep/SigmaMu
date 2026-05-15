"""This module defines exception types"""

# external
from dataclasses import dataclass, field
from pint.errors import DimensionalityError, UndefinedUnitError


class DataFlowError(RuntimeError):
    """An error used when the data flow of objects is not configured correctly.
    """

class IterativeProcessFailed(RuntimeError):
    """Base-class for exceptions in case an iterative process failed."""

class IterativeProcessInterrupted(IterativeProcessFailed):
    """Exception raised when an iterative process, such as a numerical solver
     has been interrupted by user intervention, normally callback functions."""

@dataclass
class NonSquareSystem(ValueError):
    """Exception raised if a system that is to be square is not."""
    variables: int
    equations: int
    name: str = field(default="system matrix")
    message: str = field(init=False)
    args: tuple = field(init=False)

    def __post_init__(self):
        ex = "equations" if self.equations > self.variables else "variables"
        delta = abs(self.variables - self.equations)
        self.message = (
            f"Non-square {self.name}: {self.variables} variables vs. "
            f"{self.equations} equations; {delta} too many {ex}."
        )
        self.args = (self.variables, self.equations, self.message)