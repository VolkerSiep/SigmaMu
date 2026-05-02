"""This module defines types of complex data structures"""

from collections.abc import Mapping, MutableMapping
from typing import Protocol, runtime_checkable, TYPE_CHECKING
from scipy.sparse import csr_array
from numpy.typing import NDArray


type Map[T] = Mapping[str, T]
"""A mapping of strings to another type"""

type NestedMap[T] = Map[T | "NestedMap[T]"]
"""A nested mapping of strings to another type"""

type MutMap[T] = MutableMapping[str, T]
"""A mutable mapping of strings to another type"""

type NestedMutMap[T] = MutMap[T | "NestedMutMap[T]"]
"""A nested mutable mapping of strings to another type"""


@runtime_checkable
class OutputIOStream(Protocol):
    """A Protocol class to describe anything that can write a string"""
    def write(self, line: str):  ...


@runtime_checkable
class LinearSolver(Protocol):
    """A Protocol class to describe a valid linear solver. A simplest
    implementation would be

    .. code-block::

        from scipy.sparse.linalg import spsolve

        class MySolver:
            def reset(self):
                pass

            def solve(self, matrix, rhs):
                return spsolve(matrix, rhs)
    """

    def reset(self):
        """In case the solver utilizes any pre-conditioning to enhance
        subsequent solving of similar matrices, this method can be used
        to signal a new run with a fresh matrix.
        """
        ...

    def solve(self, matrix: csr_array, rhs: NDArray) -> NDArray:
        """This is the actual solve method with the common signature as for
        ``scipy.sparse.linalg.spsolve``.
        """
        ...
