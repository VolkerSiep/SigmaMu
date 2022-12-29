"""This module contains some utility functionality specifically for the
process modelling functionality"""
from enum import Enum, auto


class ModelStatus(Enum):
    """Enumeration class to define status of model."""

    INTERFACE = auto()
    """The model is currently calling its :meth:`interface` method"""

    DEFINE = auto()
    """The model is currently calling its :meth:`define` method"""

    FUNCTION = auto()
    """The model is currently creating its local function representation"""

    READY = auto()
    """The model is defined and ready for being instantiated. For a model
    instance, this status describes that the model is ready to be integrated
    into a parent context, that is by connecting parameters and materials, and
    by using its properties"""

    FINALISED = auto()
    """The model instance is finalised, and it's local function is already
    integrated into the overall model context. Further connections cannot be
    established, and querying properties now returns the symbol in the overall
    model context."""

    def assure(self, desired: "ModelStatus", what: str):
        """Assure that current status is also the ``desired``
        status. Throw an exception otherwise, including the phrase ``what``
        describing the action::

            >>> stat = ModelStatus.DEFINE
            >>> stat.assure(ModelStatus.READY, "doing fun stuff")
            Traceback (most recent call last):
                ...
            RuntimeError: Status READY required for doing fun stuff
        """
        if not self is desired:
            raise RuntimeError(f"Status {desired.name} required for {what}")