from abc import abstractmethod
from typing import Any
from inspect import currentframe, getargvalues

from .types import Map


class Configurable:
    """In particular the numerical solvers offer a number of options that can be
    provided by the constructor but also explicitly set. This class offers such
    functionality and supports validation of the individual options given."""

    __DEFAULT_VALIDATION = {"f": lambda x: True, "msg": ""}

    def __init__(self, exclude: list[str] = None):
        """
        The constructor uses the ``inspect`` module to obtain the arguments of
        the subclass constructor and collects them into the ``options``
        attribute after validation.

        Arguments given in ``exclude`` are not included here.
        """
        if exclude is None:
            exclude = []
        frame = currentframe().f_back
        args, _, _, values = getargvalues(frame)
        self._args = args[1:]
        self.options = {}
        for k in self._args:
            if k not in exclude:
                self.set_option(k, values[k])

    def set_option(self, name: str, value: Any):
        """
        Set the option with given ``name`` to the given ``value``. A
        ``KeyError`` is raised if there is no argument with given name, and a
        ``ValueError`` is raised if the value is not accepted.

        Only options defined as arguments in the constructor can be addressed
        with this method.
        """
        if name not in self._args:
            raise KeyError(f"{name} is not a possible option")
        validation = self._arg_validations.get(name, self.__DEFAULT_VALIDATION)
        if "replace_none" in validation and value is None:
            value = validation["replace_none"]
        if not validation["f"](value):
            raise ValueError(f"Option {name} is invalid: {validation["msg"]}")
        self.options[name] = value

    @property
    @abstractmethod
    def _arg_validations(self) -> Map[Map[Any]]:
        """The subclass must here provide a dictionary that maps each option
        name to a dictionary with the following keys:

        - ``f``: The validation function that is called with the value of the
          option and is to return ``True`` if the value is acceptable.
        - ``msg``: A specific error message to be included in the ``ValueError``
          thrown if ``f`` returned ``False``
        - ``replace_none`` (optional): A replacement value in case the value
          is ``None``. This can be used for mutable or complex default
          arguments.
        """
        ...

    @staticmethod
    def _validate_between(low: float , high: float):
        """This function can be used to define an argument validation
         between ``low`` and ``high``. It returns a dict with appropriate
         entries for ``f`` and ``msg`` as described for :meth:`_arg_validations`
         """
        return {
            "f": lambda x: low <= x <= high,
            "msg": f"must be between {low} and {high}"
        }
