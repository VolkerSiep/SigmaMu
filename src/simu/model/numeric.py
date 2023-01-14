"""This module implements functionality concerning the numerical handling
of the top model instance."""
from typing import TYPE_CHECKING

if TYPE_CHECKING:  # avoid circular dependencies just for typing
    from .base import ModelProxy


class NumericHandler:
    """This class implements the function object describing the top level
    model."""

    def __init__(self, model: "ModelProxy"):
        self.model = model
        self.__make_function()

    def __make_function(self):
        """Create a function that has the following arguments, each of them as
        a flat dictionary:

            - Model Parameters
            - Material States
            - Thermodynamic Parameters

        The result of the function consists likewise of

            - Properties
            - Residuals

        All the data is to be collected from the model and all child model
        proxies. For child models, only the free parameters are to be
        collected.
        """

    # TODO:
    # query state, model parameter, thermodynamic parameter, property and
    # residual names

    # all the following as quantity dictionaries:
    #   - export / import current states and parameters (thermo & model)
    #   - export current properties

    # in utilities, provide easy functions to serialise quantities or
    # dictionaries of them, e.g. "1.3 cm**2/s" should be parsed nicely for
    # yaml or json. Format with f"{a:~.16g}" to serialise.
