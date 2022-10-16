class Model:
    """
    This class enables to create a numerical object that hosts the state,
    the thermodynamic models, the residuals, model and thermodynamic parameters.
    It also provides the functionality to provide Jacobian matrices, relax and
    update the state variables and update parameters (thermo and model).

    It does not implement solving algorithms.

    .. todo::
        There must be some kind of magic that assigns the correct physical
        dimensions to the properties defined by the thermodynamic models.

    The model defines all properties as quantities with physical dimensions.


    """

    def __init__(self):
        pass