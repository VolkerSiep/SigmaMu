from simu.core.utilities import SymbolQuantity, base_unit


# auxiliary functions
def sym(name, units):
    """Return a scalar symbol of given name and units"""
    return SymbolQuantity(name, base_unit(units))


def vec(name, size, units):
    """Return a vector symbol of given name, units, and size"""
    return SymbolQuantity(name, base_unit(units), size)