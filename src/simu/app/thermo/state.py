# external modules
from casadi import vertsplit, vertcat

# internal modules
from simu.core.thermo.state import StateDefinition, register
from simu.core.utilities.quantity import Quantity, base_magnitude, base_unit
from simu.core.utilities.types import MutMap


@register
class HelmholtzState(StateDefinition):
    """This definition interprets the state as being temperature, volume,
    and mole numbers. Accordingly, it defines:

    ======== ============================
    Property Description
    ======== ============================
    ``T``    Temperature
    ``V``    Volume
    ``n``    Mole vector
    ======== ============================
    """
    def prepare(self, result: MutMap[Quantity], flow: bool = False):
        state = result["_state"].magnitude
        result["T"], result["V"], *n_vec = vertsplit(state, 1)
        result["n"] = vertcat(*n_vec)
        s = "/s" if flow else ""
        for name, unit in [("T", "K"), ("V", f"m**3{s}"), ("n", f"mol{s}")]:
            result[name] = Quantity(result[name], base_unit(unit))

    def reverse(self, state):
        return [base_magnitude(state.temperature), None] + \
            list(base_magnitude(state.mol_vector))

    def declare_vector_keys(self, species):
        return {"n": list(species.keys())}


@register
class GibbsState(StateDefinition):
    """This definition interprets the state as being temperature, pressure,
    and mole numbers. Accordingly, it defines:

    ======== ============================
    Property Description
    ======== ============================
    ``T``    Temperature
    ``p``    Pressure
    ``n``    Mole vector
    ======== ============================
    """
    def prepare(self, result: MutMap[Quantity], flow: bool = False):
        state = result["_state"].magnitude
        result["T"], result["p"], *n_vec = vertsplit(state, 1)
        result["n"] = vertcat(*n_vec)
        q_unit = "mol/s" if flow else "mol"
        for name, unit in [("T", "K"), ("p", "Pa"), ("n", q_unit)]:
            result[name] = Quantity(result[name], base_unit(unit))

    def reverse(self, state):
        return [base_magnitude(state.temperature),
                base_magnitude(state.pressure)] + \
                list(base_magnitude(state.mol_vector))

    def declare_vector_keys(self, species):
        return {"n": list(species.keys())}
