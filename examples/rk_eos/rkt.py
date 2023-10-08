"""Common code for the examples in this folder"""

# external modules
from yaml import load, SafeLoader
from casadi import SX

# internal modules
from simu.thermo import ThermoFactory, ThermoFrame, HelmholtzState
from simu.thermo import all_contributions


class ThermoNode(dict):
    """Class sharing a frame with other nodes of same type, and keeping
    the symbolic results by being a dictionary"""

    def __init__(self, frame: ThermoFrame):
        self.frame = frame
        state = SX.sym('x', 2 + len(self.frame.species))
        # TODO: Here need to call frame with parameter symbols instead.
        dict.__init__(self, zip(self.frame.property_names, self.frame(state)))


class MyThermoFactory(ThermoFactory):
    """This factory holds a ``ThermoFrame`` object for each defined model,
    ready to be deployed in a node."""

    def __init__(self):
        """Default and only constructor"""
        ThermoFactory.__init__(self)
        self.register(*all_contributions)
        self.register_state_definition(HelmholtzState)

        with open("parameters.yml", encoding='UTF-8') as file:
            parameters = load(file, SafeLoader)
        with open("structures.yml", encoding='UTF-8') as file:
            structures = load(file, SafeLoader)
        with open("frame_definitions.yml", encoding='UTF-8') as file:
            config = load(file, SafeLoader)

        def create_frame(cfg):
            structure = cfg.pop("structure")
            param = cfg.pop("parameters")
            cfg.update(structures[structure])
            frame = self.create_frame(cfg)
            frame.parameters.set_struct_values(parameters[param])
            return frame

        # register frames for all defined configurations
        self.frames = {name: create_frame(cfg)
                       for name, cfg in config.items()}

    def create_node(self, config_name: str) -> ThermoNode:
        """Based on the name of the thermodynamic model (ThermoFrame),
        create a wrapper object that hosts the properties for the instance."""
        return ThermoNode(self.frames[config_name])

def relax(frames: dict, result: dict, delta_x: list) -> float:
    """find relaxation factor"""
    gamma = 0.9  # relaxation distance
    alpha = 1 / gamma
    idx = 0
    for name, frame in frames.items():
        res = result[name]
        length = res["state"].rows()
        new_alpha = frame.relax(res.values(), delta_x[idx:idx + length])
        alpha = min(alpha, new_alpha)
        idx += length
    return float(gamma * alpha)
