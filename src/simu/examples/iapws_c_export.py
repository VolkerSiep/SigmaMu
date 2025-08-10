from yaml import safe_load
from casadi import SX, vertcat, jacobian, Function

from simu import (SpeciesDefinition, parse_quantities_in_struct,
                  flatten_dictionary)
from simu.app import RegThermoFactory, DATA_DIR

PHASE = "Gas"


def create_iapws(phase: str):
    fac = RegThermoFactory()
    contributions = [
        "MolecularWeight", "ReducedStateIAPWS", "StandardStateIAPWS",
        "IdealGasIAPWS", "Residual1IAPWS", "Residual2IAPWS",
        "Residual3IAPWS", "Residual4IAPWS", f"{phase}IAPWSIdealMix"
    ]

    config = {
        "species": ["H2O"],
        "state": "HelmholtzState",
        "contributions": contributions
    }
    species_def =  {"H2O": SpeciesDefinition("H2O")}
    frame = fac.create_frame(species_def, config)
    with open(DATA_DIR / "parameters_iapws.yml") as file:
        params = parse_quantities_in_struct(safe_load(file)["data"])
    params = {n: v for n, v in params.items() if n in contributions}
    return frame, params


def main():
    frame, params = create_iapws(PHASE)
    T, V, n = [SX.sym(i) for i in "TVn"]
    state = vertcat(T, V, n)

    prop = frame(state, params, squeeze_results=False)
    prop = flatten_dictionary(prop["props"])
    prop = {key: value.magnitude for key, value in prop.items()}

    prop["a_t"] = -prop["S"]
    prop["a_v"] = -prop["p"]
    prop["a_n"] = prop["mu"]
    prop["a_tt"] = jacobian(prop["a_t"], T)
    prop["a_tv"] = jacobian(prop["a_t"], V)
    prop["a_tv"] = prop["a_tv"]
    prop["a_tn"] = jacobian(prop["a_n"], T)
    prop["a_nt"] = prop["a_tn"].T
    prop["a_vn"] = jacobian(prop["a_n"], V)
    prop["a_nv"] = prop["a_tn"].T
    prop["a_nn"] = jacobian(prop["a_n"], n)

    for i in "TVn":
        del prop[i]

    f = Function(f"A_iapws_{PHASE.lower()}", [T, V, n], list(prop.values()),
                 ["T", "V", "n"], list(prop.keys()))

    opts = {"with_header": True}
    f.generate(f"iapws_c_export_{PHASE.lower()}.c", opts)


if __name__ == '__main__':
    main()
