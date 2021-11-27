# -*- coding: utf-8 -*-

# external modules
from casadi import dot, log, vertsplit, vertcat, sum1

# internal modules
from .contribution import ThermoContribution
from ..constants import R_GAS


class HelmholtzState(ThermoContribution):
    """x"""
    def define(self, res, par):
        res["T"], res["V"], *res["n"] = vertsplit(res["state"], 1)
        res["n"] = vertcat(*res["n"])

class GibbsState(ThermoContribution):
    """x"""
    def define(self, res, par):
        res["T"], res["p"], *res["n"] = vertsplit(res["state"], 1)
        res["n"] = vertcat(*res["n"])


class H0S0ReferenceState(ThermoContribution):
    """x"""
    @property
    def parameter_structure(self):
        t_s = ThermoContribution._tensor_structure
        return {"dh_form": t_s(self.species),
                "s_0": t_s(self.species),
                "T_ref": None,
                "p_ref": None}

    def define(self, res, par):
        vector = self._vector
        s_0 = vector(par["s_0"])
        dh_form = vector(par["dh_form"])
        res["S"] = dot(s_0, res["n"])
        res["mu"] = dh_form - res["T"] * s_0
        res["T_ref"] = par["T_ref"]
        res["p_ref"] = par["p_ref"]


class LinearHeatCapacity(ThermoContribution):
    """x"""
    @property
    def parameter_structure(self):
        t_s = ThermoContribution._tensor_structure
        return t_s(["cp_a", "cp_b"], self.species)

    def define(self, res, par):
        vector = self._vector
        T, n = res["T"], res["n"]
        T_ref = res["T_ref"]
        d_T, f_T = T - T_ref, T / T_ref

        cp_a, cp_b = vector(par["cp_a"]), vector(par["cp_b"])
        d_h = (cp_a + 0.5 * d_T * cp_b) * d_T
        d_s = (cp_a - cp_b * T_ref) * log(f_T) + cp_b * d_T
        res["S"] += dot(d_s, n)
        res["mu"] += d_h - T * d_s

    def relax(self, state, delta_state, parameters):
        T, d_T = state[0], delta_state[0]
        return -T / d_T if d_T < 0 else 100


class StandardState(ThermoContribution):
    """x"""
    def define(self, res, par):
        # tag current chemical potential and entropy as standard state
        res["S_std"] = res["S"]
        res["p_std"] = res["p_ref"]
        res["mu_std"] = res["mu"]


class IdealMix(ThermoContribution):
    """x"""
    def define(self, res, par):
        T, n = res["T"], res["n"]
        N = sum1(n)
        x = n / N
        gtn = R_GAS * log(x)
        res["S"] -= dot(n, gtn)
        res["mu"] += T * gtn

    def relax(self, state, delta_state, parameters):
        n, d_n = state[2:], delta_state[2:]
        cand = [-n_i / dn_i for n_i, dn_i in zip(n, d_n) if dn_i < 0]
        return min(cand) if cand else 100


class GibbsIdealGas(ThermoContribution):
    """x"""
    def define(self, res, par):
        T, p, n, p_ref = res["T"], res["p"], res["n"], res["p_ref"]
        N = sum1(n)
        gtn = R_GAS * log(p / p_ref)

        res["S"] -= N * gtn
        res["V"] = N * R_GAS * T / p
        res["mu"] += T * gtn

    def relax(self, state, delta_state, parameters):
        p, d_p = state[1], delta_state[1]
        return -p / d_p if d_p < 0 else 100


class HelmholtzIdealGas(ThermoContribution):
    def define(self, res, par):
        pass  # TODO: derive and implemenet

    def relax(self, state, delta_state, parameters):
        V, d_V = state[1], delta_state[1]
        return -V / d_V if d_V < 0 else 100

    def initial_state(self, T, p, n, parameters):
        V = sum(n) * R_GAS * T / p
        return [T, V] + list(n)