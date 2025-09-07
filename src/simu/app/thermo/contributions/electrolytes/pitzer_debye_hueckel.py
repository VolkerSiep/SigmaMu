from simu import (ThermoContribution, registered_contribution, Quantity,
                  N_A, E_0, EPS_0, K_B, R_GAS, PI, sqrt, log, qvertcat)
from simu.core.utilities.types import MutMap

_M0 = Quantity(1.0, "mol/kg")


@registered_contribution
class ElectrolyteBasics(ThermoContribution):
    r"""

    """
    def define(self, res):
        n = res["n"]

        # find index of water
        species_def = self.species_definitions
        water_name = self.options.get("water_name", "H2O")
        h2o_idx = self.species.index(water_name)
        m_h2o = n[h2o_idx] * species_def[water_name].molecular_weight

        res["charges"] = c = qvertcat(*[s.charge for s in species_def.values()])
        res["molalities"] = m = n / (m_h2o * _M0)
        res["I"] = m.T @ (c ** 2) / 2


@registered_contribution
class PitzerDebyeHueckel(ThermoContribution):
    r"""

    """
    def define(self, res):
        temp, ionic_strength = res["T"], res["I"]
        rho_w = self.par_scalar("RHO_WATER", "kg/m**3")
        eps_r = self.par_scalar("EPS_WATER", "dimless")
        b = self.par_scalar("B", "dimless")

        a_gamma = (
            sqrt(2 * PI * N_A * _M0 * rho_w)
            * (EPS_0 ** 2 / (4 * PI * EPS_0 * eps_r * K_B * temp)) ** (3 / 2))
        bsqi = b * (sq_i := sqrt(ionic_strength))
        res["_f_pdh"] = f = -4 / 3  * a_gamma * sq_i ** 3 * log(1 + bsqi) / bsqi

        # TODO: update in mu and S, using chain rule (see iapws)