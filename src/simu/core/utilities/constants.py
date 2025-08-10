"""Module containing constants as quantities.
"""

# stdlib
from math import pi as _pi

# internal
from .quantity import Quantity as _Q

PI = _Q(_pi)  #: The famous pi
ALPHA = _Q(7.2973525693e-3)  #: Fine structure constant
EPS_0 = _Q("8.8541878128e-12 F/m")  #: Vacuum electric permittivity
E_0 = _Q("1.602176634e-19 C")  #: Elementary charge
F = _Q("96485.33212 C/mol")  #: Faraday constant
GAMMA_G = _Q("6.67430e-11 m**3/(kg*s**2)")  #: Newtonian gravitational constant
H_PLANCK = _Q("6.62607015e-34 J*s")  #: Planck constant
K_B = _Q("1.380649e-23 J/K")  #: Boltzmann constant
MU_0 = _Q("1.25663706212 N/A**2")  #: Vacuum magnetic permittivity
M_E = _Q("9.1093837139e-31 kg")  #: Electron mass
M_N = _Q("1.67492750056e-27 kg")  #: Neutron mass
M_P = _Q("1.67262192595e-27 kg")  #: Proton mass
N_A = _Q("6.02214076e23  1/mol")  #: Avogadro constant
R_B = _Q("5.29177210544e-11 m")  #: Bohr radius
R_GAS = _Q("8.31446261815324 J/(mol*K)")  #: Molar gas constant
R_INF = _Q("10973731.568157 1/m")  #: Rydberg constant
SIGMA = _Q("5.670374419e-8 W/(m**2*K**4)")  #: Stefan-Boltzmann constant
STD_GRAVITY = _Q("9.80665 m/s**2")  #: Standard acceleration of gravity
V_LIGHT = _Q("299792458.0 m/s")  #: Speed of light