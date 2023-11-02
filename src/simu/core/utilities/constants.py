#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Module containing constants as quantities.
"""

from math import pi as _pi
from .quantity import Quantity as Qty

PI = Qty(_pi)  # the famous pi
R_GAS = Qty("8.31446261815324 J/(mol*K)")  # gas constant
V_LIGHT = Qty("299792458.0 m/s")  # speed of light
H_PLANCK = Qty("6.62607015e-34 J*s")  # Planck constant
GAMMA_G = Qty("6.67430e-11 m**3/(kg*s**2)")  # Newtonian gravitational constant
N_A = Qty("6.02214076e23  1/mol")  # Avogadro constant
ALPHA = Qty(7.2973525693e-3)  # Fine structure constant
SIGMA = Qty("5.670374419e-8 W/(m**2*K**4)")  # Stefan-Boltzmann constant
F = Qty("96485.33212 C/mol")  # Faraday constant
EPS_0 = Qty("8.8541878128e-12 F/m")  # Vacuum electric permittivity
MU_0 = Qty("1.25663706212 N/A**2")  # Vacuum magnetic permittivity
E_0 = Qty("1.602176634e-19 C")  # Elementary charge
K_B = Qty("1.380649e-23 J/K")  # Boltzmann constant

STD_GRAVITY = Qty("9.80665 m/s**2")  # standard acceleration of gravity
