#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Extract p(T, V) data for ammonia at constant V and calculate back the
alpha function (for different V, i.e. densities).
Then fit the sub-critical alpha function and observe the super-critical
extrapolation.

Observation: The curve doesn't fit at all even for sub-critical temperatures,
and no polar parameter can make it. We adapted C to fit the volume at T_c,
but that does not improve the fit otherwise. This applies both for moderate
pressure levels (V = 1e-3 m3/mol) and for pressures near critical
(1/V = 13211 mol/m3).
"""

from numpy import array, sqrt, exp, linspace
from matplotlib import pyplot

# gas constant
R_GAS = 8.31446261815324

# ammonia fixed properties
T_C = 405.40  # Critical temperature [K]
P_C = 113.33e5  # Critical pressure [Pa]
OMEGA = 0.25601  # acentric factor [-]
ETA = 0.063  # polar parameter [-]
C = 1e-5  # volume shift parameter [m3/mol]

# constant properties
V = 1e-3  # m3/mol

# data from NIST chemistry webbook as [T [K], p[bar]]
# https://webbook.nist.gov/chemistry/fluid/
T, P = array([
    [330.00, 22.355],
    [335.00, 22.917],
    [340.00, 23.470],
    [345.00, 24.014],
    [350.00, 24.550],
    [355.00, 25.080],
    [360.00, 25.604],
    [365.00, 26.122],
    [370.00, 26.635],
    [375.00, 27.143],
    [380.00, 27.648],
    [385.00, 28.148],
    [390.00, 28.645],
    [395.00, 29.139],
    [400.00, 29.630],
    [405.00, 30.118],
    [410.00, 30.604],
    [415.00, 31.087],
    [420.00, 31.568],
    [425.00, 32.047],
    [430.00, 32.524],
    [435.00, 33.000],
    [440.00, 33.473],
    [445.00, 33.945],
    [450.00, 34.416],
    [455.00, 34.885],
    [460.00, 35.353],
    [465.00, 35.819],
    [470.00, 36.284],
    [475.00, 36.748],
    [480.00, 37.211],
    [485.00, 37.673],
    [490.00, 38.134],
    [495.00, 38.594],
    [500.00, 39.053],
    [510.00, 39.969],
    [520.00, 40.881],
    [530.00, 41.791],
    [540.00, 42.698],
    [550.00, 43.602],
    [560.00, 44.504],
    [570.00, 45.403],
    [580.00, 46.301],
    [590.00, 47.197],
    [600.00, 48.091],
    [610.00, 48.983],
    [620.00, 49.874],
    [630.00, 50.764],
    [640.00, 51.652],
    [650.00, 52.538],
    [660.00, 53.423],
    [670.00, 54.308],
    [680.00, 55.191],
    [690.00, 56.073],
    [700.00, 56.953]]).T

V2 = 1 / 13211.8  # m3/mol,  = critical volume
T2, P2 = array([
    [410.00, 122.54],
    [420.00, 142.58],
    [430.00, 162.81],
    [440.00, 183.17],
    [450.00, 203.62],
    [460.00, 224.13],
    [470.00, 244.67],
    [480.00, 265.21],
    [490.00, 285.75],
    [500.00, 306.27],
    [510.00, 326.75],
    [520.00, 347.20],
    [530.00, 367.59],
    [540.00, 387.93],
    [550.00, 408.22],
    [560.00, 428.44],
    [570.00, 448.61],
    [580.00, 468.71],
    [590.00, 488.75],
    [600.00, 508.72],
    [610.00, 528.62],
    [620.00, 548.47],
    [630.00, 568.24],
    [640.00, 587.96],
    [650.00, 607.61],
    [660.00, 627.19],
    [670.00, 646.72],
    [680.00, 666.18],
    [690.00, 685.59],
    [700.00, 704.93]]).T

V = 1e-2  # m3/mol
T, P = array([
    [280.00, 2.2558],
    [300.00, 2.4319],
    [320.00, 2.6055],
    [340.00, 2.7775],
    [360.00, 2.9483],
    [380.00, 3.1185],
    [400.00, 3.2880],
    [420.00, 3.4572],
    [440.00, 3.6261],
    [460.00, 3.7946],
    [480.00, 3.9630],
    [500.00, 4.1312],
    [520.00, 4.2993],
    [540.00, 4.4672],
    [560.00, 4.6350],
    [580.00, 4.8027],
    [600.00, 4.9703],
    [620.00, 5.1379],
    [640.00, 5.3054],
    [660.00, 5.4729],
    [680.00, 5.6403],
    [700.00, 5.8076],
    [720.00, 5.9750]]).T

def alpha_BMS(tau):
    # my own c and d parameters for smooth second derivative
    m = 0.48 + (1.57 - 0.17 * OMEGA) * OMEGA
    c = m + 0.3 * ETA
    d = 1 + c + 4 * ETA / c
    return exp(c / d * (1 - tau ** (d / 2)))


def alpha_BM(tau):
    m = 0.48 + (1.57 - 0.17 * OMEGA) * OMEGA
    d = 1 + m / 2 + 0.3 * ETA / 2
    c = 1 - 1 / d
    return exp(c * (1 - tau ** d))


def alpha_sub(tau):
    m = 0.48 + (1.57 - 0.17 * OMEGA) * OMEGA
    return 1 + m * (1 - sqrt(tau)) - ETA * (1 - sqrt(tau)) * (0.7 - tau)


def alpha_data():
    x = 2 ** (1 / 3) - 1
    omega_a = 1 / (9 * x)
    omega_b = x / 3
    # calculate B-parameter
    b = omega_b * R_GAS * T_C / P_C

    # convert pressure to Pa
    p = 1e5 * P

    a = (R_GAS * T / (V - b + C) - p) * (V + C) * (V + b + C)
    alpha_factor = omega_a * R_GAS ** 2 * T_C ** 2 / P_C
    return (a / alpha_factor) ** 0.5

def main():
    T_r = T / T_C
    pyplot.plot(T_r, alpha_data(), "k.", label="NIST data")
    T_r = linspace(0.001, 5, num=1000)
    pyplot.plot(T_r, alpha_sub(T_r), label="Mathias alpha function")
    pyplot.plot(T_r, alpha_BM(T_r), label="Boston-Mathias extrapolation")
    pyplot.plot(T_r, alpha_BMS(T_r), label="Modified extrapolation")
    pyplot.grid()
    pyplot.legend(loc="best")
    pyplot.xlabel("$T / T_c$ [-]")
    pyplot.xlim([0, None])
    pyplot.ylim([0, None])
    pyplot.ylabel(r"$\sqrt{\alpha}$ [-]")
    pyplot.show()



if __name__ == "__main__":
    main()
