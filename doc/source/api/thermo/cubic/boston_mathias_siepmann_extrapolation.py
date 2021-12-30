#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Compare the AspenTech published Mathias extrapolation of the Boston-Mathias
alpha function with one self-derived, using the same mathematical frame, but
with a smooth first and second derivative at the critical temperature.
"""

from numpy import exp, linspace, sqrt
import pylab

XMIN, XMAX = 0.5, 6.0

omega = 0.25
eta = 0.063  # freely invented (= 0.065 for NH3, 0.13 for H2O)

m = 0.48 + (1.57 - 0.17 * omega) * omega
d = 1 + m / 2 + 0.3 * eta/2
c = 1 - 1 / d

# my own c and d parameters for smooth second derivative
c2 = m + 0.3 * eta
d2 = 1 + c2 + 4 * eta / c2
c2 = c2 / d2   # transform back from my parameters
d2 = d2 / 2

tau = linspace(XMIN, XMAX)

alpha_1 = 1 + m * (1 - sqrt(tau)) - eta * (1 - sqrt(tau)) * (0.7 - tau)
alpha_2 = exp(c* (1 - tau ** d))
alpha_3 = exp(c2* (1 - tau ** d2))

dadt_1 = -m / 2 / sqrt(tau) + eta / 2 / sqrt(tau) * (0.7 - tau) + eta * (1 - sqrt(tau))
dadt_2 =  -exp(c* (1 - tau ** d)) * c * d * tau ** (d - 1)
dadt_3 =  -exp(c2* (1 - tau ** d2)) * c2 * d2 * tau ** (d2 - 1)

fig = pylab.figure(figsize=(6,4), dpi=300)
ax1 = fig.subplots()
ax2 = ax1.twinx()

ax1.plot(tau, alpha_1, "k-", label="subcritical")
ax1.plot(tau, alpha_2, "b-", label="supercritical BM")
ax1.plot(tau, alpha_3, "g-", label="supercritical BMS")

ax2.plot(tau, dadt_1, "k--", label="subcritical")
ax2.plot(tau, dadt_2, "b--", label="supercritical BM")
ax2.plot(tau, dadt_3, "g--", label="supercritical BMS")
ax1.grid()
ax1.legend(loc=4)
ax1.set_xlabel("Reduced temperature [-]")
ax1.set_ylabel("square root of alpha [-]")
ax2.set_ylabel("first derivative [-]")
ax1.set_ylim([0, None])
ax1.set_xlim([XMIN, XMAX])
pylab.savefig("boston_mathias_siepmann_extrapolation.png", bbox_inches="tight")
pylab.show()
