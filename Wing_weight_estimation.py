"""
input: zero fuel weight, wingspan, half chord sweep, ultimate load factor, wing area, thickness at root
MTOW, aspect ratio, thickness over chord, taper ratio, design mach number at sea level
output: two types of wing weight calculation
responsible: Jorn

"""

#Wwing = 0.0051(Wdg + Nz)0.557 S0.649 w A0.5 (t/c)−0.4 root (1 + λ)0.1 ·(cosΛ)−1.0 S0.1 csw
# With Wdg the design gross weight, Nz the maximum load factor, SW the wing surface,
# A the aspect ratio (t/c)root the thickness ratio at the root, λ the taper ratio,
# Λ the wing sweep angle, and Scsw the control surface area.
from math import *
import numpy as np


# def W_wing(Wdg, Nz, Sw, A, t_over_c_root, taper, sweep, Scsw):
#     # Wdg = design gross weight, Nz = max load factor, Sw = wing surface, A = aspect ratio, Scsw = control surface area
#     wing_weight = 0.0051 * (Wdg + Nz)**0.557 * Sw**0.649 A**0.5 * t_over_c_root**-0.4 *(1+taper)**0.1 * (cos(radians(sweep)))**-1 * Scsw**0.1
#     return wing_weight

def W_wing(W_zfw, b, half_sweep, n_ult, S, t_r):
    x = 0.0017 * W_zfw * (b / (cos(half_sweep)))**0.75
    y = 1 + sqrt(6.3 * cos(half_sweep) / b)
    z = (n_ult**0.55) * (b * S / (t_r * W_zfw * cos(half_sweep)))**0.3
    factor = 1 + 0.02-0.05 - 0.05 # + 0.02 for spoilers. - 0.05 for wing mounted engines. -0.05 for non wing mounted gear. fowler flaps: add 0.02
    return x*y*z *factor

# weight = W_wing(44245+18055, 76, np.radians(27), 3.5*1.5, 760, 2.16)
# print (weight)
# weight_kg = W_wing(20069+8190, 23.2, np.radians(27), 3.5*1.5, 70.6, 0.66)
# print(weight_kg)


def gd_wing(MTOW, AR, half_sweep, n_ult, S, t_over_c, taper, mach_h):
    upper = 0.00428*(S**0.48) * AR * (mach_h**0.43) * ((MTOW * n_ult)**0.84) * (taper**0.14)
    lower = ((100*t_over_c)**0.76)*(np.cos(half_sweep)**1.54)
    factor = 1 + 0.02-0.05 - 0.05 # + 0.02 for spoilers. - 0.05 for wing mounted engines. -0.05 for non wing mounted gear. fowler flaps: add 0.02

    return upper*factor/lower

# gdweight = gd_wing(75000, 7.62, np.radians(27), 3.5*1.5, 760, 0.1, 0.3, 0.5)
# print (gdweight)