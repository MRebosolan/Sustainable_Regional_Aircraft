#Wwing = 0.0051(Wdg + Nz)0.557 S0.649 w A0.5 (t/c)−0.4 root (1 + λ)0.1 ·(cosΛ)−1.0 S0.1 csw
# With Wdg the design gross weight, Nz the maximum load factor, SW the wing surface,
# A the aspect ratio (t/c)root the thickness ratio at the root, λ the taper ratio,
# Λ the wing sweep angle, and Scsw the control surface area.
from math import *
import numpy as np

def W_wing(W_zfw, b, half_sweep, n_ult, S, t_r):
    x = 0.0017 * W_zfw * (b / (cos(half_sweep)))**0.75
    y = 1 + sqrt(6.3 * cos(half_sweep) / b)
    z = (n_ult**0.55) * (b * S / (t_r * W_zfw * cos(half_sweep)))**0.3
    return x*y*z

weight = W_wing(44245+18055, 76, np.radians(27), 3.5*1.5, 760, 2.16)
print (weight)
weight_kg = W_wing(20069+8190, 23.2, np.radians(27), 3.5*1.5, 70.6, 0.66)
print(weight_kg)


def gd_wing(MTOW, AR, half_sweep, n_ult, S, t_over_c, taper, mach_h):
    upper = 0.00428*(S**0.48) * AR * (mach_h**0.43) * ((MTOW * n_ult)**0.84) * (taper**0.14)
    lower = ((100*t_over_c)**0.76)*(np.cos(half_sweep)**1.54)
    return upper/lower

gdweight = gd_wing(75000, 7.62, np.radians(27), 3.5*1.5, 760, 0.1, 0.3, 0.5)
print (gdweight)    