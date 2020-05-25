#Wwing = 0.0051(Wdg + Nz)0.557 S0.649 w A0.5 (t/c)−0.4 root (1 + λ)0.1 ·(cosΛ)−1.0 S0.1 csw
# With Wdg the design gross weight, Nz the maximum load factor, SW the wing surface,
# A the aspect ratio (t/c)root the thickness ratio at the root, λ the taper ratio,
# Λ the wing sweep angle, and Scsw the control surface area.
from math import radians


def W_wing(Wdg, Nz, Sw, A, t_over_c_root, taper, sweep, Scsw):
    # Wdg = design gross weight, Nz = max load factor, Sw = wing surface, A = aspect ratio, Scsw = control surface area
    wing_weight = 0.0051 * (Wdg + Nz)**0.557 * Sw**0.649 A**0.5 * t_over_c_root**-0.4 *(1+taper)**0.1 * (cos(radians(sweep)))**-1 * Scsw**0.1
    return wing_weight