#Wwing = 0.0051(Wdg + Nz)0.557 S0.649 w A0.5 (t/c)−0.4 root (1 + λ)0.1 ·(cosΛ)−1.0 S0.1 csw
# With Wdg the design gross weight, Nz the maximum load factor, SW the wing surface,
# A the aspect ratio (t/c)root the thickness ratio at the root, λ the taper ratio,
# Λ the wing sweep angle, and Scsw the control surface area.

def W_wing(Wdg, Nz, Sw, A, t_over_c_root, taper, sweep, Scsw):
    wing_weight = 0.00
