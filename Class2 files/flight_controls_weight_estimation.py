#flight control


def flight_controls(MTOW):
    kfc = 0.64 # powered flight controls
    factor = 1+ 0.2 + 0.15 # factor for leading edge + spoilers
    return kfc*factor *MTOW**(2/3)

def flight_controls_gd(rho_zero, V_dive, MTOW):
    q_d = V_dive*V_dive*rho_zero / 2
    factor = 1+ 0.2 + 0.15
    return 56.01* factor * ((MTOW*q_d/100000)**0.576)