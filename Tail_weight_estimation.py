from math import radians, cos
#Torenbeek commercial transport airplanes weight estimation method

def hor_tail_weight(Kh, Sh, VD, half_chord_sweep_hor):
    #Already verified via hand calculation
    #Kh=1 for fixed incidence stabilizers, 1.1 for variable incidence stabilizers
    #VD= dive speed in KEAS

    hor_tail_weight= Kh * Sh * ((3.81*(Sh**0.2 * VD)) / (1000 * (cos(radians(half_chord_sweep_hor)))**0.5)-0.287)
    return hor_tail_weight

def hor_tail_weight(Kh, Sh, VD, half_chord_sweep_hor):
    #Already verified via hand calculation
    #Kh=1 for fixed incidence stabilizers, 1.1 for variable incidence stabilizers
    #VD= dive speed in KEAS

    hor_tail_weight= Kh * Sh * ((3.81*(Sh**0.2 * VD)) / (1000 * (cos(radians(half_chord_sweep_hor)))**0.5)-0.287)
    return hor_tail_weight