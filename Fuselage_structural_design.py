import input


P_c = input.P_c
d = input.widthf

def pressure_vessel_stresses(P=P_c, d=d, t):
    hoop_stress = (P*d)/t
    longitudinal_stress = 0.5 * hoop_stress
    return hoop_stress, longitudinal_stress

