import input
import numpy as np

P_c = input.P_c
r = input.widthf/2

t_options = np.linspace(0, 0.30, 1000)[1:]


def pressure_vessel_stresses(t, P=P_c, r=r):
    hoop_stress = (P*r)/t
    longitudinal_stress = 0.5 * hoop_stress
    return hoop_stress, longitudinal_stress

hoop_stresses=[]
longitudinal_stresses=[]

for t in t_options:
    hoop, long = pressure_vessel_stresses(t)
    hoop_stresses.append(hoop)
    longitudinal_stresses.append(long)





