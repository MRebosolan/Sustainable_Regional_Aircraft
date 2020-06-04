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

# location of the systems from nose
x_mg = 15 #main gear from nose in m
x_ng = 4 #nose gear from nose in m
x_wg = 7 #wing group from nose in m
x_em = 18 #empennage group from nose in m
mtow = 35000 #maximum take of weight
w_wg = 2300 #weight of wing group
w_emp = 23000 #weight of empennage group






