import input
import numpy as np

P_c = input.P_c
r = input.widthf/2

t_options = np.linspace(0, 0.30, 1000)[1:]

#LOADS ACTING ON FUSELAGE: pressure vessel stresses, fuselage's own weight (assumed uniformly distributed
#wing weight, empennage weight, landing gear reaction forces.



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
x_emp = 18 #empennage group from nose in m

#weight of systems required
mtow = 35000 #maximum take of weight
w_wg = 2300 #weight of wing group
w_emp = 23000 #weight of empennage group

#fuselage length required
fuselage_length = 40




#Reaction force at the main gear
Reaction_MG = (w_emp(x_emp-x_ng)+w_wg(x_wg-x_ng)+(mtow-w_emp-w_wg)*(fuselage_length/2-x_ng))/(x_mg-x_ng)
#Reaction force at the nose landing gear
Reaction_NG = mtow-Reaction_MG






