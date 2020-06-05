import input
import numpy as np
import matplotlib.pyplot as plt

"""
Written by Matteo & Manuel
LOADS ACTING ON FUSELAGE: pressure vessel stresses, fuselage's own weight (assumed uniformly distributed
wing weight, empennage weight, landing gear reaction forces.

#ASSUMPTIONS
- Uniform longitudinal weight distribution
- Tail lift applied at end of aircraft length
"""

P_c = input.P_c
r = input.widthf/2
z_vl = input.bv/2 #point of action of vertical tail lift, estimated at half vert.tail height
t_options = np.linspace(0, 0.30, 1000)[1:]
W_cruise = 0.985 * input.MTOW
ac_length = 30 #dummy value
x_ac = input.x_ac
lh = input.lh

x_array = np.linspace(0, ac_length, 1000)



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


def shear_flow_due_to_empennage(rudder_load, enclosed_area, z_vl=z_vl):
    #rudder load positive towards right wing
    torque = -rudder_load*z_vl #right hand positive
    shear_flow = torque/(2*enclosed_area) #right hand positive
    return shear_flow

def internal_shear_and_moment_longitudinal(x, W= 322111, lf=ac_length, x_ac=x_ac, xh=ac_length):
    w_dist = W/lf
    L = ((0.5*lf)/(lf-x_ac))*W
    Lh = W-L
    if x>x_ac:
        d1, d2= 1,x-x_ac
    else:
        d1, d2=0,0
    if x>=xh:
        d3,d4 = 1, x-xh
    else:
        d3,d4=0,0

    shear_at_x = -w_dist*x + L*d1 + Lh*d3
    moment_at_x = -w_dist*(0.5*x**2) + L*d2 + Lh*d4
    return shear_at_x, moment_at_x

moments=[]
shears=[]
for x in x_array:
    moments.append(internal_shear_and_moment_longitudinal(x)[1])
    shears.append(internal_shear_and_moment_longitudinal(x)[0])

plt.plot(x_array, moments)
plt.plot(x_array, shears)
plt.show()

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
#Reaction_MG = (w_emp(x_emp-x_ng)+w_wg(x_wg-x_ng)+(mtow-w_emp-w_wg)*(fuselage_length/2-x_ng))/(x_mg-x_ng)
#Reaction force at the nose landing gear
#Reaction_NG = mtow-Reaction_MG






