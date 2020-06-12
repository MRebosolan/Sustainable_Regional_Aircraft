import input
import numpy as np
import matplotlib.pyplot as plt
from math import pi

"""
Written by Matteo & Manuel
LOADS ACTING ON FUSELAGE: pressure vessel stresses, fuselage's own weight (assumed uniformly distributed +)
wing lift, horizontal tail lift

#ASSUMPTIONS
- Uniform longitudinal weight distribution
- Tail lift applied at end of aircraft length
- Cylindrical, thin walled fuselage
- Effect of vertical tail neglected
"""

P_c = input.P_c
P_cruise = input.P_cruise
r = input.widthf/2
z_vl = input.bv/2 #point of action of vertical tail lift, estimated at half vert.tail height
t_options = np.linspace(0, 0.30, 1000)[1:]
W_cruise = 0.985 * input.MTOW
ac_length = input.lf #dummy value
x_ac = 12 #estimate
lh = input.lh
widthf = input.widthf
n_ult = input.n_ult

x_array = np.linspace(0, ac_length, 1000)

#----------FUSELAGE INTERNAL MOMENTS/SHEAR FORCES-----------------



def internal_shear_and_moment_longitudinal(x, W= W_cruise, lf=ac_length, x_ac=x_ac, xh=ac_length):
    w_dist = W/lf
    L = ((0.5*lf)/(lf-x_ac))*W
    Lh = W-L
    if x > x_ac:
        d1, d2= 1,x-x_ac
    else:
        d1, d2=0,0
    if x >= xh:
        d3,d4 = 1, x-xh
    else:
        d3,d4=0,0

    shear_at_x = -w_dist*x + L*d1 + Lh*d3     #positive downwards
    moment_at_x = -w_dist*(0.5*x**2) + L*d2 + Lh*d4  #sagging positive
    return shear_at_x, moment_at_x

moments=[]
shears=[]
for x in x_array:
    moments.append(internal_shear_and_moment_longitudinal(x)[1] * n_ult)
    shears.append(internal_shear_and_moment_longitudinal(x)[0] * n_ult)





#----------FUSELAGE INTERNAL STRESSES--------------
#derived from internal moments/shears + pressure vessel stresses


def pressure_vessel_stresses(t, R, P_c, P_cruise):
    P = P_c - P_cruise
    hoop_stress = (P*R)/t
    longitudinal_stress = 0.5 * hoop_stress
    max_shear_stress = 0.5 * hoop_stress #Mohr's circle
    return hoop_stress, longitudinal_stress, max_shear_stress


def shear_stress_due_to_empennage_torque(rudder_load, enclosed_area, z_vl, t):
    #rudder load positive towards right wing
    torque = -rudder_load*z_vl #right hand positive
    shear_flow = torque/(2*enclosed_area) #right hand positive
    shear_stress = shear_flow/t
    return shear_stress


def area_moments_of_inertia(R, t):
    Iyy = pi*t*R**3
    Izz = pi*t*R**3
    return Iyy, Izz

def max_internal_bending_stress(R, Iyy, My):
    along_z = (My/Iyy)*R    #positive on top of fuselage
    return along_z

# bending_stresses_bottom = []
# bending_stresses_top = []
# for m in moments:
#     Iyy = area_moments_of_inertia(2.12, 0.10)
#     bending_stresses_bottom.append(max_internal_bending_stress(2.12, Iyy, m))
#     bending_stresses_top.append(-max_internal_bending_stress(2.12, Iyy, m))


def max_internal_vertical_shear(Iyy, Vz, R, t):
    #Verified by hand calculation
    #downward positive
    ri = R-t
    Q=(2/3)*(R**3-ri**3)
    tau_max = (Vz*Q)/(Iyy*t*2)
    return tau_max


def max_fuselage_stresses(t, R, x_array, P_c, P_cruise, z_vl, rudder_load, moments_y, shears_z):

    enclosed_area = pi*R**2
    bending_stresses_bottom = []
    bending_stresses_top = []
    shear_stresses = []
    hoop_stress, axial_stress, max_shear_stress = pressure_vessel_stresses(t, R, P_c, P_cruise)
    torque_shear_stress = shear_stress_due_to_empennage_torque(rudder_load, enclosed_area, z_vl, t)

    Iyy, Izz = area_moments_of_inertia(R, t)
    for x in x_array:
        x_index = np.where(x_array == x)[0][0]
        bending_stress_bottom = max_internal_bending_stress(R, Iyy, moments_y[x_index])
        bending_stress_top = -max_internal_bending_stress(R, Iyy, moments_y[x_index])
        shear_stress = max_internal_vertical_shear(Iyy, shears_z[x_index], R, t)
        bending_stresses_bottom.append(bending_stress_bottom)
        bending_stresses_top.append(bending_stress_top)
        shear_stresses.append(shear_stress)

    bending_stresses_top = np.array(bending_stresses_top) + axial_stress
    bending_stresses_bottom = np.array(bending_stresses_bottom) + axial_stress
    shear_stresses = np.array(shear_stresses) + torque_shear_stress
    #mohr_shear_stresses_top = []
    #mohr_shear_stresses_bottom = []
    #for s1 in bending_stresses_top:
    #    if s1>0:
    #        mohr_shear = max(hoop_stress, s1)/2
    #    else:
    #        mohr_shear = (hoop_stress-s1)/2
    #    mohr_shear_stresses_top.append(mohr_shear)
    #
    #for s2 in bending_stresses_bottom:
    #    if s1>0:
    #        mohr_shear = max(hoop_stress, s2)/2
    #    else:
    #        mohr_shear = (hoop_stress-s2)/2
    #    mohr_shear_stresses_bottom.append(mohr_shear)



    return bending_stresses_bottom, bending_stresses_top, shear_stresses

bmb, bmt, sh1 =(max_fuselage_stresses(0.001, widthf/2, x_array, P_c, P_cruise, z_vl, 0, moments, shears))

idx = moments.index(min(moments))


plt.plot(x_array, bmb)
plt.plot(x_array, bmt)
plt.plot(x_array, sh1)

# plt.plot(x_array, moments)
# plt.plot(x_array, shears)
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






