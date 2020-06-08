# Wingbox for bending stress
import numpy as np
import input
import Wingbox_design
import matplotlib.pyplot as plt
from math import *

x_array = Wingbox_design.x_array[2:]
moment_array = Wingbox_design.moment_array
moment_array2 = Wingbox_design.moment_array2

# coordinates for the airfoil
yield_stress_material = 200  # start of with steel for now

# Read from file

f1 = open('airfoil2.txt', 'r')
lines1 = f1.readlines()

xcoord1 = []
ycoord1 = []
ycoord2 = []

for i in range(0, 103):
    xcoord1.append(float(lines1[i].split()[0]))
    ycoord1.append(float(lines1[i].split()[1]))
    ycoord2.append(-1 * float(lines1[i].split()[1]))

xcoord1 = xcoord1[::-1]
ycoord1 = ycoord1[::-1]
ycoord2 = ycoord2[::-1]

# plt.figure(1)
# plt.grid(True, which="major", color="#999999")
# plt.grid(True, which="minor", color="#DDDDDD", ls="--")
# plt.minorticks_on()
# plt.plot(xcoord1, ycoord1, color='r')
# plt.plot(xcoord1, ycoord2, color='r')
# plt.xlim(0, 1)
# plt.ylim(-0.3, 0.3)
# plt.text(0.0, 0.0, 'LE')
# plt.text(1.0, 0.0, 'TE')
# plt.ylabel('y/c [-]')
# plt.xlabel('x/c [-]')
#
# # plt.show()
#------------------------------------------------------------------------------------------------------------------
#INPUTS

chord_length = 2  # chord length in meters
t_d = 0.01 #THICkNESS OF AIRFOIL
number_booms = 4 #NUMBER OF POINTS (10 POINTS = 16 BOOMS,4p=4b 5p=6b )#  (ON TOP SIDE FOR NOW)
moment_cs = 450000 #MOMENT OF CROSS SECTION

#------------------------------------------------------------------------------------------------------------------

def boom_moi(t_d,number_booms,moment_cs,chord_length):


number_points = number_booms/2 + 2
n = 1000 / number_points
boom_locationx = []
boom_locationy = []
boom_area = []
boom_moi = []
stress_boom_upper = []
stress_boom_lower = []



#make airfoil symetrical and remove negative values at end

for i in range(len(xcoord1[:-3])):
    if int(xcoord1[i] * 1000) % n == 0:
        boom_locationx.append(xcoord1[i])
        boom_locationy.append(ycoord1[i])


#compute the boom area and moment of inertia

for i in range(len(boom_locationx)):
    if i == 1:
        boom_12dx = (boom_locationx[i - 1] - boom_locationx[i]) * chord_length
        boom_12dy = (boom_locationy[i - 1] - boom_locationy[i]) * chord_length
        boom_23dx = (boom_locationx[i + 1] - boom_locationx[i]) * chord_length
        boom_23dy = (boom_locationy[i + 1] - boom_locationy[i]) * chord_length

        b_1 = 2*np.sqrt((boom_12dx) ** 2 + (boom_12dy) ** 2)
        b_2 = np.sqrt((boom_23dx) ** 2 + (boom_23dy) ** 2)

        area_boom = t_d * b_1 / 6 + t_d * b_2 * (2 + boom_locationy[i + 1] / boom_locationy[i]) / 6
        moi_boom = area_boom * (boom_locationy[i] * chord_length) ** 2

        boom_area.append(area_boom)
        boom_moi.append(moi_boom)
    elif i == len(boom_locationx)-1:
        boom_12dx = (boom_locationx[i - 1] - boom_locationx[i]) * chord_length
        boom_12dy = (boom_locationy[i - 1] - boom_locationy[i]) * chord_length
        boom_23dx = (1 - boom_locationx[i]) * chord_length
        boom_23dy = boom_locationy[i] * chord_length

        b_1 = np.sqrt((boom_12dx) ** 2 + (boom_12dy) ** 2)
        b_2 = 2* np.sqrt((boom_23dx) ** 2 + (boom_23dy) ** 2)

        area_boom = t_d * b_1 * (2 + boom_locationy[i - 1] / boom_locationy[i]) / 6 + t_d * b_2 / 6
        moi_boom = area_boom * (boom_locationy[i] * chord_length) ** 2

        boom_area.append(area_boom)
        boom_moi.append(moi_boom)
    if i>1 and i< len(boom_locationx)-1:
        boom_12dx = (boom_locationx[i - 1] - boom_locationx[i]) * chord_length
        boom_12dy = (boom_locationy[i - 1] - boom_locationy[i]) * chord_length
        boom_23dx = (boom_locationx[i + 1] - boom_locationx[i]) * chord_length
        boom_23dy = (boom_locationy[i + 1] - boom_locationy[i]) * chord_length

        b_1 = np.sqrt((boom_12dx)**2+(boom_12dy)**2)
        b_2 = np.sqrt((boom_23dx)**2+(boom_23dy)**2)

        area_boom = t_d*b_1*(2+boom_locationy[i-1]/boom_locationy[i])/6+t_d*b_2*(2+boom_locationy[i+1]/boom_locationy[i])/6
        moi_boom = area_boom*(boom_locationy[i]*chord_length)**2

        boom_area.append(area_boom)
        boom_moi.append(moi_boom)

# total moment of inertia of the structure (2x top and bottom)
moi_boom_total = sum(2*boom_moi)

for i in range(len(boom_locationx)-1):
    if i>0:
        boom_stress = moment_cs*boom_locationy[i]*chord_length / moi_boom_total
        stress_boom_upper.append(boom_stress)
        stress_boom_lower.append(-boom_stress)

print(stress_boom_lower)
# Find the shear flows






# Finally look at torsion

LE_sweep = input.LE_sweep
T_to = input.Tto
y_array = Wingbox_design.generate_spanwise_locations(1000)
x_array = y_array * np.sin(LE_sweep)
lift_array_along_y = Wingbox_design.generate_lift_data_points(y_array)
weight_array_along_y = Wingbox_design.generate_weight_data_points(y_array)


y_lift = Wingbox_design.trapezoidal_integration(y_array, y_array*lift_array_along_y)/Wingbox_design.trapezoidal_integration(y_array, lift_array_along_y)
x_lift = y_lift*sin(LE_sweep)
lift = Wingbox_design.trapezoidal_integration(y_array, lift_array_along_y)
y_weight= Wingbox_design.trapezoidal_integration(y_array, y_array*weight_array_along_y)/Wingbox_design.trapezoidal_integration(y_array, weight_array_along_y)
x_weight = y_weight*sin(LE_sweep)
weight = Wingbox_design.trapezoidal_integration(y_array, weight_array_along_y)
y_engine = Wingbox_design.x_engine_root
x_engine = y_engine*sin(LE_sweep)
engine_weight = Wingbox_design.w_engine
z_engine = input.z_engine


def reaction_moment_about_y(lift, weight, engine_weight, x_lift, x_weight, x_engine):
    M = lift*x_lift - weight*x_weight - engine_weight*x_engine - T_to*-z_engine
    return M

R_z, M_y = Wingbox_design.R_z, reaction_moment_about_y(lift, weight, engine_weight, x_lift, x_weight, x_engine)

def internal_y_bending_moment(x, x_array=x_array, y_array=y_array, lift_array=lift_array_along_y, w_engine=engine_weight,\
    weight_array=weight_array_along_y, x_engine = x_engine, R_z = R_z, M_y = M_y):
    #right hand positive with y forward positive

    x = min(x_array, key=lambda y:abs(y-x))
    x_index = np.where(x_array == x)[0][0]
    x_array = x_array[:x_index]
    y_array = y_array[:x_index]
    lift_array = lift_array[:x_index]
    weight_array = weight_array[:x_index]
    x_lift = Wingbox_design.trapezoidal_integration(x_array, x_array*lift_array)/Wingbox_design.trapezoidal_integration(x_array, lift_array)
    lift = Wingbox_design.trapezoidal_integration(y_array, lift_array)
    x_weight = Wingbox_design.trapezoidal_integration(x_array, x_array*weight_array)/Wingbox_design.trapezoidal_integration(x_array, weight_array)
    weight = Wingbox_design.trapezoidal_integration(y_array, weight_array)
    if x > x_engine:
        engine_distance = x - x_engine
    else:
        engine_distance = 0
    moment_at_x = M_y + R_z*x + lift*(x-x_lift) - w_engine*engine_distance - weight*(x-x_weight)
    return moment_at_x, lift, weight, x_lift, x_weight
M_y_array=[]
for x in x_array[2:]:
    M_y_array.append(internal_y_bending_moment(x)[0])

print(lift, weight, engine_weight, x_lift, x_weight, x_engine, M_y, R_z, x_array[-1])
print(internal_y_bending_moment(x_array[-1]))
<<<<<<< HEAD
plt.plot(x_array[2:], M_y_array )
# plt.show()

=======
plt.plot(y_array[2:], M_y_array )
plt.show()
>>>>>>> 542423d13928d25969d6114ef66bea4502ee60e3
