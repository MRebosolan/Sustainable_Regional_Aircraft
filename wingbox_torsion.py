# Wingbox for bending stress
import numpy as np
import input
import Wingbox_design
import matplotlib.pyplot as plt
from math import *
#the CIA murdered JFK
# coordinates for the airfoil


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
# plt.show()
#------------------------------------------------------------------------------------------------------------------
#INPUTS

chord_length = 2  # chord length in meters
<<<<<<< HEAD
t_d = 0.01 #THICkNESS OF AIRFOIL
number_booms = 10 #NUMBER OF POINTS (10 POINTS = 16 BOOMS,4p=4b 5p=6b )#  (ON TOP SIDE FOR NOW)
=======
t_d = .02 #THICkNESS OF AIRFOIL
number_booms = 20 #NUMBER OF POINTS (10 POINTS = 16 BOOMS,4p=4b 5p=6b )#  (ON TOP SIDE FOR NOW)
>>>>>>> d12cbb86a0c909806e7021594cb7780993542bf9
moment_cs = 450000 #MOMENT OF CROSS SECTION
Kc, Ks = 6.4, 8.2 #clamping coeeficients
E = 71.7E9

#------------------------------------------------------------------------------------------------------------------

def boom_moi(moment_cs, chord_length, shear_cs, t_d=t_d, number_booms=number_booms, stringer_area = 0):

#returns moments of inertia and boom normal stresses based on load, chord length and number of booms

    number_points = number_booms
    boom_area = []
    boom_moi = []
    stress_boom_upper = []
    stress_boom_lower = []
    boom_deltashear = []

    #make airfoil symmetrical and remove negative values at end


    boom_locationx = np.linspace(0, xcoord1[-4], int(number_points))
    boom_locationy = []

    for x in boom_locationx:
        x = min(xcoord1, key=lambda z: abs(z - x))
        x_index = xcoord1.index(x)
        y1 = ycoord1[x_index]
        boom_locationy.append(y1)


    #compute the boom area and moment of inertia

    for i in range(len(boom_locationx)):
        if i == 1:
            boom_12dx = (boom_locationx[i - 1] - boom_locationx[i]) * chord_length
            boom_12dy = (boom_locationy[i - 1] - boom_locationy[i]) * chord_length
            boom_23dx = (boom_locationx[i + 1] - boom_locationx[i]) * chord_length
            boom_23dy = (boom_locationy[i + 1] - boom_locationy[i]) * chord_length

            b_1 = 2*np.sqrt((boom_12dx) ** 2 + (boom_12dy) ** 2)
            b_2 = np.sqrt((boom_23dx) ** 2 + (boom_23dy) ** 2)

            area_boom = stringer_area + t_d * b_1 / 6 + t_d * b_2 * (2 + boom_locationy[i + 1] / boom_locationy[i]) / 6
            moi_boom = area_boom * (boom_locationy[i] * chord_length) ** 2

            boom_area.append(area_boom)
            boom_moi.append(moi_boom)
        if i == len(boom_locationx)-1:
            boom_12dx = (boom_locationx[i - 1] - boom_locationx[i]) * chord_length
            boom_12dy = (boom_locationy[i - 1] - boom_locationy[i]) * chord_length
            boom_23dx = (1 - boom_locationx[i]) * chord_length
            boom_23dy = boom_locationy[i] * chord_length

            b_1 = np.sqrt((boom_12dx) ** 2 + (boom_12dy) ** 2)
            b_2 = 2* np.sqrt((boom_23dx) ** 2 + (boom_23dy) ** 2)

            area_boom = stringer_area + t_d * b_1 * (2 + boom_locationy[i - 1] / boom_locationy[i]) / 6 + t_d * b_2 / 6
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

            area_boom = stringer_area + t_d*b_1*(2+boom_locationy[i-1]/boom_locationy[i])/6+t_d*b_2*(2+boom_locationy[i+1]/boom_locationy[i])/6
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

    for i in range(len(boom_locationx)-1):
        if i > 0:
            deltashear_boom = shear_cs * boom_area[i] * boom_locationy[i] / moi_boom_total
            boom_deltashear.append(deltashear_boom)

    print(boom_deltashear)
    shear_flow1 = []
    shear_flow2 = []

    shearflow = 0
    shear_center = 0.4
    shear_flow1.append(shearflow)

    for i in range(len(boom_deltashear)):
<<<<<<< HEAD
        if boom_locationx[i] > shear_center:
            shearflow = shearflow + boom_deltashear[i]
            shear_flow1.append(shearflow)
        if boom_locationx[i] < shear_center:
            shearflow = shearflow + boom_deltashear[i]
            shear_flow2.append(shearflow)

    shear_flow = shear_flow2 + shear_flow1
    shear_stress_upper = np.array(shear_flow) / t_d
    shear_stress_lower = -shear_stress_upper
=======
        if boom_locationx[i] >= 1/3:
            shearflow = shearflow + boom_deltashear[i]
            shear_flow1.append(shearflow)
    shearflow = 0
    for i in range(len(boom_deltashear)):
        if boom_locationx[i] < 1 / 3:
            shearflow = shearflow - boom_deltashear[i]
            shear_flow2.append(shearflow)

    shear_flow = shear_flow2[::-1] + shear_flow1
    shear_stress_upper = np.array(shear_flow) / t_d
    #shear_stress_lower = -shear_stress_upper

    #stringer pitch
>>>>>>> d12cbb86a0c909806e7021594cb7780993542bf9

    stringer_pitch = (boom_locationx[2]-boom_locationx[1])*chord_length

    return moi_boom_total, max(stress_boom_upper), min(stress_boom_lower), shear_stress_upper[0], shear_stress_upper[-1], stringer_pitch


#-----------------ITERATE OVER WINGSPAN-----------------

Cr = input.Cr
Ct = input.Ct
wing_length = input.b * 0.5


def generate_chord_array(y_array, Cr=Cr, Ct=Ct, b=wing_length):
    chords = []
    def generate_chord_lengths(y, Cr=Cr, Ct=Ct, b=wing_length):
        c = Cr - y*(Cr-Ct)/b
        return c
    for y in y_array:
        chords.append(generate_chord_lengths(y))
    return chords


spanwise_array = Wingbox_design.generate_spanwise_locations(1000)[2:]
moments_around_x = []
moments_around_z = []
shears_y = []
chords = generate_chord_array(spanwise_array)
moi_boom_along_span=[]
upper_stress_along_span=[]
lower_stress_along_span = []
shear_stresses = []
pitches = []

for y in spanwise_array:
    moments_around_x.append(Wingbox_design.internal_x_bending_moment(y))
    moments_around_z.append(Wingbox_design.internal_z_bending_moment(y))
    shears_y.append(Wingbox_design.internal_vertical_shear_force(y))


for i, y in enumerate(spanwise_array):
    moment_around_x = moments_around_x[i-1]
    shear_y = shears_y[i-1]
    chord = chords[i]
    moi_boom, stress_boom_upper, stress_boom_lower, shear_stress_upper, shear_stress_lower, pitch = \
        boom_moi(moment_around_x, chord, shear_y)
    moi_boom_along_span.append(moi_boom)
    upper_stress_along_span.append(stress_boom_upper)
    lower_stress_along_span.append(stress_boom_lower)
    shear_stress = max(abs(shear_stress_upper), abs(shear_stress_lower))
    shear_stresses.append(shear_stress)
    pitches.append(pitch)

# plt.plot(spanwise_array, moments_around_x)
# plt.plot(spanwise_array, shears_y)

# plt.plot(spanwise_array, moi_boom_along_span)
plt.plot(spanwise_array, upper_stress_along_span)
plt.plot(spanwise_array, lower_stress_along_span)
plt.plot(spanwise_array, shear_stresses)
#plt.plot(spanwise_array, pitches)
#plt.show()

print(chords[0] * max(ycoord1))

# compression_buckling_stress = Kc * E * (t_d/max(pitches))**2
# shear_buckling_stress = Ks * E * (t_d/max(pitches))**2
# print(compression_buckling_stress, shear_buckling_stress)
# print(pitches)

#-----------------------------------------------------------------------------------------------------------

# Find the shear flows
# boom_deltashear = []
# for i in range(len(boom_locationx)):
#     deltashear_boom = shear_cs*boom_area[i]*boom_locationy[i]/moi_boom_total
#     boom_deltashear.append(deltashear_boom)
#
# # print(boom_deltashear)
#
#
# shear_flow1 = []
# shear_flow2 = []
# shear_stress = []
#
# shearflow = 0
# shear_flow1.append(shearflow)
#
# for i in range(len(boom_locationx[int(number_boom/4):])):
#     shearflow = shearflow + boom_deltashear[int(number_boom/4)+i]
#     shear_flow1.append(shearflow)
#
# for i in range(len(boom_locationx[:int(number_boom/4)])-1):
#     shearflow = shearflow + boom_deltashear[i]
#     shear_flow2.append(shearflow)
#
# shear_flow = shear_flow2+shear_flow1
# shear_stress = np.array(shear_flow)/t_f
#

#-----------------------------------------------------------------------------------------------------------


# Finally look at torsion

#M_y_array is torque (moment about y) and can be expressed as function of y (wingspan) or x (longitudinal distance
#from center of wing root)

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


def reaction_moment_about_y(lift, weight, engine_weight, x_lift, x_weight, x_engine, z_engine):
    M = lift*x_lift - weight*x_weight - engine_weight*x_engine + T_to*z_engine
    return M

R_z, M_y = Wingbox_design.R_z, reaction_moment_about_y(lift, weight, engine_weight, x_lift, x_weight, x_engine, z_engine)

def internal_y_bending_moment(x, x_array=x_array, y_array=y_array, lift_array=lift_array_along_y, w_engine=engine_weight,\
    weight_array=weight_array_along_y, x_engine = x_engine, R_z = R_z, M_y = M_y, T_to=T_to, z_engine=z_engine):
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
    moment_at_x = M_y + R_z*x + lift*(x-x_lift) - w_engine*engine_distance - weight*(x-x_weight) - T_to*z_engine
    return moment_at_x, lift, weight, x_lift, x_weight

M_y_array=[]
for x in x_array[2:]:
    M_y_array.append(internal_y_bending_moment(x)[0])

<<<<<<< HEAD
# print(lift, weight, engine_weight, x_lift, x_weight, x_engine, M_y, R_z, x_array[-1])
# print(internal_y_bending_moment(x_array[-1]))
=======

>>>>>>> d12cbb86a0c909806e7021594cb7780993542bf9
#plt.plot(y_array[2:], M_y_array )
#plt.show()

