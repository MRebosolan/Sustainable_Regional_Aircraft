import input
import matplotlib.pyplot as plt
import numpy as np
from math import *
import Class_2_estimation as cl2


wingloading = input.wingloading #import wingloading

b = cl2.b #import wingspan
S = cl2.S #import wing area
AR = b**2 / S #aspect ratio
taper = input.taper #input taper
t_r = input.t_r # input t_r
SMC = b/AR #standar mean chord
widthf = input.widthf
wing_length = 0.5*(b)
wing_weight = 3000*9.81
T_to = input.Tto

total_lift = wingloading/wing_length
total_lift_half = total_lift/2


lift_duniform = total_lift_half/wing_length
lift_dtriangle = 2*total_lift_half/wing_length


w_wingd = 3000*9.81 # wing weight distribution
w_engine = 1300*9.81 # engine weight 1300kg
x_engine = 2*wing_length/3 # engine distance from tip

#-------UNIFORM LIFT DISTRIBUTION CALCULATIONS------------------

x_loc = []
y_shear_uniform = []
y_moment_uniform = []

for i in range(1,100):
    x = wing_length*(i/100)
    if x<x_engine:
        Shear_Force = (lift_duniform - w_wingd)*x
        Moment = (lift_duniform*x**2 - w_wingd*x**2)/2
    else:
        Shear_Force = (lift_duniform - w_wingd)*x - w_engine
        Moment = (lift_duniform - w_wingd)*x**2/2 - (x-x_engine)*w_engine
    x_loc.append(x)
    y_shear_uniform.append(Shear_Force)
    y_moment_uniform.append(Moment)

#-------TRIANGULAR LIFT DISTRIBUTION CALCULATIONS------------------

y_shear_triangle = []
y_moment_triangle = []

for i in range(1,100):
    x = wing_length*(i/100)
    if x<x_engine:
        Shear_Force = lift_dtriangle*x**2/(2*wing_length) - w_wingd*x
        Moment = lift_dtriangle*x**3/(6*wing_length) - w_wingd*x**2/2
    else:
        Shear_Force = lift_dtriangle*x**2/(2*wing_length) - w_wingd*x - w_engine
        Moment = lift_dtriangle*x**3/(6*wing_length) - w_wingd*x**2/2 - (x-x_engine)*w_engine

    y_shear_triangle.append(Shear_Force)
    y_moment_triangle.append(Moment)


#-------ELLIPTICAL/GENERAL LIFT DISTRIBUTION CALCULATIONS------------------


def generate_spanwise_locations(n, b=wing_length):
    x_array = np.linspace(0, b, n, endpoint=False)
    return x_array[1:]

def general_lift_d (x, a=11.43427, b=15926.5):
    #Insert formula for lift distribution here
    #current obtained by fitting elliptical distribution to given wing loading, a=half wing span - half fuselage width
    loading_at_x = np.sqrt((1-(x**2/a**2))*b**2)
    return loading_at_x

def general_weight_d (x, slope=-450.248, intercept=5148.26):
    #Insert weight distribution. Current is linear for wing weight of 300*9.81
    weight_at_x = intercept + slope*x
    return weight_at_x

def generate_lift_data_points(x_array):
    lift_array = []
    for i in x_array:
        lift_array.append(general_lift_d(i))
    return lift_array

def generate_weight_data_points(x_array):
    weight_array = []
    for i in x_array:
        weight_array.append(general_weight_d(i))
    return weight_array

def trapezoidal_integration(x_array, y_array):
    integral_value = 0
    for i in range(len(x_array)-1):
        integral_value += (x_array[i+1]-x_array[i])/2 * (y_array[i+1] + y_array[i])
    return integral_value

x_array = generate_spanwise_locations(1000)

x_engine_root = wing_length/3
x_weight = wing_length/3
lift_array = generate_lift_data_points(x_array)
x_lift = trapezoidal_integration(x_array, x_array*lift_array)/trapezoidal_integration\
    (x_array, lift_array)
weight_array = generate_weight_data_points(x_array)
L_wing = trapezoidal_integration(x_array, lift_array)
W_wing = trapezoidal_integration(x_array, weight_array)


def wing_root_reaction_forces (L_wing, x_lift, W_wing, x_weight, W_engine, x_engine, T_to):
    #Drag reaction forces not included yet
    R_y = W_wing + W_engine - L_wing  #upwards positive
    M_x = x_lift*L_wing - x_weight*W_wing - x_engine*W_engine #left hand positive
    R_x = T_to #aft-ward positive
    M_z = T_to * x_engine #left hand positive
    return (R_y, M_x, R_x, M_z)

R_y, M_x, R_x, M_z = wing_root_reaction_forces(L_wing, x_lift, wing_weight, x_weight, w_engine, x_engine_root, T_to)


def internal_x_bending_moment(x, x_array=x_array, lift_array=lift_array, w_engine=w_engine,\
    weight_array=weight_array, x_engine = x_engine_root, R_y = R_y, M = M_x):
    #counterclockwise positive

    x = min(x_array, key=lambda y:abs(y-x))
    x_index = np.where(x_array == x)[0][0]
    x_array = x_array[:x_index]
    lift_array = lift_array[:x_index]
    weight_array = weight_array[:x_index]
    x_lift = trapezoidal_integration(x_array, x_array*lift_array)/trapezoidal_integration(x_array, lift_array)
    lift = trapezoidal_integration(x_array, lift_array)
    x_weight = trapezoidal_integration(x_array, x_array*weight_array)/trapezoidal_integration(x_array, weight_array)
    weight = trapezoidal_integration(x_array, weight_array)
    if x > x_engine:
        engine_distance = x - x_engine
    else:
        engine_distance = 0
    moment_at_x = M + R_y*x + lift*(x-x_lift) - w_engine*engine_distance - weight*(x-x_weight)
    return moment_at_x

def internal_z_bending_moment(x, R_x=R_x, T_to=T_to, x_engine=x_engine_root, M_z=M_z):
    if x>x_engine:
        d=x-x_engine
    else:
        d=0
    moment_at_x = M_z -R_x*x + T_to*d
    return moment_at_x


def internal_vertical_shear_force(x, x_array=x_array, lift_array=lift_array, w_engine=w_engine,\
    weight_array=weight_array, x_engine = x_engine_root, R_y = R_y):
    #downward positive

    x = min(x_array, key=lambda y: abs(y - x))
    x_index = np.where(x_array == x)[0][0]
    x_array = x_array[:x_index]
    lift_array = lift_array[:x_index]
    weight_array = weight_array[:x_index]
    lift = trapezoidal_integration(x_array, lift_array)
    weight = trapezoidal_integration(x_array, weight_array)
    if x > x_engine:
        n=1
    else:
        n=0
    shear_at_x = R_y + lift - weight - w_engine * n
    return shear_at_x

def internal_longitudinal_shear_force(x, R_x=R_x, T_to=T_to, x_engine=x_engine_root):
    #forward positive
    if x>x_engine:
        return 0
    else:
        return R_x

moment_array = []
shear_array = []
moment_array2=[]
shear_array2=[]

for i in x_array[2:]:
    moment_array.append(internal_x_bending_moment(i))
    shear_array.append(internal_vertical_shear_force(i))
    moment_array2.append(internal_z_bending_moment(i))
    shear_array2.append(internal_longitudinal_shear_force(i))

plt.plot(x_array[2:], moment_array)
plt.plot(x_array[2:], shear_array)
plt.plot(x_array[2:], moment_array2)
plt.plot(x_array[2:], shear_array2)




#-------BENDING STRESS CALCULATIONS------------------
base_wb = 0.05 # wingbox base
height_wb = 0.05 # wingbox height
thickness_wb = [0.5,1,1.5] # create array of thickness
material_wb = [100,200,300] # insert all material types
y = 0.05  #y in meters
# moment of inertia of a rectangle
bending_moment = 890 # link to functions above feeding an input

for i in thickness_wb:
    area_wb = 2*i*(base_wb+height_wb) #area of wingbox
    volume_wb = area_wb * wing_length  # assume that the wingbox covers the majority of the wing length
    for j in material_wb:
        mass_wb = volume_wb*j

for i in thickness_wb:
    moi_rectangle_x = i*base_wb*height_wb**2/3
    moi_rectangle_y = i*base_wb**2*height_wb/3
    hoop_x = - bending_moment*(height_wb/2)/moi_rectangle_x
    hoop_y = - bending_moment*(base_wb/2)/moi_rectangle_y


plt.show()
