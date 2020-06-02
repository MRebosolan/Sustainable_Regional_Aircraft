import input
import matplotlib.pyplot as plt
import numpy as np
from math import *


wingloading = input.wingloading #import wingloading

b = input.b #import wingspan
S = input.S #import wing area
AR = input.AR #aspect ratio
taper = input.taper #input taper
t_r = input.t_r # input t_r
SMC = b/AR #standar mean chord
widthf = input.widthf
wing_length = 0.5*(b-widthf)

total_lift = wingloading/wing_length
total_lift_half = total_lift/2

lift_duniform = total_lift_half/wing_length
lift_dtriangle = 2*total_lift_half/wing_length


w_wingd = 3000 # wing weight distribution
w_engine = 1300 # engine weight 1300kg
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
        Moment = (lift_duniform*x**2 - w_wingd*x**2)/2 - (x-x_engine)*w_engine
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


#-------ELLIPTICAL LIFT DISTRIBUTION CALCULATIONS------------------

def elliptical_lift_d (x, a=10.0135, b=40957.2):
    #obtained by fitting elliptical distribution to given wing loading, a=half wing span - half fuselage width
    loading_at_x = sqrt((1-(x**2/a**2))*b**2)
    return loading_at_x

def generate_spanwise_locations(n, b=wing_length):
    x_array = [0]
    step = b/n
    for i in range (n):
        x_array.append(x_array[-1]+step)
    return np.array(x_array)

def generate_lift_data_points(x_array):
    lift_array = []
    for i in x_array:
        lift_array.append(elliptical_lift_d(i))
    return lift_array

def trapezoidal_integration(x_array, y_array):
    integral_value = 0
    for i in range (len(x_array)-1):
        integral_value += (x_array[i+1]-x_array[i])/2 * (y_array[i+1] + y_array[i])
    return integral_value



x_lift = 4.2499 #application point of lift force


def wing_root_reaction_forces (L_wing, x_lift, W_wing, x_weight, W_engine, x_engine):
    R_y = W_wing + W_engine - L_wing  #upwards positive
    M = x_lift*L_wing - x_weight*W_wing - x_engine*W_engine #clockwise positive
    return (R_y, M)


#section 1-2 wingtip to engine

plt.plot(x_loc, y_moment_uniform)
plt.plot(x_loc, y_shear_uniform)
plt.plot(x_loc, y_moment_triangle)
plt.plot(x_loc, y_shear_triangle)

plt.show()
