import input
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

lift_d = SMC*wingloading #not yet adapted to triangle
w_wingd = 1000 # wing weight distribution

def elliptical_lift_d (x, a=10.0135, b=40957.2):
    #obtained by fitting elliptical distribution to given wing loading, a=half wing span - half fuselage width
    loading_at_x = sqrt((1-(x**2/a**2))*b**2)
    return loading_at_x


x_lift = 4.2499 #application point of lift force


def wing_root_reaction_forces (L_wing, x_lift, W_wing, x_weight, W_engine, x_engine):
    R_y = W_wing + W_engine - L_wing  #upwards positive
    M = x_lift*L_wing - x_weight*W_wing - x_engine*W_engine #clockwise positive
    return (R_y, M)


#section 1-2 wingtip to engine




#section 2-3 engine to wing root