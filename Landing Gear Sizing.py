import numpy as np
import input
"""
Responsible person: Tobias

This code requires as inputs:
    MTOW           | Maximum take-off weight
    g              | gravitational acceleration
    theta          | Clearance angle 
    x_cg           | X location of the center of gravity
    z_cg           | Z location of the center of gravity
    x_cg_fwrd      | Most forward location in x direction of cg
    x_cg_aft       | Most aft location in x direction of cg
    x_lg           | X location of main landing gear

This code gives as outputs:
    Angle beta     | The line through the contactpoint of the wheel and ground that goes through the c.g. should make an angle of beta > theta with the line normal to the wheel-ground contact point. It is required to be strictly larger
    
    


Gear strut must be fully extended: 
    The line through the contactpoint of the wheel and ground that goes
    through the c.g. should make an angle of beta > theta with the line normal
    to the wheel-ground contact point. It is required to be strictly larger
Note: When determining the position of the main landing gear, it is important
    to keep in mind that the more aft the landing gear is, the larger down
    force is required to be produces by the horizontal tail surface during the
    rotation phase for take-off

"""
#Variables that will not change
MTOW = input.MTOW
g = input.g
theta = input.theta #tip-back angle ~15 degrees
z_cg  =  input.z_cg
x_cg  = input.x_cg
x_cg_fwrd = input.x_cg_fwrd
x_cg_aft = input.x_cg_aft

#Variables that can/will be changed
x_lg =   #Obtained from weight &balance when looking at c.g. excursion

def beta_check(theta,x_lg,x_cg,z_cg,x_cg_fwrd,x_cg_aft):
    beta = np.arccos(z_cg/np.sqrt((x_lg-x_cg)**2+z_cg**2))
    beta_fwrd = np.arccos(z_cg/np.sqrt((x_lg-x_cg_fwrd)**2+z_cg**2))
    beta_aft = np.arccos(z_cg/np.sqrt((x_lg-x_cg_aft)**2+z_cg**2)) 
    return beta, beta_fwrd,beta_aft

beta, beta_fwrd,beta_af = beta_check(theta,x_lg,x_cg,z_cg,x_cg_fwrd,x_cg_aft)


if beta <= theta:
    print ('Beta is not larger than theta; the main landing gear must be moved further back')
elif beta_fwrd <= theta:
    print ('Beta_fwrd is not larger than theta; the main landing gear must be moved further back')

elif beta_aft <= theta:
    print ('Beta_aft is not larger than theta; the main landing gear must be moved further back')
else:
    print ('the x position of the landing gear is large enough')

#weight distribution check:

def gear_loading_check(MTOW,g):
    W = 2*F_M+F_N
    