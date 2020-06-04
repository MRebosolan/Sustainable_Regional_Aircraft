import numpy as np
import input
import Class_2_estimation as Cl2
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
MTOW = Cl2.MTOW
g = input.g
theta = (np.pi/180)*input.theta     #tip-back angle ~15 degrees
z_cg  =  input.z_cg
x_cg  = input.x_cg
x_cg_fwrd = input.x_cg_fwrd
x_cg_aft = input.x_cg_aft

######## Dummy variables to test the program, as some have not yet
# Been determined ################################################

import numpy as np
safetymargin_theta = 1*(np.pi/180)
theta = 15*np.pi/180
z_cg = 4
x_cg_aft = 15
x_cg = x_cg_aft - np.tan(theta)*z_cg
MTOW = 35000
g = 9.81
#################################################

def main_lg_loc(theta=theta,z_cg=z_cg,x_cg_aft=x_cg_aft,safetymargin_theta=safetymargin_theta):
    x_main_lg = x_cg_aft + z_cg*np.tan(theta+safetymargin_theta)
    return x_main_lg
x_main_lg = main_lg_loc()


def nose_lg_loc(x_main_lg= x_main_lg, x_cg=x_cg,MTOW=MTOW,g=g):
    Force_on_nose_lg = []
    d = 0.005
    for distance in np.arange(0+d,15,d):
        F_nose_lg = MTOW*g*(x_main_lg-x_cg)/distance
        #Force_on_nose_lg.append(F_nose_lg)
        if F_nose_lg <= 0.08*MTOW*g or F_nose_lg >= 0.15*MTOW*g:
            continue
        else:
            Force_on_nose_lg.append(F_nose_lg)
    return Force_on_nose_lg
Force_on_nose_lg = nose_lg_loc()


print('The x-location of the main landing gear w.r.t. the nose is',x_main_lg,)
print('The minimum x-distance from the nose equals', np.min(Force_on_nose_lg))
print('The maximum x-distance from the nose equals', np.max(Force_on_nose_lg))