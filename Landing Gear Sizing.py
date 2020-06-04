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

This code gives as outputs:
    The main landing gear position
    The minimum and maximum distance of the nose landing gear w.r.t. 
        the nose 
    The lateral position (y_lg) of the landing gear taking into account
        the turnover angle (has to be less than 55 degrees)
    
    
Remark(s):
    For now dummy variables have been used, as it is not yet certain 
    which programs provide which parameters.

"""
#Variables that will not change
MTOW = Cl2.MTOW
g = input.g
theta = (np.pi/180)*input.theta     #tip-back angle ~15 degrees
z_cg  =  input.z_cg
x_cg  = input.x_cg
x_cg_fwrd = input.x_cg_fwrd
x_cg_aft = input.x_cg_aft

#Determining the minimum required surface area
S = input.S
rho_0 = input.rho0
rho_to = rho_0
CLmax = input.CLmax_land 
Vmin = np.sqrt(MTOW*g * 2 /(S * rho_to * CLmax))
Vlof = 1.05*Vmin
Cl_htail = 0.0 #tbd
x_ac_htail = 0.0 # distance from aerodynamic centre to nose of htail airfoil

######## Dummy variables to test the program, as some have not yet
# Been determined ################################################

import numpy as np
safetymargin_theta = 1*(np.pi/180)
theta = 15*np.pi/180
z_cg = 3.0
x_cg_aft = 8
x_cg = x_cg_aft-0.2
MTOW = 35000
g = 9.81

rho_to = 1.225
S = 80
CLmax = 2.2
Cl_htail = 0.5
x_ac_htail = 17
Vmin = np.sqrt(MTOW*g*2/(S*rho_to*CLmax))
Vlof = Vmin*1.05
#################################################

def main_lg_loc(theta=theta,z_cg=z_cg,x_cg_aft=x_cg_aft,safetymargin_theta=safetymargin_theta):
    x_main_lg = x_cg_aft + z_cg*np.tan(theta+safetymargin_theta)
    return x_main_lg
x_main_lg = main_lg_loc()
#print (x_main_lg)

def nose_lg_loc(x_main_lg= x_main_lg, x_cg=x_cg,MTOW=MTOW,g=g):
    dist = []
    d = 0.005
    for distance in np.arange(0+d,15,d):
        F_nose_lg = MTOW*g*(x_main_lg-x_cg)/(x_cg-distance)
        #Force_on_nose_lg.append(F_nose_lg)
        if F_nose_lg <= 0.08*MTOW*g or F_nose_lg >= 0.15*MTOW*g:
            continue
        else:
            dist.append(distance)
    return dist
dist = nose_lg_loc()

#print('The x-location of the main landing gear w.r.t. the nose is', np.round(x_main_lg,3),'m')
#print('Nose gear: The minimum x-distance from the nose equals', np.round(np.min(dist),3),'m')
#print('Nose gear: The maximum x-distance from the nose equals', np.round(np.max(dist),3),'m')

dist_max = np.round(np.max(dist),3) #the minimum distance between the nose and main landing gear
dist_min = np.round(np.min(dist),3) #the maximum distance between the nose and main landing gear

#dist_lg =  # this will be the actual value of the distance, above values are used to model the ranges of lateral positions of the main landing gear

"""
Psi, the tipover angle must be smaller than 55 degrees
"""

def lat_pos_lg(z_cg,dist_min=dist_min,dist_max=dist_max,x_main_lg=x_main_lg):
    y_lg_list = []
    d = 0.005
    for y_lg in np.arange(d,4+d,d):
        lg_distance = [x_main_lg-dist_max,x_main_lg-dist_min] #distance between the nose and main landing gear
        for i in range(len(lg_distance)):      
            alpha = np.arctan2(y_lg,lg_distance[i])
            c = (x_cg-lg_distance[i])*np.sin(alpha)
            psi = np.arctan2(z_cg,c)
            if psi >= 55/180*np.pi:
                continue
            else:
                y_lg_list.append(y_lg)  
    return y_lg_list, lg_distance

y_lg_list, lg_distance = lat_pos_lg(z_cg)

def req_htail_area(x_main_lg,Cl_htail=Cl_htail,x_ac_htail=x_ac_htail,x_cg = x_cg,rho_to=rho_to,Vlof=Vlof,MTOW=MTOW,g=g):
    htail_area = (x_main_lg-x_cg)/(x_ac_htail-x_main_lg)*MTOW*g/(0.5*rho_to*Vlof**2*Cl_htail)
    return htail_area
htail_area = req_htail_area(x_main_lg)
print ('Required htail area:', np.round(htail_area,3))
print ('Lateral distance of the landing gear:',y_lg_list)
