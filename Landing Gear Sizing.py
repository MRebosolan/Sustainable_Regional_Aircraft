import numpy as np
import input
import Class_2_estimation as Cl2
"""
Responsible person: Tobias | FOR NOW!!!! run from line 55

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
g = input.g
MTOW = g*Cl2.MTOM
theta = (np.pi/180)*input.theta     #tip-back angle ~15 degrees
z_cg  =  0.5*input.hf                  #z loc cg. WRT fuselage!
x_cg  = input.x_cg                  #x location of the cg
x_cg_fwrd = input.x_cg_fwrd         #x location of most forward cg
x_cg_aft = input.x_cg_aft           #x-location of most aft cg
x_empennage = input.x_empennage     #x-location where fuselage diameter 
#starts decreasing (clearance angle), wrt c.g. OR wrt nose
#If it is wrt the the c.g. change x_empennage

#Determining the minimum required surface area
S = Cl2.S
rho_0 = input.rho0
CLmax = input.CLmax_land 
Vmin = np.sqrt(MTOW*g * 2 /(S * rho_0 * CLmax))
Vlof = 1.05*Vmin
Cl_htail = 0.0 #tbd
x_ac_htail = 0.0 # distance from aerodynamic centre to nose of htail airfoil

safetymargin_theta = 1*(np.pi/180)
rho_to = rho_0
######## Dummy variables to test the program, as some have not yet
# Been determined ################################################

import input
import numpy as np
safetymargin_theta = 1*(np.pi/180)
theta = 15*np.pi/180
z_cg = 0.5*(input.hf)   #z loc cg. WRT fuselage!
x_cg_aft = 10
x_cg = x_cg_aft-0.2
MTOW = 26000
g = 9.81
x_empennage = 17
rho_to = 1.225
S = 80
CLmax = 2.2
Cl_htail = 1.5
x_ac_htail = 22
Vmin = np.sqrt(MTOW*g*2/(S*rho_to*CLmax))
Vlof = Vmin*1.05

#print (Vmin,Vlof)
#################################################

def main_lg_loc(x_empennage=x_empennage,theta=theta,z_cg=z_cg,x_cg_aft=x_cg_aft,safetymargin_theta=safetymargin_theta):
    d = 0.01
    for z_fus_ground in np.arange(0,5,d):
        
        z_main_lg1 = z_fus_ground + z_cg
        x_main_lg1 = x_cg_aft + z_main_lg1*np.tan(theta+safetymargin_theta)
        if np.tan(z_fus_ground/(x_empennage-x_main_lg1)) >= 15*np.pi/180:
            z_main_lg = np.round(z_main_lg1,4)
            x_main_lg = np.round(x_main_lg1,4)
            z_f_ground = z_fus_ground
            break
        else:
            continue
    print ()
    print ('[Make sure that x_empennage is wrt the nose!]')
    return x_main_lg,z_main_lg, z_f_ground    
x_main_lg, z_main_lg, z_f_ground  = main_lg_loc()

print ()
print('The x-location of the main landing gear w.r.t. the nose is', np.round(x_main_lg,4),'m')
print ()
print ('The distance from the ground to c.g. equals:', z_main_lg,'[m]')
print ()
print ('The Clearance angle is:',np.round(np.tan((z_f_ground)/(x_empennage-x_main_lg))*180/np.pi,3),'[deg]')
print ()

def nose_lg_loc(x_main_lg= x_main_lg, x_cg=x_cg,MTOW=MTOW,g=g):
    dist = []
    d = 0.005
    for distance in np.arange(-10+d,x_cg,d):
        F_nose_lg = MTOW*g*(x_main_lg-x_cg)/(x_cg-distance)
        #Force_on_nose_lg.append(F_nose_lg)
        if F_nose_lg <= 0.08*MTOW*g or F_nose_lg >= 0.15*MTOW*g:
            continue
        else:
            dist.append(distance)
    return dist
dist = nose_lg_loc()

print('Nose gear: The minimum x-distance from the nose equals', np.round(np.min(dist),4),'[m]')
print ()
print ('In case the min distance (or both) is negative, it is not an error. It simply calculates the lower and upper limits to comply with the requirements')
print()
print('Nose gear: The maximum x-distance from the nose equals', np.round(np.max(dist),4),'[m]')
print ()

#dist_max = np.round(np.max(dist),3) #the minimum distance between the nose and nose landing gear
#dist_min = np.round(np.min(dist),3) #the maximum distance between the nose and nose landing gear
#dist_lg =  # this will be the actual value of the distance, above values are used to model the ranges of lateral positions of the main landing gear

"""
Psi, the tipover angle must be smaller than 55 degrees
"""


def lat_pos_lg(z_main_lg=z_main_lg,dist=dist,x_main_lg=x_main_lg,x_cg_aft=x_cg_aft):
    y_lg_list = []
    d = 0.005
    for y_lg in np.arange(d,4+d,d):
        b_n = [x_cg_aft-np.min(dist),x_cg_aft-np.max(dist)]   #distance from most forward nose lg to aft cg.
        for i in range(len(b_n)):
            alpha = np.arctan2(y_lg,b_n[i])
            c = b_n[i]*np.sin(alpha)
            psi = np.arctan2(z_main_lg,c)
            if psi < 55/180*np.pi and b_n[i] <= x_main_lg:
                y_lg_list.append(y_lg)
            else:
                continue  
    return y_lg_list, b_n

y_lg_list, b_n = lat_pos_lg(z_main_lg)

def req_htail_area(x_main_lg,Cl_htail=Cl_htail,x_ac_htail=x_ac_htail,x_cg = x_cg,rho_to=rho_to,Vlof=Vlof,MTOW=MTOW,g=g):
    htail_area = (x_main_lg-x_cg)/(x_ac_htail-x_main_lg)*MTOW*g/(0.5*rho_to*Vlof**2*Cl_htail)
    return htail_area
htail_area = req_htail_area(x_main_lg)

print ('The minimum lateral distance of the landing gear:',np.round(min(y_lg_list),3),'[m]')
print ('This means the main landing gear stick out',np.round(min(y_lg_list)-z_cg,3),'meters from the fuselage' )
print ()
print ('Required htail area:', np.round(htail_area,3),'[m^2]')

"""
Strength of landing gear: take the shock of landing into account
Size for strength
Size tires
Size strut?
"""