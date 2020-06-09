import numpy as np
import input
import Class_2_estimation as Cl2
import scissor_plot_wing_shift as sc_shift
from CG_excursion_wing_shift import tailcone_length
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
#Variables that still need to be properly coupled to other code:
z_cg  =  0.5*input.hf               # z loc cg. WRT fuselage!
#tailcone_length    # x-location where fuselage diameter 
x_tailcone = input.lf - tailcone_length                                    # starts decreasing (clearance angle), wrt c.g. OR wrt nose
                                    # If it is wrt the the c.g. change x_empennage; [Start of the aft cone]
Cl_htail = input.cl_htail_max                     # tbd

#Variables that will not change
g = input.g
MTOW = g*Cl2.MTOM
theta = np.radians(input.theta)                 # tip-back angle ~15 degrees
x_cg  = sc_shift.cg_loaded_nose     # x location of the cg
x_cg_fwrd = sc_shift.cg_fwd*sc_shift.MAC + sc_shift.xlemac        # x location of most forward cg
x_cg_aft = sc_shift.cg_aft*sc_shift.MAC + sc_shift.xlemac      # x-location of most aft cg
x_ac_htail = sc_shift.x_ac_h_nose   # distance from aerodynamic centre to nose of htail airfoil
print (x_ac_htail,x_cg_fwrd)
S = Cl2.S
rho_0 = input.rho0
rho_to = rho_0
CLmax = input.CLmax_land 
Vmin = np.sqrt(MTOW*g * 2 /(S * rho_to * CLmax))
Vlof = 1.05*Vmin
safetymargin_theta = np.radians(1)
htail_sweep = input.half_chord_sweep_hor      # Sweep of the horizontal tail

# ######## Dummy variables to test the program, as some have not yet
# # Been determined ################################################
# import Class_2_estimation as Cl2
# import input
# import numpy as np
# #import scissor_plot_wing_shift as scissor_w_shift


# safetymargin_theta = np.radians(1)
# theta = np.radians(15)
# z_cg = 0.5*(input.hf)   #z loc cg. WRT fuselage!
# x_cg_aft = 10
# x_cg = x_cg_aft-0.2
# g = input.g
# MTOW = g*Cl2.MTOM

# rho_to = 1.225
# S = 80
# CLmax = 2.2
# Cl_htail = 1.5
# x_ac_htail = 22
# Vmin = np.sqrt(MTOW*2/(S*rho_to*CLmax))
# Vlof = Vmin*1.05
# htail_sweep = input.half_chord_sweep_hor      # Sweep of the horizontal tail

def main_lg_loc(x_tailcone=x_tailcone,theta=theta,z_cg=z_cg,x_cg_aft=x_cg_aft,safetymargin_theta=safetymargin_theta):
    d = 0.01
    for z_fus_ground in np.arange(0,5,d):
        
        z_main_lg1 = z_fus_ground + z_cg
        x_main_lg1 = x_cg_aft + z_main_lg1*np.tan(theta+safetymargin_theta)
        if np.tan(z_fus_ground/(x_tailcone-x_main_lg1)) >= np.radians(15):
            z_main_lg = np.round(z_main_lg1,4)
            x_main_lg = np.round(x_main_lg1,4)
            z_f_ground = z_fus_ground
            break
        else:
            continue

    return x_main_lg,z_main_lg, z_f_ground    
x_main_lg, z_main_lg, z_f_ground  = main_lg_loc()

print ()
print('The x-location of the main landing gear w.r.t. the nose is', np.round(x_main_lg,4),'m')
print ()
print ('The distance from the ground to c.g. equals:', z_main_lg,'[m]')
print ()
print ('The Clearance angle is:',np.round(np.tan((z_f_ground)/(x_tailcone-x_main_lg))*180/np.pi,3),'[deg]')
print ()

def nose_lg_loc(x_main_lg= x_main_lg, x_cg=x_cg,MTOW=MTOW,g=g):
    dist = []
    d = 0.001
    for distance in np.arange(0,x_cg,d):
        F_nose_lg = MTOW*(x_main_lg-x_cg)/(x_cg-distance)
        #Force_on_nose_lg.append(F_nose_lg)
        if F_nose_lg <= 0.08*MTOW or F_nose_lg >= 0.15*MTOW:
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

dist_max = np.round(np.max(dist),4) #-1 if it is too much forward; it does not have a large effect on the lateral position; that is to say, the effect is not severe.   #the minimum distance between the nose and nose landing gear
dist_min = np.round(np.min(dist),3) #the maximum distance between the nose and nose landing gear
#dist_lg =  # this will be the actual value of the distance, above values are used to model the ranges of lateral positions of the main landing gear

def lat_pos_lg(z_main_lg=z_main_lg,dist=dist,x_main_lg=x_main_lg,x_cg_aft=x_cg_aft):
    y_lg_list = []
    d = 0.005
    for y_lg in np.arange(d,4+d,d):
        b_n = [x_cg_aft-dist_min,x_cg_aft-dist_max]   #distance from most forward nose lg to aft cg.
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

def req_htail_area(x_main_lg,Cl_htail=Cl_htail,x_ac_htail=x_ac_htail,x_cg = x_cg,rho_to=rho_to,Vlof=Vlof,MTOW=MTOW,g=g,htail_sweep=htail_sweep):
    
    htail_area = -(x_main_lg-x_cg)/(x_ac_htail-x_main_lg)*MTOW/(0.5*rho_to*(Vlof/np.cos(htail_sweep))**2*Cl_htail)
    return htail_area
htail_area = req_htail_area(x_main_lg)

print ('The minimum lateral distance of the landing gear:',np.round(min(y_lg_list),3),'[m]')
print ('This means the main landing gear stick out',np.round(min(y_lg_list)-z_cg,3),'meters from the fuselage' )
print ()
#if np.round(htail_area,3) < scissor_w_shift.min_Sh_over_S:
    #print('The htail surface area is large enough for take-off rotation')
#else:
    #print('Error! The htail surface area is too small for take-off rotation')


N_w = 2       # Number of wheels per strut
N_nw = 2      # Number of nosewheels
N_mw = 2      # Number of main landing wheels | Rounded to nearest multiple of 2, Wto/120000
N_struts = 2  # Number of struts for MAIN lg, if Number of wheels <= 12, 2 struts
LCN = 40      # Load classification number, the ACN of CRJ700 = 23 max value for rigid pavement; this value should be taken as small as possible.

def tire_pressure(LCN=LCN):
    ptire_max = 430*np.log(LCN)-680   # maximum allowable tire pressure
    return ptire_max
ptire_max = tire_pressure()
print (np.round(ptire_max,4),'kPa',np.round(ptire_max,4)*0.145037738,'psi')
#print ('This is actually quite comparable to the CRJ700, according to their airport planning manual')
mg_x_cg = x_main_lg-x_cg # distance from the main lg to the cg
ng_x_cg = x_cg-np.round(np.max(dist),4) # distance from the nose_lg to the cg
mw_nw_d = x_main_lg-np.round(np.max(dist),4)

def static_loads_lg(MTOW=MTOW,N_mw=N_mw,N_struts=N_struts,mg_x_cg=mg_x_cg,ng_x_cg=ng_x_cg,z_cg=z_cg,dist_max=dist_max,mw_nw_d=mw_nw_d):
    P_mw = (MTOW*ng_x_cg)/(N_struts*mw_nw_d) # [N]
    P_nw = MTOW-P_mw*N_struts                # [N]
    ESWL_n = P_nw/1.33/g                     # [kg] equivalent single wheel load twin dual | NOSE
    ESWL_m = P_mw/1.33/g                     # [kg] equivalent single wheel load twin dual | MAIN
    P_mw_stat = 1.07*P_mw/2                  # [N] maximum static load on main gear per wheel (2 wheels 2 struts = 4 wheels)
    P_nw_des = 1.07*P_nw/2                   # [N] maximum static load on nose gear per wheel (2 wheels 1 strut = 2 wheels)
    ax_over_g = 0.35                         # 0.45 for advanced anti-sid brakes
    P_n_dyn_t = MTOW*(mg_x_cg+ax_over_g/z_cg)/((N_nw)*(mw_nw_d))  # [N] Maximum dynamic load per nose gear tire
    P_nw_des_stat = P_n_dyn_t/1.5            #[lbs]
    V_tire_max_to = 1.1*Vlof
    V_bar = 1.3*Vlof/1.05/np.sqrt(2)
    V_tire_max_land = 1.2*V_bar   
    return P_mw,P_nw, ESWL_n, ESWL_m,P_nw_des, P_n_dyn_t,P_nw_des_stat,P_mw_stat,V_tire_max_to, V_tire_max_land
P_mw,P_nw, ESWL_n, ESWL_m,P_nw_des, P_n_dyn_t,P_nw_des_stat,P_mw_stat,V_tire_max_to,V_tire_max_land = static_loads_lg()

print ()
print ('The max. tire speed is',max(V_tire_max_to,V_tire_max_land)*2.2369,'MPH')
print ()
print (1.07*P_mw_stat*0.2248,1.07*max(P_nw_des,P_nw_des_stat)*0.2248,1.07*P_n_dyn_t*0.2248,'Maximum static load on mg, maximum static load on ng and maximum dynamic load ng, respectively in [lbs]')
print ()
print ('ESWL nose and ESWL for main, respectively:', ESWL_n,ESWL_m,'kg')

# DUMMY
D_o = 26         # [in]outside tire diameter, also: Dt
load_radius = 11.2 # obtain from table section 2.4.5 roskam book IV
s_t = D_o - 2*(load_radius) #[TBD]
#
w_td = 10 #[fps] touch down rate in feet per second; for FAR 25.723 certified aircraft, otherwiswe take 10
N_g = 2.0 #1.5-2 Landing gear load factor
eta_t = 0.47 # tire energy absorption efficiency
eta_s = 0.8 # show absorption efficiency; OLEO-pneumatic

def energy_absorption(g=g,w_td=w_td,eta_t=eta_t,eta_s=eta_s,s_t=s_t,N_struts=N_struts,P_mw=1.07*0.2248*P_mw/np.arctan2(z_main_lg,min(y_lg_list)),N_g=N_g,P_nw=1.07*0.2248*P_nw,P_n_dyn_t=1.07*0.2248*P_n_dyn_t):
    E_t_max = 0.5*P_mw*N_struts*w_td**2/g                   # [ft*lbs] maximum kinetic energy which needs to be absorbed #all touchdown energy absorbed by main landing gear
    s_s = ((E_t_max/(N_struts*P_mw*N_g))-eta_t*s_t)/eta_s # [ft] #Note that the gear loads are now converted from [N] to [lbs]
    s_s_des = (s_s + 1/12)*0.3048                # [meters] Design shock absorber length
    d_s = (0.041 + 0.0025*(P_mw)**0.5)*0.3048 # meters
    s_s_des_nose = (((0.5*P_nw/g*w_td**2/(P_n_dyn_t*2*N_g))-eta_t*s_t)/eta_s +1/12)*0.3048 #Stroke length nose landing gear, converted to meters
    d_s_n = (0.041 + 0.0025*(P_nw)**0.5)*0.3048 # meters
    return s_s_des, d_s, s_s_des_nose,d_s_n
s_s_des, d_s, s_s_des_nose,d_s_n = energy_absorption()
print ()
print ('The required stroke length of the shock absorption is:',s_s_des,'and the required diameter equals:',d_s,'Both in meters')
print ('Shock absorber stroke nose landing gear equals:',s_s_des_nose, 'meters, with a diameter of',d_s_n,'meters. diameter for main landing gear is taken.')
print ()
print ('TODO: estimate volume for storage, and tire selection, material selection')


# Strength of landing gear: take the shock of landing into account
# Size strut(s)
# Oleo-pneumatic shock absorbers, for their high energy absorption efficiency (tires 0.47, olea 0.8)
# Touchdown rate w_touchd = 12 feet per second (transports FAR 25.723 certified )
# Roskam source 2 page 94: OLEO-pneumatic shock absorber explanation
# same source: linearly decrease in tire presure and loadpage 51
# Just make an assumption on what tire is good enough, I have absolutely no clue what the linear relationship is...