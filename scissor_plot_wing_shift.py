# -*- coding: utf-8 -*-
"""
Created on Fri Jun  5 09:04:23 2020

@author: Gebruiker
"""

# -*- coding: utf-8 -*-
''' 
file to be used for scissor plots. Being worked out atm
responsible: Jorn & Rick
'''

import numpy as np
import matplotlib.pyplot as plt
import input
from math import sqrt, pi, tan, atan
import Class_2_estimation as cl2
from Xacregression_scissor import Xacregression, Xacregression_app
import Aero
import CG_excursion_wing_shift as shift

#Inputs invarient with wing position
MAC = input.MAC
S =  cl2.S
b = cl2.b
widthf = input.widthf
fuselage_lenght = input.lf #m
hf = input.hf
A_fuselage = input.A_fuselage
speedratio = input.tail_speedratio**2 #already squared
cl0 = input.cl0
cm0 = input.cm0
Chord_root = input.Cr
area_fuselage = widthf * Chord_root
mach = input.mach_cruise
v_approach = input.v_approach
mach_app = input.mach_app
n = 0.95                                            #airfoil efficiency
stabilizer_sweep = input.half_chord_sweep_hor
AR_tail = input.AR_h
AR = input.AR
sweep = input.quarter_sweep
rho = input.rho
OEW = cl2.OEWINPUT
payload = input.W_payload
MLW = cl2.M_zfw_kg                      ### STILL SIZE MLW
mass = MLW                              #kg, mlw
MTOW = cl2.MTOW_kg
v_cruise = input.V_C_estimate                #m/s
rho_cruise = input.rho_c
kn = -4                                         # Nacelles are mounted before wing LE, change in case of fuselage mounted engines
bn = input.bn
ct = input.Ct
cr = input.Cr
taper = input.taper
sweep_LE = input.LE_sweep
C_lh_max = input.cl_htail_max                                 #adjustable tail, from SEAD slides
lh_fix = input.lh                   #distance between wing and horizontal tail aerodynamic centers
ln = 0.25 * MAC - input.x_engine_start        # distance between front of engine to quarter chord mac
z_position_horizontal = input.z_position_horizontal
z_position_wing = input.z_position_wing
zero_lift_angle = input.zero_lift_angle
e_tail = input.e_tail                   #Oswald efficiency factor
x_lemac_Cr = input.x_lemac_rootchord         #x location of leading edge mac measured from root chord [m]

#Inputs varient with wing position

x_start_Cr = shift.x_start_Cr
lemac = shift.x_lemac
xh = input.x_lemac_rootchord_h + input.x_rootchord_h + 0.25*input.c_mac_h
lh = [xh - i - 0.25*MAC for i in x_start_Cr]


#Functions
def chord_along_span(Cr, Ct, b, y):
    c = Cr - (Cr - Ct) / (b / 2) * y
    return c

def swf(widthf, outboard_flap):
    b_imag = outboard_flap - widthf 
    swf = 2 * b_imag * (chord_along_span(cr, ct, b, widthf) + chord_along_span(cr, ct, b, outboard_flap)) / 2
    return swf

def trimdrag(cm_ac, tail_armh, horizontal_area, CL, xcg, xac, MAC):
    # Moment_ac = 0.5* rho_cruise *v_cruise**2 * cm_ac * MAC *S
    
    # Lift_tail = Moment_ac/tail_armh
    # CL_h = Lift_tail/(0.5* rho_cruise *v_cruise**2  * horizontal_area)
    
    
    CL_h = ( cm_ac + ((xcg - xac)*CL))* S * MAC/ (horizontal_area*tail_armh)
    
    k = 1 / (np.pi*AR_tail *e_tail)
    
    
    
    Dtrim = 0.5* rho_cruise *v_cruise**2 *speedratio * horizontal_area * CL_h**2 * k
    
    return Dtrim, CL_h
#Import cg ranges loading diagram due to wing shifting
cg_fwd_lst = shift.cg_fwd_excursion_lst
cg_aft_lst = shift.cg_aft_excursion_lst

#cg_fwd_lst = [0.3]
#cg_aft_lst = [1.9]


def main_lg_loc(x_cg_aft, theta=15, z_cg=input.z_cg,safetymargin_theta=0.01, x_tailcone=input.lf - shift.tailcone_length):
    d = 0.01
    for z_fus_ground in np.arange(0,5,d):
        
        z_main_lg1 = z_fus_ground + z_cg
        x_main_lg1 = x_cg_aft + z_main_lg1*np.tan(np.radians(theta+safetymargin_theta))
        if np.arctan(z_fus_ground/(x_tailcone-x_main_lg1)) >= np.radians(15):
            z_main_lg = np.round(z_main_lg1,4)
            x_main_lg = np.round(x_main_lg1,4)
            z_f_ground = z_fus_ground
            break
        else:
            x_main_lg = 999
            continue

    return x_main_lg     

def req_htail_area(x_cg_aft,x_ac_htail, x_cg_front, front_weight,x_ac, momentcoefficient, Cl_htail = input.cl_htail_max,rho_to=rho, Vlof=59.1, MAC = MAC): 
    x_main_LG = main_lg_loc(x_cg_aft)
    
    CL_ground = input.CL0takeoff 
    cg_contribution = (x_main_LG-x_cg_front)*front_weight*9.81
    lift_con = (- 0.5 *rho_to * Vlof * Vlof * S * CL_ground)*(x_main_LG - x_ac)
    moment_con = - momentcoefficient*.5*rho_to*(Vlof**2)*S*MAC
    
    htail_moment =  -(0.5*rho_to*(Vlof**2)*Cl_htail*(x_ac_htail-x_main_LG))
    
    htail_area = (cg_contribution + lift_con + moment_con)/htail_moment
    sh_over_s = htail_area/S
    return sh_over_s, x_main_LG



print('Read off acutal values for cprime_c from SEAD lecture 5 slides 18-20 once wing is designed')
print('Change to *1.6 in DCLmax if double slotted flaps are used, see slide 35 ADSEE II')
def scissor_wing_shift():
    Sh_min_lst = []
    rotation = []
    for i in range(len(x_start_Cr)):
    #for i in range(len(cg_fwd_lst)):
        
        #Minor calculations with input parameters
        CL = input.Clmax_Lnd            #approach CL
        #CL_cruise = 2*MTOW*9.81/(rho_cruise*(v_cruise**2)*S) #approach CL
        l_fn = x_start_Cr[i] + widthf * np.tan(sweep_LE)
        beta = (1-(mach**2))**0.5
        beta_tail = (1-speedratio*(mach**2))**0.5
        # C_h = horizontal_area*tail_armh/(S*MAC)
    
    
        #Datcom method to compute lift curve slopes
        clalpha_datcom =  2*np.pi*AR/(2+((4+ ((AR*beta/n)**2)*(1+ (np.tan(sweep)**2)/beta**2))**0.5)) 
        clalpha_tail =  2*np.pi*AR_tail/(2+((4+ ((AR_tail*beta_tail/n)**2)*(1+ (np.tan(stabilizer_sweep)**2)/beta_tail**2))**0.5))
        #clalpha_tail_degrees = np.radians(clalpha_tail)
        clalpha_acless_wing = clalpha_datcom*(1+2.15*(widthf/b))*(S-area_fuselage)/S
        clalpha_acless_fuselage = (np.pi*widthf**2)/(2*S)# not sure bout this
        clalpha_acless = clalpha_acless_wing+clalpha_acless_fuselage
        CLalpha_Ah = clalpha_acless



        beta = sqrt(1-mach*mach)
        beta_A = AR*beta
        sweepbeta = np.degrees(sweep)/beta
        beta_app = sqrt(1-mach_app*mach_app)
        beta_A_app = AR*beta_app
        sweepbeta_app = np.degrees(sweep)/beta_app
        beta_low = (1-mach_app**2)**0.5
        #beta_low_tail = (1- speedratio*mach_app**2)**0.5
        
        clalpha_datcom_lowspeed =  2*np.pi*AR/(2+((4+ ((AR*beta_low)**2)*(1+ (np.tan(sweep)**2)/beta_low**2))**0.5))
        clalpha_acless_lowspeed = clalpha_datcom_lowspeed*(1+2.15*(widthf/b))*(S-area_fuselage)/S + (np.pi*widthf**2)/(2*S)
        #clalpha_tail_lowspeed =  2*np.pi*AR_tail/(2+((4+ ((AR_tail*beta_low_tail/n)**2)*(1+ (np.tan(stabilizer_sweep)**2)/beta_low_tail**2))**0.5))

        # Look in SEAD lecture 4 slide 33 to get xac_w from these parameters
        xac_w = Xacregression(beta_A, taper, sweepbeta) # for Mach = 0.78 (cruise)
        xac_w2 = Xacregression_app(beta_A_app, taper, sweepbeta_app) # for Vappr = 66.36 (approach/landing)
        
        
        xac_f1 = -1.8*widthf*hf*l_fn / (S * MAC * clalpha_acless_lowspeed) #due to nose, destabilizing
        xac_f2 = 0.273*widthf*(S/b)*(b-widthf)*tan(sweep) / ((1+taper)*MAC*MAC*(b+2.15*widthf)) # liftloss of intersection with wing, stabilizing
        xac_n = 2*kn*bn*bn*ln / (S*MAC*clalpha_acless_lowspeed)
        xac = xac_w2+xac_f1+xac_f2+xac_n
        
        xac_f1_cruise = -1.8*widthf*hf*l_fn / (S * MAC * CLalpha_Ah) #due to nose, destabilizing
        xac_f2_cruise = xac_f2
        xac_n_cruise = 2*kn*bn*bn*ln / (S*MAC*CLalpha_Ah)
        xac_cruise = xac_w+xac_f1_cruise+xac_f2_cruise+xac_n_cruise

        tail_armh = lh[i] + MAC * (0.25-xac_cruise)
        
        r = 2*tail_armh/b
        K_ea = (0.1124+0.1265*sweep+0.1766*sweep**2)/(r**2) + 0.1024/r +2 #ADSEE LECTURE 4 SLIDE 43
        K_0 = 0.1124/(r**2)+0.1024/r +2 #ADSEE LECTURE 4 SLIDE 43
        
        theta = np.tanh((z_position_horizontal - z_position_wing)/tail_armh)
        hypotenuse = tail_armh/np.cos(theta)
        tail_wing_distance = hypotenuse*np.cos(theta+zero_lift_angle)
        m_tv = tail_wing_distance *2/b
        
        part_a = (r*0.4876)/((r**2+m_tv**2)*((r**2+m_tv**2+0.6319)**0.5)) #ADSEE LECTURE 4 SLIDE 43
        part_b = 1+ ((r**2)/(r**2+0.7915+5.0734*m_tv**2))**0.3113 
        part_c = 1-((m_tv**2)/(1+m_tv**2))**0.5
        
        downwash = (K_ea/K_0)*((part_a)+part_b*part_c)*clalpha_datcom/(np.pi*AR)
        

        # assert    0.01 < downwash/(4/(AR+2)) < 2, 'downwash value not within expected range for T-Tail'
        
        
        #downwash_lowspeed = (K_ea/K_0)*((part_a)+part_b*part_c)*clalpha_datcom_lowspeed/(np.pi*AR)

####################### CONTROL


        mu1 = 0.18
        mu2 = 0.5
        mu3 = 0.05
        cprime_c = 1.2 
        
        
        DClmax = cprime_c*1.3 # Based on adsee 2
        

        outboard_flap = widthf + Aero.x2



        Swf = swf(widthf, outboard_flap)

        CL0_flapped = cl0+0.9*DClmax*(Swf/S)*0.975

        cm_wing = cm0 *(AR *np.cos(sweep)**2)/(AR + 2*np.cos(sweep))
        cm_fus = -1.8 * (1 - 2.5*widthf/fuselage_lenght)*(A_fuselage*fuselage_lenght*CL0_flapped/(4*S*MAC*clalpha_acless_lowspeed))
        cm_fus_cruise = -1.8 * (1 - 2.5*widthf/fuselage_lenght)*(A_fuselage*fuselage_lenght*CL0_flapped/(4*S*MAC*clalpha_acless))
        DCm025 = mu2*(-mu1*DClmax*cprime_c-(CL+DClmax*(1-Swf/S))*0.125*cprime_c*(cprime_c-1)) + 0.7*AR*mu3*DClmax*tan(sweep) / (1+2/AR) - CL * (0.25 - xac / MAC)
        
        DCm025_takeoff = mu2*(-mu1*DClmax*cprime_c-(input.CL0takeoff+DClmax*(1-Swf/S))*0.125*cprime_c*(cprime_c-1)) + 0.7*AR*mu3*DClmax*tan(sweep) / (1+2/AR) - input.CL0takeoff * (0.25 - xac / MAC)

        cm_flaps = DCm025 -CL*(0.25 - xac/MAC)
        cm_flaps_takeoff = DCm025_takeoff -input.CL0takeoff*(0.25 - xac/MAC)
        cm_nac = 0
        cm_ac = cm_wing + cm_flaps + cm_fus + cm_nac
        cm_cruise = cm_wing+cm_fus_cruise
        cm_ac_takeoff = cm_wing + cm_flaps_takeoff + cm_fus + cm_nac



        #ShS = np.arange(0.0,0.605,0.005)
        ShS = np.arange(0.0,50,0.005)
   
        stabilityxcg_cruise = xac_cruise + ShS*(clalpha_tail/clalpha_acless)*(1-downwash)*speedratio*tail_armh/MAC -0.05
        controlxcg = xac - cm_ac/CL + ShS*(C_lh_max/CL)*(tail_armh/MAC)*speedratio +0.05

        cg_stab = stabilityxcg_cruise[::-1]         #Reverse array for looping over it
        cg_cont = controlxcg[::-1]
        ShS = ShS[::-1]
        for j in range(len(ShS)):

              if cg_cont[-1] - cg_stab[-1] < 0:
                Sh_min_lst.append([10,0,0,0,0,0, 0, 0])         #append a zero if this condition is not met
                print('At a root chord position of', x_start_Cr[i],' [m], the scissor plot shows no intersection')
                break 
              if cg_cont[j] - cg_fwd_lst[i] > 0 or cg_stab[j] - cg_aft_lst[i] < 0:   #in this case, the cg range does not meet the stability or contorllability requirements
                if j==0: 
                    Sh_min = ShS[j]*S
                    Sh_min_lst.append([ShS[j],x_start_Cr[i], cg_stab[j], cg_aft_lst[i], cg_cont[j], cg_fwd_lst[i], 0, cg_cont, cg_stab,i, cm_ac_takeoff, xac, cm_cruise , xac_cruise, CL])
                    print(cg_fwd_lst[i])
                    break
                else:
                    Sh_min = ShS[j-1]*S
                    Sh_min_lst.append([ShS[j-1],x_start_Cr[i], cg_stab[j-1], cg_aft_lst[i], cg_cont[j-1], cg_fwd_lst[i], 0, cg_cont, cg_stab,i, cm_ac_takeoff, xac, cm_cruise , xac_cruise, CL])
                    break
              else:
                continue
        xlemac = input.x_lemac_rootchord + x_start_Cr[i]
        rotationshs = req_htail_area(cg_aft_lst[i]*MAC + xlemac , xh, cg_fwd_lst[i]*MAC + xlemac, MTOW-1750, shift.x_ac[i], cm_ac_takeoff)
        rotation.append([rotationshs[0], rotationshs[1],cg_aft_lst[i]*MAC + xlemac,shift.x_ac[i] + shift.lh_fix, cg_fwd_lst[i]*MAC + xlemac,  MTOW, shift.x_ac[i], cm_ac_takeoff])
    
    
    
    combi = []
    for i in range(len(rotation)):
        combi.append(max([np.array(rotation)[i,0], np.array(Sh_min_lst)[i,0]]))
            
    lowest = combi.index(min(combi))#+Sh_min_lst.index(min(Sh_min_lst)))/2)
    x = 19
    # minimum = Sh_min_lst[51 ]
    # min_Sh_over_S = combi[51] * 1# (1+input.horizontal_margin)
    
    minimum = min(Sh_min_lst)
    min_Sh_over_S = minimum[0]
    
    Sh_min = min_Sh_over_S * S
    x_Cr_opt_nose = minimum[1]
    cg_stab_lim = minimum[2] 
    cg_aft = minimum[3] 
    cg_cont_lim = minimum[4] 
    cg_fwd = minimum[5] 
    print(cg_fwd)
    controlplot = minimum[7] 
    stabilityplot = minimum[8]
    index = minimum[9]
    momentcoefficient = minimum[10]
    aerodynamiccenter = minimum[11]
    cm_cruise = minimum[12]
    xac_cruise = minimum[13]
    CLalpha_Ah = minimum[14]
    
    tailarm = lh[index] + MAC * (0.25-xac_cruise)
    
    Dtrim = trimdrag(cm_cruise, tail_armh, Sh_min,  CL, cg_aft, xac_cruise, MAC), trimdrag(cm_cruise, tail_armh, Sh_min,  CL, cg_fwd, xac_cruise, MAC)

    

    rotation = np.array(rotation)

    
    return Sh_min_lst, min_Sh_over_S, x_Cr_opt_nose, cg_stab_lim, cg_aft, cg_cont_lim, cg_fwd, Dtrim, Sh_min, controlplot, stabilityplot, ShS,index, momentcoefficient, aerodynamiccenter, rotation, combi




def scissorplot(stabilityplot,controlplot, ShS, frontcg, aftcg, Sh_over_S):
    plt.close()
    plt.figure()
    plt.plot(stabilityplot*100 +5,ShS, color = 'grey', label = 'Neutral stability')
    plt.plot(stabilityplot*100,ShS, color = 'b', label = 'Stability aft limit')
    plt.plot(controlplot*100,ShS, color = 'orange', label = 'Control margin')
    plt.plot(controlplot*100 -5,ShS, color = 'grey', label = 'Control fwd limit')

    plt.plot([frontcg*100,aftcg*100], [Sh_over_S, Sh_over_S], color = 'r', marker = '|', label = 'CG range')
    plt.ylim(-0.05,0.5)
    plt.xlim(0, 100)
    plt.grid()
    plt.xlabel("Xcg/MAC [%]")
    plt.ylabel("Sh/S [-]")
    plt.legend(loc = 'upper right')
    # plt.title('Tail Tank and Drop Tank')
    plt.show()
    



Sh_min_lst, min_Sh_over_S, x_Cr_opt_nose, cg_stab_lim, cg_aft, cg_cont_lim, cg_fwd, Dtrim, Sh_min, controlplot, stabilityplot, ShS, index, momentcoefficient, aerodynamiccenter, rotation, combi = scissor_wing_shift()
scissorplot(stabilityplot, controlplot, ShS, cg_fwd, cg_aft, min_Sh_over_S)
print(min_Sh_over_S, 'sh over s')
print(x_Cr_opt_nose, 'root location')

#todo: check capability of horizontal tail for providing negative lift to sufficiently rotate the aircraft at take-off

cg_fully_loaded = shift.cg_loaded_lst[index]



xlemac = x_Cr_opt_nose + input.x_lemac_rootchord
cg_loaded_nose =  cg_fully_loaded #this was first converted, but should actually be like this 
x_ac_h_nose = shift.x_ac[index] + shift.lh_fix
x_ac_v_nose = shift.x_ac[index] + shift.lv_fix
x_ac_wing_approx = shift.x_ac[index]
print(cg_loaded_nose)
print("Check that the final cg position after fuel loading is the same as above value for cg_loaded_nose, check with loading diagram after manually changing all input parameters")