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


#Inputs varient with wing position

x_start_Cr = shift.x_start_Cr
lemac = shift.x_lemac



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
beta = (1-(mach**2))**0.5
beta_tail = (1-speedratio*(mach**2))**0.5
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
ln = input.x_LEMAC_nose + 0.25*MAC - input.x_engine_start# distance between front of engine to quarter chord mac
ct = input.Ct
cr = input.Cr
taper = input.taper

sweep_LE = input.LE_sweep
C_lh_max = -0.8                                 #adjustable tail, from SEAD slides
lh = input.lh
z_position_horizontal = input.z_position_horizontal
z_position_wing = input.z_position_wing
zero_lift_angle = input.zero_lift_angle
e_tail = input.e_tail                   #Oswald efficiency factor
x_lemac_Cr = input.x_lemac_rootchord #x location of leading edge mac measured from root chord [m]



#Minor calculations with input parameters
CL = 2*mass*9.81/(rho*(v_approach**2)*S)            #approach CL
CL_cruise = 2*MTOW*9.81/(rho_cruise*(v_cruise**2)*S) #approach CL
l_fn = x_start_Cr + widthf * np.tan(sweep_LE)
# C_h = horizontal_area*tail_armh/(S*MAC)


#Datcom method to compute lift curve slopes
clalpha_datcom =  2*np.pi*AR/(2+((4+ ((AR*beta/n)**2)*(1+ (np.tan(sweep)**2)/beta**2))**0.5))
clalpha_tail =  2*np.pi*AR_tail/(2+((4+ ((AR_tail*beta_tail/n)**2)*(1+ (np.tan(stabilizer_sweep)**2)/beta_tail**2))**0.5))
clalpha_tail_degrees = np.radians(clalpha_tail)
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
beta_low_tail = (1- speedratio*mach_app**2)**0.5

clalpha_datcom_lowspeed =  2*np.pi*AR/(2+((4+ ((AR*beta_low)**2)*(1+ (np.tan(sweep)**2)/beta_low**2))**0.5))
clalpha_acless_lowspeed = clalpha_datcom_lowspeed*(1+2.15*(widthf/b))*(S-area_fuselage)/S + (np.pi*widthf**2)/(2*S)
clalpha_tail_lowspeed =  2*np.pi*AR_tail/(2+((4+ ((AR_tail*beta_low_tail/n)**2)*(1+ (np.tan(stabilizer_sweep)**2)/beta_low_tail**2))**0.5))

# Look in SEAD lecture 4 slide 33 to get xac_w from these parameters
xac_w = Xacregression(beta_A, taper, sweepbeta) # for Mach = 0.78 (cruise)
xac_w2 = Xacregression_app(beta_A_app, taper, sweepbeta_app) # for Vappr = 66.36 (approach/landing)


xac_f1 = -1.8*widthf*hf*l_fn / (S * MAC * clalpha_acless_lowspeed) #due to nose, destabilizing
xac_f2 = 0.273*widthf*(S/b)*(b-widthf)*tan(sweep) / ((1+taper)*MAC*MAC*(b+2.15*widthf)) # liftloss of intersection with wing, stabilizing
xac_n = 2*kn*bn*bn*ln / (S*MAC*clalpha_acless_lowspeed)
xac = xac_w2+xac_f1+xac_f2+xac_n

xac_f1_cruise = -1.8*widthf*hf*l_fn / (S * MAC * CLalpha_Ah) #due to nose, destabilizing
xac_f2_cruise = 0.273*widthf*(S/b)*(b-widthf)*tan(sweep) / ((1+taper)*MAC*MAC*(b+2.15*widthf)) # liftloss of intersection with wing, stabilizing
xac_n_cruise = 2*kn*bn*bn*ln / (S*MAC*CLalpha_Ah)
xac_cruise = xac_w+xac_f1+xac_f2+xac_n

tail_armh = lh + MAC * (0.25-xac_cruise)

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

assert    0.1 < downwash/(4/(AR+2)) < 1, 'downwash value not within expected range for T-Tail'


downwash_lowspeed = (K_ea/K_0)*((part_a)+part_b*part_c)*clalpha_datcom_lowspeed/(np.pi*AR)







####################### CONTROL


mu1 = 0.18
mu2 = 1.1
mu3 = 0.04
cprime_c = 1.2 
print('Read off acutal values from SEAD lecture 5 slides 18-20 once wing is designed')

DClmax = cprime_c*1.3 # Based on adsee 2
print('Change to *1.6 if double slotted flaps are used, see slide 35 ADSEE II')


def chord_along_span(Cr, Ct, b, y):
    c = Cr - (Cr - Ct) / (b / 2) * y
    return c


outboard_flap = widthf + Aero.x2

def Swf(widthf, outboard_flap):
    b_imag = outboard_flap - widthf
    swf = 2 * b_imag * (chord_along_span(cr, ct, b, widthf) + chord_along_span(cr, ct, b, outboard_flap)) / 2
    return swf

Swf = Swf(widthf, outboard_flap)

CL0_flapped = cl0+0.9*DClmax*(Swf/S)*0.975

cm_wing = cm0 *(AR *np.cos(sweep)**2)/(AR + 2*np.cos(sweep))
cm_fus = -1.8 * (1 - 2.5*widthf/fuselage_lenght)*(A_fuselage*fuselage_lenght*CL0_flapped/(4*S*MAC*clalpha_acless_lowspeed))
DCm025 = mu2*(-mu1*DClmax*cprime_c-(CL+DClmax*(1-Swf/S))*0.125*cprime_c*(cprime_c-1)) + 0.7*AR*mu3*DClmax*tan(sweep) / (1+2/AR) - CL * (0.25 - xac / MAC)
cm_flaps = DCm025 -CL*(0.25 - xac/MAC)
cm_nac = 0
cm_ac = cm_wing + cm_flaps + cm_fus + cm_nac




ShS = np.arange(0.0,0.605,0.005)
stabilityxcg_cruise = xac_cruise + ShS*(clalpha_tail/clalpha_acless)*(1-downwash)*speedratio*tail_armh/MAC
controlxcg = xac - cm_ac/CL + ShS*(C_lh_max/CL)*(tail_armh/MAC)*speedratio


def scissorplot(stabilityxcg_cruise,controlxcg, ShS, frontcg, aftcg, Sh_over_S  ):
    plt.figure()
    plt.close()
    plt.plot(stabilityxcg_cruise*100,ShS, color = 'grey', label = 'Neutral stability')
    plt.plot(stabilityxcg_cruise*100 -5,ShS, color = 'b', label = 'Stability aft limit')
    plt.plot(controlxcg*100,ShS, color = 'orange', label = 'Control fwd limit')
    plt.plot([frontcg,aftcg], [Sh_over_S, Sh_over_S], color = 'r', marker = '|')
    plt.grid()
    plt.xlabel("Xcg/MAC [%]")
    plt.ylabel("Sh/S [-]")
    plt.legend(loc = 'lower left')
    plt.title('CS100')
    plt.show()
    
scissorplot(stabilityxcg_cruise,controlxcg, ShS, frontcg, aftcg, Sh_over_S  )

Moment_ac = 0.5* rho_cruise *v_cruise**2 * cm_ac * MAC

Lift_tail = Moment_ac/tail_armh
CL_h = Lift_tail/(0.5* rho_cruise *v_cruise**2  * horizontal_area)
k = 1 / (np.pi*AR_tail *e_tail)

Dtrim = abs(0.5* rho_cruise *v_cruise**2 *speedratio * horizontal_area * CL_h * k)





#todo: check capability of horizontal tail for providing negative lift to sufficiently rotate the aircraft at take-off
