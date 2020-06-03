# -*- coding: utf-8 -*-
''' 
file to be used for scissor plots. Being worked out atm
responsible: Jorn
'''

import numpy as np
import matplotlib.pyplot as plt
from math import sqrt, pi, tan, atan
import input
import Class_2_estimation as cl2


MAC = input.MAC
lemac = input.X_LEMAC_nose
tail_armh = input.lh

wing_area =  cl2.S#m^2
span = cl2.b
horizontal_area = input.Sh
span_tail = input.bh
 
widthf = input.widthf
fuselage_lenght = input.lf #m


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
n = 0.95
stabilizer_sweep = input.half_chord_sweep_hor

AR_tail = span_tail**2/horizontal_area
AR = input.AR
sweep = input.half_sweep

C_h = horizontal_area*tail_armh/(wing_area*MAC)


C_lh_max = -0.8 #adjustable tail



clalpha_datcom =  2*np.pi*AR/(2+((4+ ((AR*beta/n)**2)*(1+ (np.tan(sweep)**2)/beta**2))**0.5))


clalpha_tail =  2*np.pi*AR_tail/(2+((4+ ((AR_tail*beta_tail/n)**2)*(1+ (np.tan(stabilizer_sweep)**2)/beta_tail**2))**0.5))
clalpha_tail_degrees = np.radians(clalpha_tail)


clalpha_acless_wing = clalpha_datcom*(1+2.15*(widthf/span))*(wing_area-area_fuselage)/wing_area
clalpha_acless_fuselage = (np.pi*widthf**2)/(2*wing_area)# not sure bout this
clalpha_acless = clalpha_acless_wing+clalpha_acless_fuselage

r = 2*tail_armh/span
K_ea = (0.1124+0.1265*sweep+0.1766*sweep**2)/(r**2) + 0.1024/r +2
K_0 = 0.1124/(r**2)+0.1024/r +2
m_tv = 5.5*34.9*2/(87*span)

part_a = (r*0.4876)/((r**2+m_tv**2)*((r**2+m_tv**2+0.6319)**0.5))
part_b = 1+ ((r**2)/(r**2+0.7915+5.0734*m_tv**2))**0.3113
part_c = 1-((m_tv**2)/(1+m_tv**2))**0.5

downwash = (K_ea/K_0)*((part_a)+part_b*part_c)*clalpha_datcom/(np.pi*AR)





beta_low = (1-mach_app**2)**0.5
beta_low_tail = (1- 0.85*mach_app**2)**0.5

clalpha_datcom_lowspeed =  2*np.pi*AR/(2+((4+ ((AR*beta_low)**2)*(1+ (np.tan(sweep)**2)/beta_low**2))**0.5))
clalpha_acless_lowspeed = clalpha_datcom_lowspeed*(1+2.15*(widthf/span))*(wing_area-area_fuselage)/wing_area + (np.pi*widthf**2)/(2*wing_area)
clalpha_tail_lowspeed =  2*np.pi*AR_tail/(2+((4+ ((AR_tail*beta_low_tail/n)**2)*(1+ (np.tan(stabilizer_sweep)**2)/beta_low_tail**2))**0.5))

downwash_lowspeed = (K_ea/K_0)*((part_a)+part_b*part_c)*clalpha_datcom_lowspeed/(np.pi*AR)



rho = input.rho
OEW = cl2.OEWINPUT
payload = input.W_payload
MLW = cl2.M_zfw_kg
mass = MLW #kg, mlw
MTOW = cl2.MTOW_kg
CL = 2*mass*9.81/(rho*(v_approach**2)*wing_area) #approach CL
v_cruise = input.V_C_estimate #m/s
rho_cruise = input.rho_c
CL_cruise = 2*MTOW*9.81/(rho_cruise*(v_cruise**2)*wing_area) #approach CL
