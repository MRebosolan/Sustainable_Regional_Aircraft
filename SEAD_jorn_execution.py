''' 
file to be used for scissor plots. Being worked out atm
responsible: Jorn
'''

import numpy as np
import matplotlib.pyplot as plt
from math import sqrt, pi, tan, atan
import input
import Class_2_estimation as cl2

lemac = 14.402
MAC = input.MAC
ac = lemac + 0.25*MAC
tail_armh = input.lh


wing_area =  cl2.S#m^2
span = cl2.b
horizontal_area = input.Sh
span_tail = input.bh
 
diameter = input.widthf
fuselage_lenght = input.lf #m

speedratio = input.tail_speedratio**2 #already squared

cl0 = 0.153333 
cm0 = -0.018
area_fuselage = diameter*14*(10/21)
mach = 0.78
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


clalpha_acless_wing = clalpha_datcom*(1+2.15*(diameter/span))*(wing_area-area_fuselage)/wing_area
clalpha_acless_fuselage = (np.pi*diameter**2)/(2*wing_area)# not sure bout this
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
clalpha_acless_lowspeed = clalpha_datcom_lowspeed*(1+2.15*(diameter/span))*(wing_area-area_fuselage)/wing_area + (np.pi*diameter**2)/(2*wing_area)
clalpha_tail_lowspeed =  2*np.pi*AR_tail/(2+((4+ ((AR_tail*beta_low_tail/n)**2)*(1+ (np.tan(stabilizer_sweep)**2)/beta_low_tail**2))**0.5))

downwash_lowspeed = (K_ea/K_0)*((part_a)+part_b*part_c)*clalpha_datcom_lowspeed/(np.pi*AR)


rho = 1.225
OEW = cl2.OEWINPUT
payload = input.W_payload
MLW = cl2.M_zfw_kg
mass = MLW #kg, mlw
MTOW = cl2.MTOW_kg
CL = 2*mass*9.81/(rho*(v_approach**2)*wing_area) #approach CL
v_cruise = 230 #m/s
rho_cruise = 0.348364
CL_cruise = 2*MTOW*9.81/(rho_cruise*(v_cruise**2)*wing_area) #approach CL




######################################################
# Aerodynamic center aircraft-less-tail calculations #
######################################################

b = span
S = wing_area
mac = MAC
AR = AR
ct = 1.69
cr = 6.8
taper = cr/ct
l_fn = 11.5
bf = diameter
hf = diameter
CLalpha_Ah = clalpha_acless
kn = -4         # Nacelles are mounted before wing LE
bn = 2.48
ln = 3.53

beta = sqrt(1-mach*mach)
beta_A = AR*beta
# Look in SEAD lecture 4 slide 33 to get xac_w from these parameters
xac_w = 0.345 # for Mach = 0.78 (cruise)
xac_w2 = 0.39 # for Vappr = 66.36 (approach/landing)
xac_f1 = -1.8*bf*hf*l_fn / (S * mac * clalpha_acless_lowspeed) #due to nose, destabilizing
xac_f2 = 0.273*bf*(S/b)*(b-bf)*tan(sweep) / ((1+taper)*mac*mac*(b+2.15*bf)) # liftloss of intersection with wing, stabilizing
xac_n = 2*kn*bn*bn*ln / (S*mac*clalpha_acless_lowspeed)
xac = xac_w2+xac_f1+xac_f2+xac_n

xac_f1_cruise = -1.8*bf*hf*l_fn / (S * mac * CLalpha_Ah) #due to nose, destabilizing
xac_f2_cruise = 0.273*bf*(S/b)*(b-bf)*tan(sweep) / ((1+taper)*mac*mac*(b+2.15*bf)) # liftloss of intersection with wing, stabilizing
xac_n_cruise = 2*kn*bn*bn*ln / (S*mac*CLalpha_Ah)
xac_cruise = xac_w+xac_f1+xac_f2+xac_n

print(xac_w)
print(xac_f1)
print(xac_f2)
print(xac_n)
print(xac)
print(xac_cruise)


# Based on SEAD lecture 5
mu1 = 0.2
mu2 = 0.95
mu3 = 0.044
cprime_c = 1.15 # Estimation
DClmax = cprime_c*1.3 # Based on adsee 2
# deltaf = 40*pi/180 # deflection angle in radians
Swf = 79.1 # Geometric estimation
DCm025 = mu2*(-mu1*DClmax*cprime_c-(CL+DClmax*(1-Swf/S))*0.125*cprime_c*(cprime_c-1)) + 0.7*AR*mu3*DClmax*tan(sweep) / (1+2/AR)
print(DCm025)

CL0_flapped = cl0+0.9*DClmax*(Swf/S)*0.975
cm_flaps = DCm025 -CL*(0.25 - xac/MAC)
cm_nac = 0
cm_wing = cm0 *(AR *np.cos(sweep)**2)/(AR + 2*np.cos(sweep))
cm_fus = -1.8 * (1 - 2.5*diameter/fuselage_lenght)*(np.pi*diameter*diameter*fuselage_lenght*CL0_flapped/(4*wing_area*MAC*clalpha_acless_lowspeed))
cm_fus_cruise =-1.8 * (1 - 2.5*diameter/fuselage_lenght)*(np.pi*diameter*diameter*fuselage_lenght*cl0/(4*wing_area*MAC*clalpha_acless))
cm_ac = cm_wing + cm_flaps + cm_fus + cm_nac
cm_ac_cruise = cm_wing+cm_fus_cruise



ShS = np.arange(0.0,0.605,0.005)
stabilityxcg_cruise = xac_cruise + ShS*(clalpha_tail/clalpha_acless)*(1-downwash)*speedratio*tail_armh/MAC
stabilityxcg = xac + ShS*(clalpha_tail_lowspeed/clalpha_acless_lowspeed)*(1-downwash_lowspeed)*speedratio*tail_armh/MAC
controlxcg = xac - cm_ac/CL + ShS*(C_lh_max/CL)*(tail_armh/MAC)*speedratio


vertical1 = 0.184 *100
vertical2 = 0.373 *100
# plt.plot(stabilityxcg_cruise*100,ShS, color = 'grey', label = 'Neutral stability')
# plt.plot(stabilityxcg_cruise*100 -5,ShS, color = 'b', label = 'Stability aft limit')
plt.close()
plt.subplot(1,1,1)
plt.plot(stabilityxcg_cruise*100,ShS, color = 'grey', label = 'Neutral stability')
plt.plot(stabilityxcg_cruise*100 -5,ShS, color = 'b', label = 'Stability aft limit')
plt.plot(controlxcg*100,ShS, color = 'orange', label = 'Control fwd limit')
# plt.axvline(vertical1, color = 'r', label = 'Front CG limit')
# plt.axvline(vertical2, color = 'magenta', label = 'Aft CG limit')
# plt.axhline(0.33, color = 'black',label = 'Surface ratio')
plt.plot([vertical1,vertical2], [0.33,0.33], color = 'r', marker = '|')
plt.grid()
plt.xlabel("Xcg/MAC [%]")
plt.ylabel("Sh/S [-]")
# plt.legend(loc = 'lower left')
plt.title('CS100')
plt.show()
    
