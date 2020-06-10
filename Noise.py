# -*- coding: utf-8 -*-
"""
Created on Sun Jun  7 09:02:25 2020

@author: Rick van Rooijen
bw = wing span [m]
Aw = wing area [m^2]
Af = total flap area [m^2]
bf = total flap span [m]
delta_f = flap deflection angle [rad]
mu_inf = dynamic viscosity of ambient air, constant
f = all center frequencies of 1/3 octave band that will be plotted
U = approach velocity [m/s]
M = approach mach number [-]
c = speed of sound [m/s]
rho = density [kg/m^3]
theta = polar directivity angle [rad]
phi = azimuth directivity angle [rad]
Pref = reference effective pressure (treshold of hearing) [N/m^2]
n = number of main landing gear wheels [-]
d = main landing gear wheel diameter [m]


P = noise power function [W]
F = 

Noise of nose landing gear discarded
"""

import numpy as np
import input
import Class_2_estimation as cl2
import matplotlib.pyplot as plt


bw = cl2.b
Aw = cl2.S
Af = 20
delta_f = np.radians(30)
mu_inf = 1.84E-5
f = np.arange(10,10010,10)
U = input.v_approach
M = input.mach_app
c = U / M
rho = input.rho0
theta = np.radians(90)
phi = np.radians(0)
Pref = (2E-5)**2 
n = 4
d = input.d_wheel_main_lg


def clean_wing():
    G_clean_wing = 0.37 * Aw / bw**2 * (rho * M * c * Aw / (mu_inf * bw))**-0.2
    L_clean_wing = G_clean_wing * bw
    a_clean_wing = 5
    P_clean_wing = 4.646E-5 * M**a_clean_wing * G_clean_wing
    S_clean_wing = [i * L_clean_wing * (1 - M * np.cos(theta)) / (M * c) for i in f]
    FS_clean_wing = [0.613 * (10 * i)**4 * ((10 * i)**1.5 + 0.5)**-4 for i in S_clean_wing]
    D_clean_wing = 4 * (np.cos(phi))**2 * (np.cos(theta/2))**2
    pe2_clean_wing = [rho * c * P_clean_wing * D_clean_wing * i / (4 * np.pi * 1**2 * (1 - M * np.cos(theta))**4) for i in FS_clean_wing]                  #per 1 m distance, compensate for 1/r^2
    SPL_clean_wing = [10 * np.log10(i / Pref) for i in pe2_clean_wing]
    return SPL_clean_wing

def flaps():
    G_flaps = Af / bw**2 * (np.sin(delta_f))**2
    L_flaps = G_flaps * bw
    a_flaps = 6
    P_flaps = 2.787E-4 * M**a_flaps * G_flaps
    S_flaps = [i * L_flaps * (1 - M * np.cos(theta)) / (M * c) for i in f]
    FS_flaps = []
    for j in range(len(S_flaps)):
        if S_flaps[j] < 2:
            FS_flaps.append(0.0480 * S_flaps[j])
        elif S_flaps[j] > 20:
            FS_flaps.append(216.49 * S_flaps[j]**-3)
        else:
            FS_flaps.append(0.1406 * S_flaps[j]**-0.55)
    D_flaps = 3 * (np.sin(delta_f) * np.cos(theta) + np.cos(delta_f) * np.sin(theta) * np.cos(phi))**2
    pe2_flaps = [rho * c * P_flaps * D_flaps * i / (4 * np.pi * 1**2 * (1 - M * np.cos(theta))**4) for i in FS_flaps]
    SPL_flaps = [10 * np.log10(i / Pref) for i in pe2_flaps]    
    return SPL_flaps
    
def landing_gear():
    G_lg = n * (d / bw)**2
    L_lg = d
    a_lg = 6
    P_lg  =  3.414E-4 * M**a_lg * G_lg 
    S_lg = [i * L_lg * (1 - M * np.cos(theta)) / (M * c) for i in f]
    FS_lg = [0.0577 * i**2 * (0.25 * i**2 + 1)**-1.5 for i in S_lg]
    D_lg = 3 / 2 * (np.sin(theta))**2
    pe2_lg = [rho * c * P_lg * D_lg * i / (4 * np.pi * 1**2 * (1 - M * np.cos(theta))**4) for i in FS_lg]
    SPL_lg = [10 * np.log10(i / Pref) for i in pe2_lg]
    return SPL_lg
    


SPL_clean_wing = clean_wing() 
SPL_flaps = flaps() 
SPL_lg = landing_gear() 

plt.figure()
plt.plot(f, SPL_clean_wing, f, SPL_flaps, f, SPL_lg)
plt.xscale('log')
plt.show()

#
##Fixed input parameters
#Pref = 20E-6                        #Pa
#c = input.MAC
#U = input.v_approach
#M_ap = input.mach_app
#a = U / M_ap
#gamma_ap = input.gamma_ap
#v = 1.48E-5             #at 15 degrees celsius
#rho = 1.225
#S = cl2.S
#CLmax_land = input.CLmax_land
#CD0_w = 0.05
#AR = input.AR
#e = input.e
#psi = input.LE_sweep
##TBD input parameters
#theta = np.radians(270)
#phi = np.radians(90)
#mu = np.radians(180)
#mu_s = np.radians(90-gamma_ap)
#beta = np.radians(0)
#Kw = -26
#Kwf = -28
#Klg = -31
#R = 120 
#CDw = CD0_w + CLmax_land **2 / (np.pi * AR * e)
#dCDlg = (input.CD0_landGD - input.CD0) + CLmax_land **2 / (np.pi * AR * e)
##Noise functions
#def OASPL_w(Kw, CDw, U, S, theta, phi, psi, Pref, a, R, c, v, rho):
#    OASPL_w = Kw +  10 * np.log10((CDw * rho**2 * U**4.8 * S * (np.cos(theta/2))**2 * np.sin(phi) * (np.cos(psi))**2) / (Pref * a * R**2 * (c / v)**0.2))
#    return OASPL_w
#
#def OASPL_lg(Klg, rho, U, dCDlg, S, mu, mu_s, beta, a, R, Pref):
#    D = (np.cos(mu_s))**2 * (np.cos(mu))**2 + (np.sin(mu_s))**2 * (np.sin(mu))**2 * (np.cos(beta))**2 + np.sin(mu_s) * np.cos(mu_s) * np.sin(mu) * np.cos(beta)
#    OASPL_lg = Klg + 10 * np.log10((rho**2 * U**6 * dCDlg * S * D) / (a**2 * R**2 * Pref))
#    return OASPL_lg
#
#wing_noise = OASPL_w(Kw, CDw, U, S, theta, phi, psi, Pref, a, R, c, v, rho)
#landing_gear_noise = OASPL_lg(Klg, rho, U, dCDlg, S, mu, mu_s, beta, a, R, Pref)
#airframe_noise = wing_noise + landing_gear_noise
#print(wing_noise)
#print(landing_gear_noise)
#print(airframe_noise)
















