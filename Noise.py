# -*- coding: utf-8 -*-
"""
Created on Sun Jun  7 09:02:25 2020

@author: Rick van Rooijen


"""

import numpy as np
import input
import Class_2_estimation as cl2

#Fixed input parameters
Pref = 20E-6                        #Pa
c = input.MAC
U = input.v_approach
M_ap = input.mach_app
a = U / M_ap
gamma_ap = input.gamma_ap
v = 1.48E-5             #at 15 degrees celsius
rho = 1.225
S = cl2.S
CLmax_land = input.CLmax_land
CD0_w = 0.05
AR = input.AR
e = input.e
psi = input.LE_sweep
#TBD input parameters
theta = np.radians(270)
phi = np.radians(90)
mu = np.radians(180)
mu_s = np.radians(90-gamma_ap)
beta = np.radians(0)
Kw = -26
Kwf = -28
Klg = -31
R = 120 
CDw = CD0_w + CLmax_land **2 / (np.pi * AR * e)
dCDlg = (input.CD0_landGD - input.CD0) + CLmax_land **2 / (np.pi * AR * e)
#Noise functions
def OASPL_w(Kw, CDw, U, S, theta, phi, psi, Pref, a, R, c, v, rho):
    OASPL_w = Kw +  10 * np.log10((CDw * rho**2 * U**4.8 * S * (np.cos(theta/2))**2 * np.sin(phi) * (np.cos(psi))**2) / (Pref * a * R**2 * (c / v)**0.2))
    return OASPL_w

def OASPL_lg(Klg, rho, U, dCDlg, S, mu, mu_s, beta, a, R, Pref):
    D = (np.cos(mu_s))**2 * (np.cos(mu))**2 + (np.sin(mu_s))**2 * (np.sin(mu))**2 * (np.cos(beta))**2 + np.sin(mu_s) * np.cos(mu_s) * np.sin(mu) * np.cos(beta)
    OASPL_lg = Klg + 10 * np.log10((rho**2 * U**6 * dCDlg * S * D) / (a**2 * R**2 * Pref))
    return OASPL_lg

wing_noise = OASPL_w(Kw, CDw, U, S, theta, phi, psi, Pref, a, R, c, v, rho)
landing_gear_noise = OASPL_lg(Klg, rho, U, dCDlg, S, mu, mu_s, beta, a, R, Pref)
airframe_noise = wing_noise + landing_gear_noise
print(wing_noise)
print(landing_gear_noise)
print(airframe_noise)
















