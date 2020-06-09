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


#Fan noise
Ae = ...            #Engine reference area [m^2]
rs = ...            #Distance from source to observer [m]
n_rot = ...         #Number of rotor blades [-]
n_stat = ...        #Number of stator vanes [-]
N_star = 0.3        #Rotational speed ,reference value [0,0.3,0.5], re a/d
Mt = np.pi * N_star            #Fan rotor tip mach number [-]
if Mt > 1.05:
    Mt = Mt / abs(1 - n_stat / n_rot)
Mm = 1                              #Design point mach number index = max(1,Md)
A_star = ...                        #Fan inlet cross-sectional area
m_dot_star = ...                    #Mass flow rate [kg/s]
Mx = m_dot_star / A_star            #Axial flow mach number
Mr = np.sqrt(Mt**2 - Mx**2)
F = 1
if Mr <= 0.9:
    F = 0.81 * Mr**-2
    
frequency = np.arange(0,10000,1)   
eta_l = int(10**(-1/20) * eta) + 1 #Lowest harmonic number that falls within the band, calculate nodes for each frequency
eta_u = int(10**(1/20) * eta)       #Highest harmonic nubmer that falls within the band
n_har = eta_u - eta_l + 1

    
def p_squared():
    rs_star = rs / np.sqrt(Ae)
    p2 = A * Pi / (4 * np.pi * rs_star**2) * D_theta * S_eta / (1 - M_ap * np.cos(theta_e))**4 
    return P2

def inlet_broadband_noise():
    Pi = 1.552E-4 * np.exp()













