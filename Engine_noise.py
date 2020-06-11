# -*- coding: utf-8 -*-
"""
Created on Tue Jun  9 10:50:38 2020

@author:Rick van Rooijen
A_star: fan inlet cross-sectional area [m^2]
r_s: distance from source to observer [m]
A_e = engine reference area [m^2], determine from table 1?
f = frequency range for which noise is plotted
N = rotational speed of the fan[Hz]!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!not rad/s
B = number of rotor blades [-]
d = fan rotor diameter [m]
theta = polar directivity angle [rad]
V = number of stator vanes [-]
m_dot = mass flow rate [kg/s]


r_s_star = normalized distance from noise source to observer
N_star = normalized rotational speed
d_star = normalized fan rotor diameter
fb = blade passing frequency
eta = frequency parameter
M_t = rotor tip mach number
delta = cut-off factor
m_dot_star = normalized mass flow rate
M_x = axial flow mach number
M_r = relative tip mach number
eta_l = lowest harmonic number present in the 1/3 Octave Band at center frequencies f
eta_u = lowest harmonic number present in the 1/3 Octave Band at center frequencies f
n = all harmonic tones contained in the 1/3 Octave Band at center frequencies f
"""
import numpy as np

#Input from input file
c = input.MAC
U = input.v_approach
M_ap = input.mach_app
a = U / M_ap
rho = input.rho0

#Input parameters
A_star = ...
r_s = ...
A_e = ...
f = np.arange(10,10010,10)
N = ...
B = ...
d = ...
theta = ...
V = ...
m_dot = ...
#Minor computations with input parameters
r_s_star = r_s / np.sqrt(A_e)
N_star = N / (a / d)
d_star = d / np.sqrt(A_e)
fb = N_star * B * a / (d_star * np.sqrt(A_e))
eta = (1 - M_ap * np.cos(theta)) / fb * f 
M_t = np.pi * N_star
delta = M_t
if M_t > 1.05:
    delta = M_t / (abs(1 - V / B))
m_dot_star = m_dot / (rho * a * A_e) 
M_x = m_dot_star / A_star
M_r = np.sqrt(M_t**2 - M_x**2)
eta_l = [int(10**(-1/20) * i) + 1 for i in eta]#Lowest harmonic number that falls within the band, calculate nodes for each frequency
eta_u = [int(10**(1/20) * i) for i in eta]       #Highest harmonic nubmer that falls within the band
n = [np.arange(eta_l[i], eta_u[i] + 1, 1) for i in range(len(f))]

def p2_fan():
    p2 = A_star * Pi_star / (4 * np.pi * rs_star**2) * D_theta * S_eta / (1 - M_ap * np.cos(theta))**4
    

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
    


    
def p_squared():
    rs_star = rs / np.sqrt(Ae)
    p2 = A * Pi / (4 * np.pi * rs_star**2) * D_theta * S_eta / (1 - M_ap * np.cos(theta_e))**4 
    return P2

def inlet_broadband_noise():
    Pi = 1.552E-4 * np.exp()