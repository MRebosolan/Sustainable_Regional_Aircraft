## -*- coding: utf-8 -*-
#"""
#Created on Tue Jun  9 18:37:28 2020
#
#@author: Rick van Rooijen
#
#



# -*- coding: utf-8 -*-
"""
Created on Tue Jun  9 18:37:28 2020

@author: Rick van Rooijen

rs = distance from noise source to observer
delta_w = turbulent boundary layer thickness
PI_TE = wing trailing edge noise accoustic power
PI_h = horizontal tail noise accoustic power
PI_v = vertical tail noise accoustic power
PI_f = flap trialing edge noise accoustic power
PI_nlg = nose landing gear noise accoustic power
PI_s_nlg = landing gear strut noise accoustic power
PI_mlg = main landing gear noise accoustic power
"""

import numpy as np
import matplotlib.pyplot as plt
import Class_2_estimation as cl2
import input 

rs = 104.815

bw = cl2.b
Aw = cl2.S
Ah = 4
Av = 1.7
Af = 2 * 5.97
flap_defl = np.radians(60)
bh = 9
bv = 5
bf = 6
mu_inf = 1.84E-5
frequency = np.arange(10,100010,10)
U = input.v_approach
M = input.mach_app
c = U / M
rho = input.rho0
theta = np.radians(90)
phi = np.radians(10)
n_nlg = 2
d_nlg = input.d_wheel_nose_lg
l_nlg = input.strut_length_nose_lg
n_mlg = 4
d_mlg = input.d_wheel_main_lg
l_mlg = input.strut_length_main_lg

Pref = 2E-5

rs_star = rs / bw
#===================================================================================
#Airframe noise
#===================================================================================
delta_w = 0.37 * Aw / bw**2 * (rho * M * c * Aw / (mu_inf * bw))**-0.2
delta_h = 0.37 * Ah / bh**2 * (rho * M * c * Ah / (mu_inf * bh))**-0.2
delta_v = 0.37 * Av / bv**2 * (rho * M * c * Av / (mu_inf * bw))**-0.2

PI_TE = 4.464E-5 * M**5 * delta_w       
PI_h = 4.464E-5 * M**5 * delta_h * (bh / bw)**2
PI_v = 4.464E-5 * M**5 * delta_v * (bv / bw)**2
PI_f = 2.787E-4 * M**6 * Af / bw**2 * (np.sin(flap_defl))**2
PI_nlg = 4.349E-4 * M**6 * n_nlg * (d_nlg / bw)**2
PI_s_nlg = 2.573E-4 * M**6 * (d_nlg / bw)**2 * l_nlg / d_nlg
PI_s_mlg = 2.573E-4 * M**6 * (d_mlg / bw)**2 * l_mlg / d_mlg     #Add this contribution twice later with log dB function!!!!!!!!!!!!!!!!!!!!!!!!!!
PI_mlg = 3.414E-4 * M**6 * n_mlg * (d_mlg / bw)**2

D_w = 4 * (np.cos(phi))**2 * (np.cos(theta/2))**2
D_h = D_w
D_v = 4 * (np.sin(phi))**2 * (np.cos(theta/2))**2
D_f = 3 * (np.sin(flap_defl) * np.cos(theta) + np.cos(flap_defl) * np.sin(theta) * np.cos(phi))**2
D_nlg = 3 / 2 * (np.sin(theta))**2
D_mlg = D_nlg
D_s_nlg = 3 * (np.sin(theta))**2 * (np.sin(phi))**2
D_s_mlg = D_s_nlg

S_w = [f * delta_w * bw / (M * c) * (1 - M * np.cos(theta)) for f in frequency]
S_h = [f * delta_h * bh / (M * c) * (1 - M * np.cos(theta)) for f in frequency]
S_v = [f * delta_v * bv / (M * c) * (1 - M * np.cos(theta)) for f in frequency]
S_f = [f * Af / (M * bf * c) * (1 - M * np.cos(theta)) for f in frequency]
S_nlg = [f * d_nlg / (M * bf * c) * (1 - M * np.cos(theta)) for f in frequency]
S_mlg = [f * d_mlg / (M * bf * c) * (1 - M * np.cos(theta)) for f in frequency]
S_s_nlg = [f * l_nlg / (M * bf * c) * (1 - M * np.cos(theta)) for f in frequency]           #Possibly change l_nlg to d_nlg
S_s_mlg = [f * l_mlg / (M * bf * c) * (1 - M * np.cos(theta)) for f in frequency]           #Possibly change l_mlg to d_mlg

F_w = [0.613 * (10 * S)**4 * ((10 * S)**1.5 + 0.5)**-4 for S in S_w]                                        
F_h = [0.613 * (10 * S)**4 * ((10 * S)**1.5 + 0.5)**-4 for S in S_h]
F_v = [0.613 * (10 * S)**4 * ((10 * S)**1.5 + 0.5)**-4 for S in S_v]
F_f = []
for i in range(len(frequency)):
    if S_f[i] < 2:
        F_f.append(0.0480 * S_f[i])
    elif S_f[i] > 20:
        F_f.append(216.49 * S_f[i]**-3)
    else:
        F_f.append(0.1406 * S_f[i]**-0.55)
F_nlg = [13.59 * S**2 * (12.5 + S**2)**-2.25 for S in S_nlg]
F_s_nlg = [5.325 * S**2 * (30 + S**8)**-1 for S in S_s_nlg]
F_mlg = [0.0577 * S**2 * (1 + 0.25 * S**2)**-1.5 for S in S_mlg]
F_s_mlg = [1.280 * S**3 * (1.06 + S**2)**-3 for S in S_s_mlg]


#Possibly make this an array by making rs an array
P2_w_lst = []
P2_h_lst = []
P2_v_lst = []
P2_f_lst = []
P2_nlg_lst = []
P2_mlg_lst = []
P2_s_nlg_lst = []
P2_s_mlg_lst = []
P2_nlg_tot_lst = []
P2_mlg_tot_lst = []
P2_aiframe_lst = []

for i in range(len(frequency)):
    P2_w_lst.append(PI_TE / (4 * np.pi * rs_star**2) * D_w * F_w[i] / (1 - M * np.cos(theta))**4)
    P2_h_lst.append(PI_h / (4 * np.pi * rs_star**2) * D_h * F_h[i] / (1 - M * np.cos(theta))**4)
    P2_v_lst.append(PI_v / (4 * np.pi * rs_star**2) * D_v * F_v[i] / (1 - M * np.cos(theta))**4)
    P2_f_lst.append(PI_f / (4 * np.pi * rs_star**2) * D_f * F_f[i] / (1 - M * np.cos(theta))**4)
    P2_nlg_lst.append(PI_nlg / (4 * np.pi * rs_star**2) * D_nlg * F_nlg[i] / (1 - M * np.cos(theta))**4)
    P2_mlg_lst.append(PI_mlg / (4 * np.pi * rs_star**2) * D_mlg * F_mlg[i] / (1 - M * np.cos(theta))**4)
    P2_s_nlg_lst.append(PI_s_nlg / (4 * np.pi * rs_star**2) * D_s_nlg * F_s_nlg[i] / (1 - M * np.cos(theta))**4)
    P2_s_mlg_lst.append(PI_s_mlg / (4 * np.pi * rs_star**2) * D_s_mlg * F_s_mlg[i] / (1 - M * np.cos(theta))**4)

P2_nlg_tot_lst = [P2_nlg_lst[i] + P2_s_nlg_lst[i] for i in range(len(frequency))]
P2_mlg_tot_lst = [P2_mlg_lst[i] + P2_s_mlg_lst[i] for i in range(len(frequency))]
P2_airframe_lst = [P2_w_lst[i] + P2_h_lst[i] + P2_v_lst[i] + P2_f_lst[i] + P2_nlg_tot_lst[i] + P2_mlg_tot_lst[i] for i in range(len(frequency))]



SPL_w_lst = [10 * np.log10(P2) + 20 * np.log10(rho * c**2 / Pref) for P2 in P2_w_lst] 
SPL_h_lst = [10 * np.log10(P2) + 20 * np.log10(rho * c**2 / Pref) for P2 in P2_h_lst]
SPL_v_lst = [10 * np.log10(P2) + 20 * np.log10(rho * c**2 / Pref) for P2 in P2_v_lst]
SPL_f_lst = [10 * np.log10(P2) + 20 * np.log10(rho * c**2 / Pref) for P2 in P2_f_lst]
SPL_nlg_lst = [10 * np.log10(P2) + 20 * np.log10(rho * c**2 / Pref) for P2 in P2_nlg_lst]
SPL_mlg_lst = [10 * np.log10(P2) + 20 * np.log10(rho * c**2 / Pref) for P2 in P2_mlg_lst]
SPL_s_nlg_lst = [10 * np.log10(P2) + 20 * np.log10(rho * c**2 / Pref) for P2 in P2_s_nlg_lst]
SPL_s_mlg_lst = [10 * np.log10(P2) + 20 * np.log10(rho * c**2 / Pref) for P2 in P2_s_mlg_lst]
SPL_nlg_tot_lst = [10 * np.log10(P2) + 20 * np.log10(rho * c**2 / Pref) for P2 in P2_nlg_tot_lst] 
SPL_mlg_tot_lst = [10 * np.log10(P2) + 20 * np.log10(rho * c**2 / Pref) for P2 in P2_mlg_tot_lst]
SPL_airframe_lst = [10 * np.log10(P2) + 20 * np.log10(rho * c**2 / Pref) for P2 in P2_airframe_lst]




SPL_w_lst = np.ma.masked_where(np.array(SPL_w_lst) < 0, np.array(SPL_w_lst))
SPL_h_lst = np.ma.masked_where(np.array(SPL_h_lst) < 0, np.array(SPL_h_lst))
SPL_v_lst = np.ma.masked_where(np.array(SPL_v_lst) < 0, np.array(SPL_v_lst))
SPL_f_nlg_lst = np.ma.masked_where(np.array(SPL_f_lst) < 0, np.array(SPL_f_lst))
SPL_nlg_lst = np.ma.masked_where(np.array(SPL_nlg_lst) < 0, np.array(SPL_nlg_lst))
SPL_mlg_lst = np.ma.masked_where(np.array(SPL_mlg_lst) < 0, np.array(SPL_mlg_lst))
SPL_s_nlg_lst = np.ma.masked_where(np.array(SPL_s_nlg_lst) < 0, np.array(SPL_s_nlg_lst))
SPL_s_mlg_lst = np.ma.masked_where(np.array(SPL_s_mlg_lst) < 0, np.array(SPL_s_mlg_lst))
SPL_nlg_tot_lst = np.ma.masked_where(np.array(SPL_nlg_tot_lst) < 0, np.array(SPL_nlg_tot_lst))
SPL_mlg_tot_lst = np.ma.masked_where(np.array(SPL_mlg_tot_lst) < 0, np.array(SPL_mlg_tot_lst))
SPL_airframe_lst = np.ma.masked_where(np.array(SPL_airframe_lst) < 0, np.array(SPL_airframe_lst))





#===================================================================================
#Engine fan noise
#===================================================================================
"""
A: fan inlet cross-sectional area [m^2]
r_s: distance from source to observer [m]
A_e = engine reference area [m^2], determine from table 1?
f = frequency range for which noise is plotted
N = rotational speed of the fan[Hz]!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!not rad/s
B = number of rotor blades [-]
d = fan rotor diameter [m]
theta = polar directivity angle [rad]
V = number of stator vanes [-]
m_dot = mass flow rate [kg/s]
s = rotor stator spacing [m]
C = mean rotor blade chord [m]
delta_T = total temperature rise across the fan [K]

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


#Input parameters
A = 4.100619584
A_e = np.pi / 4 #from table default value
N = 50
B = 20
d = 2.285
V = 15
m_dot = 273.34078
s = 0.5
C = 2.5
T_amb = 273.15 
delta_T = 33.334





#Minor computations with input parameters
A_star = A / A_e
r_s_star = rs / np.sqrt(A_e)
N_star = N / (c / d)
d_star = d / np.sqrt(A_e)
fb = N_star * B * c / (d_star * np.sqrt(A_e))
eta = [(1 - M * np.cos(theta)) / fb * f for f in frequency] 
M_t = np.pi * N_star
M_m = max(1,M_t)
delta = M_t
delta_T_star = delta_T / T_amb
if M_t > 1.05:
    delta = M_t / (abs(1 - V / B))
    
#Inlet broadband noise    
    
i_m = 0    
j_m = 0
if delta <= 1.05:
    j_m = 1
m_dot_star = m_dot / (rho * c * A_e) 
M_x = m_dot_star / A_star
M_r = np.sqrt(M_t**2 - M_x**2)
eta_l = [int(10**(-1/20) * i) + 1 for i in eta]#Lowest harmonic number that falls within the band, calculate nodes for each frequency
eta_u = [int(10**(1/20) * i) for i in eta]       #Highest harmonic nubmer that falls within the band

n = []
for i in range(len(frequency)):
    if eta_l[i] > eta_u[i]:
        n.append([0,0])
    else:
        #node = np.arange(eta_l[i], eta_u[i] + 1, 1) 
        n.append([eta_l[i], eta_u[i]])

s_star = s / C
k_m = 0
if s_star > 1:
    k_m = 1
F = 1
if M_r > 0.9:
    F = 0.81 * M_r**-2
a_kl = np.array([[0.5,0.5],[0.5,0]])
l_m = 1
sigma = 2.2

PI_star_inlet_broadband = 1.552E-4 * s_star**-a_kl[k_m][l_m] * M_m**2 * m_dot_star / A_star * delta_T_star**2 * F

S_inlet_broadband = [0.116 * np.exp(-0.5 * (np.log(f / 2.5) / np.log(sigma))**2) for f in eta]
        
D = 10**-0.87

P2_inlet_broadband = [A_star * PI_star_inlet_broadband / (4 * np.pi * r_s_star**2) * D * S_eta / (1 - M**2 * np.cos(theta))**4 for S_eta in S_inlet_broadband]


SPL_inlet_broadband = np.array([10 * np.log10(P2) + 20 * np.log10(rho * c**2 / Pref) for P2 in P2_inlet_broadband])


#Inlet rotor-stator interaction tones
G_ij = np.array([[1,0.580], [0.625,0.205]])

a_kl = np.array([[1,1],[1,0]])
F = 2.053 * M_m**-2.31 * M_r**5
if M_r <= 0.72:
    F = 0.397 * M_m**-2.31
if 0.866 * M_m**0.462 < M_r:
    F = 0.315 * M_m**3.69 * M_r**-8

S_1ij = np.array([[0.499,0.136],[0.799,0.387]])
S_nij = np.array([[0.250,0.432],[0.101,0.307]])
S_inlet_rotor = []
count = 0
for y in range(len(frequency)):
    for z in range(len(n[y])):
        if n[y][z] == 0:
            continue
        else:
            tones = np.arange(n[y][0], n[y][1] + 1, 1)
            for x in range(len(tones)):
                if tones[x] == 1:
                    count += S_1ij[i_m,j_m]
                else:
                    count += S_nij[i_m,j_m] * 10**(-0.3*(tones[x]-2))
    S_inlet_rotor.append(count)
    count = 0        
D = 10**-0.85            
a_kl = np.array([[1,1],[1,0]])

PI_star_inlet_rotor = 2.683E-4 * G_ij[i_m,j_m] * s_star**-a_kl[k_m,l_m] * M_m**4.31 * m_dot_star / A_star * delta_T_star**2 * F

P2_inlet_rotor = [A_star * PI_star_inlet_rotor / (4 * np.pi * r_s_star**2) * D * S_eta / (1 - M**2 * np.cos(theta))**4 for S_eta in S_inlet_rotor]

SPL_inlet_rotor = np.array([10 * np.log10(P2) + 20 * np.log10(rho * c**2 / Pref) for P2 in P2_inlet_rotor])


#Inlet flow distortion tones
a_kl = np.array([[1,1],[1,0]])

S_inlet_distortion = []
count = 0
for y in range(len(frequency)):
    for z in range(len(n[y])):
        if n[y][z] == 0:
            continue
        else:
            tones = np.arange(n[y][0], n[y][1] + 1, 1)
            for x in tones:
                count += 9 * 10**-float(x)
    S_inlet_distortion.append(count)
    count = 0        

PI_star_inlet_distortion = 1.488E-4 * s_star**-a_kl[k_m,l_m] * M_m**4.31 * m_dot_star / A_star * delta_T_star**2 * F

P2_inlet_distortion = [A_star * PI_star_inlet_distortion / (4 * np.pi * r_s_star**2) * D * S_eta / (1 - M**2 * np.cos(theta))**4 for S_eta in S_inlet_distortion]

SPL_inlet_distortion = np.array([10 * np.log10(P2) + 20 * np.log10(rho * c**2 / Pref) for P2 in P2_inlet_distortion])

#Combination tone noise
G_ij = np.array([[1,1],[0.316,0.316]])

F_8 = 10**(-6.75*(1.61-M_r))
if M_r < 1:
    F_8 = 0
if M_r > 1.61:
    F_8 = 10**(-1.21*(M_r-1.61))
F_4 = 10**(-14.75*(1.322-M_r))
if M_r < 1:
    F_4 = 0
if M_r > 1.322:
    F_4 = 10**(-1.33*(M_r-1.322))
F_2 = 10**(-31.85*(1.146-M_r))
if M_r < 1:
    F_2 = 0
if M_r > 1.146:
    F_2 = 10**(-1.41*(M_r-1.146))
D = 10**-0.38

S_8 = []
S_4 = []
S_2 = []
for f in eta:
    if f <= 0.125:
        S_8.append(0.405 * (8 * f)**5)
    else:
        S_8.append(0.405 * (8 * f)**-float(-3))
for f in eta:
    if f <= 0.25:
        S_4.append(0.520 * (4 * f)**5)
    else:
        S_4.append(0.520 * (4 * f)**-float(-5))
for f in eta:
    if f <= 0.5:
        S_2.append(0.332 * (2 * f)**3)
    else:
        S_2.append(0.332 * (2 * f)**-float(-3))   
K_8 = 6.109E-4
K_4 = 2.030E-3
K_2 = 2.525E-3

PI_8 = K_8 * G_ij[i_m,j_m] * m_dot_star / A_star * delta_T_star**2 * F_8
PI_4 = K_4 * G_ij[i_m,j_m] * m_dot_star / A_star * delta_T_star**2 * F_4
PI_2 = K_2 * G_ij[i_m,j_m] * m_dot_star / A_star * delta_T_star**2 * F_2

P2_8 = [A_star * PI_8 / (4 * np.pi * r_s_star**2) * D * S_eta / (1 - M**2 * np.cos(theta))**4 for S_eta in S_8]
P2_4 = [A_star * PI_4 / (4 * np.pi * r_s_star**2) * D * S_eta / (1 - M**2 * np.cos(theta))**4 for S_eta in S_4]
P2_2 = [A_star * PI_2 / (4 * np.pi * r_s_star**2) * D * S_eta / (1 - M**2 * np.cos(theta))**4 for S_eta in S_2]
P2_combination_tone = [P2_8[i] + P2_4[i] + P2_2[i] for i in range(len(frequency))]

SPL_combination_tone = np.array([10 * np.log10(P2) + 20 * np.log10(rho * c**2 / Pref) for P2 in P2_combination_tone])


#Discharge broadband noise
a_kl = np.array([[0.5,0.5],[0.5,0]])
G_ij = np.array([[1,1],[2,2]])
F = 1
if M_r > 1:
    F = M_r**-2
D = 10**-0.04

PI_star_outlet_broadband = 3.206E-4 * G_ij[i_m,j_m] * s_star**-a_kl[k_m][l_m] * M_m**2 * m_dot_star / A_star * delta_T_star**2 * F

S_outlet_broadband = [0.116 * np.exp(-0.5 * (np.log(f / 2.5) / np.log(sigma))**2) for f in eta]
        
P2_outlet_broadband = [A_star * PI_star_outlet_broadband / (4 * np.pi * r_s_star**2) * D * S_eta / (1 - M**2 * np.cos(theta))**4 for S_eta in S_outlet_broadband]

SPL_outlet_broadband = np.array([10 * np.log10(P2) + 20 * np.log10(rho * c**2 / Pref) for P2 in P2_outlet_broadband])

#Discharge rotor-stator interaction tones
a_kl = a_kl = np.array([[1,1],[1,0]])
G_ij = np.array([[1,0.580],[2.50,0.820]])
D = 10**-0.05
S_1ij = np.array([[0.499,0.136],[0.799,0.387]])
S_nij = np.array([[0.250,0.432],[0.101,0.307]])
S_outlet_rotor = S_inlet_rotor
       
PI_star_outlet_rotor = 2.643E-4 * G_ij[i_m,j_m] * s_star**-a_kl[k_m][l_m] * M_m**2 * m_dot_star / A_star * delta_T_star**2 * F

P2_outlet_rotor = [A_star * PI_star_outlet_rotor / (4 * np.pi * r_s_star**2) * D * S_eta / (1 - M**2 * np.cos(theta))**4 for S_eta in S_outlet_rotor]

SPL_outlet_rotor = np.array([10 * np.log10(P2) + 20 * np.log10(rho * c**2 / Pref) for P2 in P2_outlet_rotor])


#Summation to obtain total fan noise
#P2_fan = [P2_inlet_broadband[i] + P2_inlet_rotor[i] + P2_inlet_distortion[i] + P2_combination_tone[i] + P2_outlet_broadband[i] + P2_outlet_rotor[i] for i in range(len(frequency))]
P2_fan = [P2_inlet_broadband[i] + P2_inlet_rotor[i] + P2_inlet_distortion[i] + P2_outlet_broadband[i] + P2_outlet_rotor[i] for i in range(len(frequency))]
P2_2fan = [2*P2_inlet_broadband[i] + 2*P2_inlet_rotor[i] + 2*P2_inlet_distortion[i] + 2*P2_outlet_broadband[i] + 2*P2_outlet_rotor[i] for i in range(len(frequency))]

SPL_fan = np.array([10 * np.log10(P2) + 20 * np.log10(rho * c**2 / Pref) for P2 in P2_fan])
SPL_2fan = np.array([10 * np.log10(P2) + 20 * np.log10(rho * c**2 / Pref) for P2 in P2_2fan])

#===================================================================================
#Engine combustion noise
#===================================================================================
"""
Ac = combuster entrance area [m^2]
mi_dot = combustor entrance mass flow rate [kg/s]
Pti = combustor entrance total pressure [N/m^2]
Ti = combustor entrance total temperature [K]
Tj = combustor exit total temperature [K]
delta_T_des = design turbine temperature extraction [K]

fp = peak frequency
"""
Ac = 0.3 #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
mi_dot = 27.334
Pti = 2437307.099
Ti = 808.3812315
Tj = 1200
delta_T_des =  691.9022561


Pinf = 101325
mi_dot_star = mi_dot / (rho * c * A_e)
Ac_star = Ac / A_e
Tj_star = Tj / T_amb
Ti_star = Ti / T_amb
Pti_star = Pti / Pinf
delta_T_des_star = delta_T_des / T_amb
fp = 400 / (1 - M * np.cos(theta))
f_over_fp = [f / fp for f in frequency]
S_combustion = []
for f in f_over_fp:
    if np.log10(f) <= -1.1:
        S_combustion.append(10**-4.3)
    elif -1.1 < np.log10(f) <= -1.0:
        S_combustion.append(10**-3.85)
    elif -1.0 <= np.log10(f) <= -0.9:
        S_combustion.append(10**-3.35)
    elif -0.9 < np.log10(f) <= -0.8:
        S_combustion.append(10**-2.85)
    elif -0.8 < np.log10(f) <= -0.7:
        S_combustion.append(10**-2.4)
    elif -0.7 < np.log10(f) <= -0.6:
        S_combustion.append(10**-2.00)
    elif -0.6 < np.log10(f) <= -0.5:
        S_combustion.append(10**-1.65)
    elif -0.5 < np.log10(f) <= -0.4:
        S_combustion.append(10**-1.35)
    elif -0.4 < np.log10(f) <= -0.3:
        S_combustion.append(10**-1.10)
    elif -0.3 < np.log10(f) <= -0.2:
        S_combustion.append(10**-0.95)
    elif -0.2 < np.log10(f) <= -0.1:
        S_combustion.append(10**-0.80)
    elif -0.1 < np.log10(f) <= 0.0:
        S_combustion.append(10**-0.80)
    elif 0.0 < np.log10(f) <= 0.1:
        S_combustion.append(10**-0.85)
    elif 0.1 < np.log10(f) <= 0.2:
        S_combustion.append(10**-0.95)
    elif 0.2 < np.log10(f) <= 0.3:
        S_combustion.append(10**-1.15)
    elif 0.3 < np.log10(f) <= 0.4:
        S_combustion.append(10**-1.40)
    elif 0.4 < np.log10(f) <= 0.5:
        S_combustion.append(10**-1.65)
    elif 0.5 < np.log10(f) <= 0.6:
        S_combustion.append(10**-1.95)
    elif 0.6 < np.log10(f) <= 0.7:
        S_combustion.append(10**-2.35)
    elif 0.7 < np.log10(f) <= 0.8:
        S_combustion.append(10**-2.70)
    elif 0.8 < np.log10(f) <= 0.9:
        S_combustion.append(10**-3.15)
    elif 0.9 < np.log10(f) <= 1.0:
        S_combustion.append(10**-3.55)
    elif 1.0 < np.log10(f) <= 1.1:
        S_combustion.append(10**-4.00)
    elif 1.1 < np.log10(f) <= 1.2:
        S_combustion.append(10**-4.40)
    elif 1.2 < np.log10(f) <= 1.3:
        S_combustion.append(10**-4.75)
    elif 1.3 < np.log10(f) <= 1.4:
        S_combustion.append(10**-5.25)
    elif 1.4 < np.log10(f) <= 1.5:
        S_combustion.append(10**-5.70)
    else:
        S_combustion.append(10**-6.20)
        
D = 10**-0.16

PI_star_combustion = 9.85E-7 * mi_dot_star / Ac_star * ((Tj_star - Ti_star) / Ti_star)**2 * Pti_star**2 * delta_T_des_star**-4

P2_star_combustion = [PI_star_combustion * Ac_star / (4 * np.pi * r_s_star**2) * D * S_eta / (1 - M * np.cos(theta))**4 for S_eta in S_combustion]
P2_2star_combustion = 2*P2_star_combustion

SPL_combustion = np.array([10 * np.log10(P2) + 20 * np.log10(rho * c**2 / Pref) for P2 in P2_star_combustion])
SPL_2combustion = np.array([10 * np.log10(P2) + 20 * np.log10(rho * c**2 / Pref) for P2 in P2_2star_combustion])

#===================================================================================
#Engine single stream circular jet noise
#===================================================================================
"""
delta_e = angle between flight vector and engine inlet axis [rad]
Aj = fully expanded jet area [m^2]
rhoj = jet density [kg/m^3]
Tj = jet total temperature [K]


w = density exponent
P = power deviation factor (deviation of accoustic power from V^8 law)
ksi = strouhal number correction factor
"""

Aj = 0.50506
rhoj =  0.7317754679
Vj = 505.26
delta_e = np.radians(3)
Tj = 402.6101202

Tj_star = Tj / T_amb
Aj_star = Aj / A_e
rhoj_star = rhoj / rho
Vj_star = Vj / c
dj_star = np.sqrt(4 * Aj_star / np.pi) 
m_theta = 1
f_star = [f * np.sqrt(A_e) / c for f in frequency]

w = 1.85
P = 10**0.365
D = 10**-0.805
ksi = 1

S_c_jet = [f * dj_star / (ksi * (Vj_star - M)) for f in f_star]
F_jet = []
for i_s in S_c_jet:
    i = np.log10(i_s)
    if i <= -2:
        F_jet.append(10**(38.2/-10))
    elif -2 < i <= -1.6:
        F_jet.append(10**(29.8/-10))
    elif -1.6 < i <= -1.3:
        F_jet.append(10**(23.5/-10))
    elif -1.3 < i <= -1.15:
        F_jet.append(10**(20.5/-10))
    elif -1.15 < i <= -1:
        F_jet.append(10**(17.8/-10))
    elif -1 < i <= -0.824:
        F_jet.append(10**(15.3/-10))
    elif -0.824 < i <= -0.699:
        F_jet.append(10**(13.5/-10))
    elif -0.699 < i <= -0.602:
        F_jet.append(10**(12.3/-10))
    elif -0.602 < i <= -0.5:
        F_jet.append(10**(11.5/-10))
    elif -0.5 < i <= -0.398:
        F_jet.append(10**(11/-10))
    elif -0.398 < i <= -0.301:
        F_jet.append(10**(10.7/-10))
    elif -0.301 < i <= -0.222:
        F_jet.append(10**(10.7/-10))
    elif -0.222 < i <= 0:
        F_jet.append(10**(11.1/-10))
    elif 0 < i <= 0.477:
        F_jet.append(10**(14.3/-10))
    elif 0.477 < i <= 1:
        F_jet.append(10**(19.5/-10))
    elif 1 < i <= 1.6:
        F_jet.append(10**(27/-10))
    else:
        F_jet.append(10**(28.25/-10))
      
        
    
PI_star_jet = 6.67E-5 * rhoj_star**w * Vj_star**8*P

P2_star_jet = [PI_star_jet * Aj_star / (4 * np.pi * r_s_star**2) * D * S_eta / (1 - M * np.cos(theta-delta_e)) * ((Vj_star - M) / Vj_star)**m_theta for S_eta in F_jet]
P2_2star_jet = 2*P2_star_jet

SPL_jet = np.array([10 * np.log10(P2) + 10 * np.log10(rho **2 * c**4 / Pref**2) for P2 in P2_star_jet])
SPL_2jet = np.array([10 * np.log10(P2) + 10 * np.log10(rho **2 * c**4 / Pref**2) for P2 in P2_2star_jet])
#===================================================================================
#Engine turbine noise
#===================================================================================
"""
A_t = turbine inlet corss-sectional area [m^2]
hti = turbine entrance total enthalpy
hsj = turbine exit static enthalpy

""" 

A_t =  1.669618547 * 0.5
#hti = ...
#hsj = ...

A_t_star = A_t / A_e
Ut_star = np.pi * N_star
R_air = 287
hti_star = 0.9 #hti / (R_air * T_amb) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
hsj_star = 0.8 #hsj / (R_air * T_amb) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#Turbine braodband noise
K = 8.589E-5
a = 1.27
b = -1.27
D = 10**0.071

S_turbine_broadband = []
for node in eta:
    no = np.log10(node)
    if no <= -0.903:
        S_turbine_broadband.append(10**-1.884)
    elif -0.903 < no < -0.796:
        S_turbine_broadband.append(10**-1.604)
    elif -0.796 < no < -0.699:
        S_turbine_broadband.append(10**-1.444)
    elif -0.699 < no < -0.602:
        S_turbine_broadband.append(10**-1.304)
    elif -0.602 < no < -0.502:
        S_turbine_broadband.append(10**-1.184)
    elif -0.502 < no < -0.398:
        S_turbine_broadband.append(10**-1.084)
    elif -0.398 < no < -0.301:
        S_turbine_broadband.append(10**-1.004)
    elif -0.301 < no < -0.201:
        S_turbine_broadband.append(10**-0.924)
    elif -0.201 < no < -0.097:
        S_turbine_broadband.append(10**-0.844)
    elif -0.097 < no < 0:
        S_turbine_broadband.append(10**-0.784)
    elif 0 < no < 0.097:
        S_turbine_broadband.append(10**-1.004)
    elif 0.097 < no < 0.204:
        S_turbine_broadband.append(10**-1.204)
    elif 0.204 < no < 0.301:
        S_turbine_broadband.append(10**-1.384)
    else:
        S_turbine_broadband.append(10**-1.924)
    
PI_star_turbine_broadband = K * ((hti_star - hsj_star) / hti_star)**a * Ut_star**b

P2_turbine_broadband = [PI_star_turbine_broadband * A_t_star / (4 * np.pi * r_s_star**2) * D * S_eta / (1 - M * np.cos(theta))**4 for S_eta in S_turbine_broadband]

SPL_turbine_broadband = np.array([10 * np.log10(P2) + 20 * np.log10(rho * c**2 / Pref) for P2 in P2_turbine_broadband])        

#Turbine pure tone noise
K = 1.162E-4
a = 1.46
b = -4.02
D = 10**-0.021

S_turbine_tone = []
count = 0
for y in range(len(frequency)):
    for z in range(len(n[y])):
        if n[y][z] == 0:
            continue
        else:
            tones = np.arange(n[y][0], n[y][1] + 1, 1)
            for x in tones:
                count += 0.6838*10**(-1*(x-1)/2)
    S_turbine_tone.append(count)
    count = 0         

PI_star_turbine_tone = K * ((hti_star - hsj_star) / hti_star)**a * Ut_star**b

P2_turbine_tone = [PI_star_turbine_tone * A_t_star / (4 * np.pi * r_s_star**2) * D * S_eta / (1 - M * np.cos(theta))**4 for S_eta in S_turbine_tone]

SPL_turbine_tone = np.array([10 * np.log10(P2) + 20 * np.log10(rho * c**2 / Pref) for P2 in P2_turbine_tone])  


P2_turbine = [P2_turbine_broadband[i] + P2_turbine_tone[i] for i in range(len(frequency))]  
P2_2turbine = [2*P2_turbine_broadband[i] + 2*P2_turbine_tone[i] for i in range(len(frequency))]     
   

SPL_turbine = np.array([10 * np.log10(P2) + 20 * np.log10(rho * c**2 / Pref) for P2 in P2_turbine])  
SPL_2turbine = np.array([10 * np.log10(P2) + 20 * np.log10(rho * c**2 / Pref) for P2 in P2_2turbine]) 

#===================================================================================
#Total engine noise
#===================================================================================

SPL_single_engine = []
count = 0
for i in range(len(frequency)):
    count = 10**(SPL_fan[i]/10) + 10**(SPL_combustion[i]/10) + 10**(SPL_jet[i]/10) + 10**(SPL_turbine[i]/10)
    SPL_single_engine.append(10*np.log10(count))
    count = 0
    
    
count = 0
SPL_both_engines = []
for i in range(len(frequency)):
    count = 10**(SPL_2fan[i]/10) + 10**(SPL_2combustion[i]/10) + 10**(SPL_2jet[i]/10) + 10**(SPL_2turbine[i]/10)
    SPL_both_engines.append(10*np.log10(count))

plt.close()
plt.figure()
#plt.plot(frequency, SPL_w_lst, label='Wing TE noise')
#plt.plot(frequency, SPL_h_lst, label='Horizontal tail noise')
#plt.plot(frequency, SPL_v_lst, label='Vertical tail noise')
#plt.plot(frequency, SPL_f_lst, label='Wing TE flap noise')
#plt.plot(frequency, SPL_nlg_tot_lst, label='Nose landing gear noise')
#plt.plot(frequency, SPL_mlg_tot_lst, label='Main landing gear noise')
#plt.plot(frequency, SPL_airframe_lst, label='Total airframe noise')

#plt.plot(frequency, SPL_inlet_broadband, label='Inlet broadband noise')
#plt.plot(frequency, SPL_inlet_rotor, label='Inlet rotor-stator interaction tones noise')            #Possibly change this to eta to shift the curves
#plt.plot(frequency, SPL_inlet_distortion, label='Inlet flow distortion noise')
##plt.plot(frequency, SPL_combination_tone, label='Combination tone noise')
#plt.plot(frequency, SPL_outlet_broadband, label='Outlet broadband noise')
#plt.plot(frequency, SPL_outlet_rotor, label='Outlet rotor-stator interaction tones noise')
plt.plot(frequency, SPL_fan, label='Fan noise')

plt.plot(frequency, SPL_combustion, label='Combustion noise')

plt.plot(frequency, SPL_jet, label='Jet noise')

#plt.plot(frequency, SPL_turbine_broadband, label='Turbine broadband noise')
#plt.plot(frequency, SPL_turbine_tone, label='Turbine tone noise')
plt.plot(frequency, SPL_turbine, label='Turbine noise')

plt.plot(frequency, SPL_single_engine, label='Total engine noise')

plt.plot(frequency, SPL_both_engines, label='Total engine noise')

plt.xscale('log')
plt.xlim([10**1.5,10**4.5])
plt.ylim([15,200])
plt.xlabel('1/3 Octave Band central frequency [Hz]')
plt.ylabel('SPL [dB]')

#plt.plot(frequency, SPL_nlg_lst, label='Nose landing gear wheel noise')
#plt.plot(frequency, SPL_mlg_lst, label='Main landing gear wheel noise')
#plt.plot(frequency, SPL_s_nlg_lst, label='Nose landing gear strut noise')
#plt.plot(frequency, SPL_s_mlg_lst, label='Main landing gear strut noise')

plt.rc('axes', labelsize=16)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=16)    # fontsize of the tick labels
plt.rc('ytick', labelsize=16)    # fontsize of the tick labels
plt.rc('legend', fontsize=16)    # legend fontsize
plt.legend()
plt.show()


print('Combination tone noise is excluded in fan noise module (recommendation)')  
print('Assumption: polar directivity angle theta of 90 degrees') 
print('Possibly change Ae reference engine area to obtain lower noise values')
print('If noise levels allow, add main landing gear contribution twice')
print('NO shock noise module due to lack of information (recommendation), no shock at the end (see propulsive design)')
