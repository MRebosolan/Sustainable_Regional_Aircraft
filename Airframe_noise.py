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
bw = cl2.b
Aw = cl2.S
Ah = 30
Av = 15
Af = 20
flap_defl = np.radians(30)
bh = 11
bv = 6
bf = 15
mu_inf = 1.84E-5
frequency = np.arange(10,100010,10)
U = input.v_approach
M = input.mach_app
c = U / M
rho = input.rho0
theta = np.radians(90)
phi = np.radians(5)
n_nlg = 2
d_nlg = input.d_wheel_nose_lg
l_nlg = input.strut_length_nose_lg
n_mlg = 4
d_mlg = input.d_wheel_main_lg
l_mlg = input.strut_length_main_lg

Pref = 2E-5
rs = 100 
rs_star = rs / bw

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
    if S_w[i] < 2:
        F_f.append(0.0480 * S_f[i])
    elif S_w[i] > 20:
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

plt.close()
plt.figure()
plt.plot(frequency, SPL_w_lst, label='Wing TE noise')
plt.plot(frequency, SPL_h_lst, label='Horizontal tail noise')
plt.plot(frequency, SPL_v_lst, label='Vertical tail noise')
#plt.plot(frequency, SPL_nlg_lst, label='Nose landing gear wheel noise')
#plt.plot(frequency, SPL_mlg_lst, label='Main landing gear wheel noise')
#plt.plot(frequency, SPL_s_nlg_lst, label='Nose landing gear strut noise')
#plt.plot(frequency, SPL_s_mlg_lst, label='Main landing gear strut noise')
plt.plot(frequency, SPL_nlg_tot_lst, label='Nose landing gear noise')
plt.plot(frequency, SPL_mlg_tot_lst, label='Main landing gear noise')
plt.plot(frequency, SPL_airframe_lst, label='Total airframe noise')
plt.xscale('log')

plt.legend()
plt.show


    