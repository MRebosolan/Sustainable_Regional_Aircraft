#input parameters
import numpy as np
import Tail_weight_estimation as tail

# MTOM = 36000 # estimate, [kg]
# MTOW = MTOM * 9.81 #N
AR = 8 # estimate, [-]
half_sweep = np.cos(27 ) #estimate, [degrees]
n_max = 3 #estimate
n_ult = 1.5* n_max
wingloading = 4000 #estimate, N/m^2
powerloading = 0.44 #thrust over weight
# S = MTOW /wingloading #m^2
t_over_c = 0.1 #estimate, []
taper = 0.4 #estimate, []
mach_h = 0.5 #estimate, []
V_dive = 300 #estimate, knots
lf = 30 #estimate, lil shorter than CRJ as 5 seat rows are used
hf = 2.5 #estimate
A_inlet = 1.17 #m
ln = 0.8129 #m 1/4 of CRJ engine length
# b = (S * AR)**0.5 #wingspan [m]
t_r= 1.2 #maximum thickness at root [m] #bullshit estimation
widthf = 4.24 #m max fuselage width
S_fgs = widthf * np.pi * lf *0.9 #fuselage gross shell area
lh = 15 #very random estimate

Kgr = 1.08
V_pax = 282.391 #m^3
lpax = 20 # estimate, meters
Npax = 75
N_fdc = 2 #probably, pilots
N_cc = 2 #probably, cabin crew
P_c =74682.5 #Pa
Sff = 7.6
N_eng = 2
N_t = 2
rho_hydrogen = 71 #g/l
H_to_ker_ratio = 0
Sh = 20.75
half_chord_sweep_hor = 20 #deg
half_chord_sweep_vert = 35 #deg

Sv = 13.36
bv = 7.57 #m good ol' CRJ700
bh = 8.54 #estimation, m
zh = bv * 0.95

LD_c = 15
LD_c2 = 17
LD_loiter = 17

V_c=230.3
V_c2=0.8*V_c
V_loiter=0.6*V_c

cj_ck = 1.6 * 10 ** (-5)  # kerosene cj
cj_c = cj_ck * 0.349 * H_to_ker_ratio + cj_ck * (1 - H_to_ker_ratio)

cj_ck2 = cj_ck * 0.9
cj_c2 = cj_ck2 * 0.349 * H_to_ker_ratio + cj_ck2 * (1 - H_to_ker_ratio)
cj_kloiter = cj_ck*0.7
cj_loiter = cj_kloiter * 0.349 * H_to_ker_ratio + cj_kloiter * (1 - H_to_ker_ratio)

W_pax=93 #includes luggage kg
W_cargo=1000 #kg
n_crew= N_fdc+N_cc
W_payload=Npax*W_pax+W_cargo
Design_range=2000#[km]
