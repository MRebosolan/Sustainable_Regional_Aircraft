#input parameters
import numpy as np

MTOM = 36000 # estimate, [kg]
MTOW = MTOM * 9.81 #N
AR = 8 # estimate, [-]
half_sweep = np.cos(27 ) #estimate, [degrees]
n_max = 3 #estimate
n_ult = 1.5* n_max
wingloading = 4000 #estimate, N/m^2
powerloading = 0.44 #thrust over weight
S = MTOW /wingloading #m^2
t_over_c = 0.1 #estimate, []
taper = 0.4 #estimate, []
mach_h = 0.5 #estimate, []
V_dive = 300 #estimate, knots
lf = 
hf =
A_inlet =
ln =
p2 =
W_zfw =
b =
t_r=
wf = #max fuselage width
S_fgs = #fuselage gross shell area
lh =
T_TO = MTOW*powerloading
# Kgr =
# Ag, Bg, Cg, Dg = #look at LG_weight_estimation.py for coefficients
V_pax =
lpax = 20 # estimate, meters
Npax = 75
N_fdc = 2 #probably, pilots
N_cc = 2 #probably, cabin crew
P_c =
Sff =
N_eng = 2
N_t =
K_fsp =
W_fuel =


 
W_pax=93 #includes luggage kg
W_cargo=1000 #kg
n_crew= N_fdc+N_cc
W_payload=Npax*W_pax+W_cargo
Design_range=2000#[km]