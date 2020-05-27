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
W_zfw = #zero f
b = #wingspan [m]
t_r=  #maximum thickness at root [m]
wf = #max fuselage width
S_fgs = #fuselage gross shell area
lh =
T_TO = MTOW*powerloading
Kgr = 1.08
V_pax =
lpax = 20 # estimate, meters
Npax = 75
N_fdc = 2 #probably, pilots
N_cc = 2 #probably, cabin crew
P_c =
Sff =
N_eng = 2
N_t =
rho_hydrogen = 33.2307 #lb/gal
W_fuel =
H_to_ker_ratio =

LD_c = 15
LD_c2 = 17
LD_loiter = 17

V_c=230.3
V_c2=0.8*V_c
V_loiter=0.6*V_c

W_pax=93 #includes luggage kg
W_cargo=1000 #kg
n_crew= N_fdc+N_cc
W_payload=Npax*W_pax+W_cargo
Design_range=2000#[km]

hydrogenratio = 0