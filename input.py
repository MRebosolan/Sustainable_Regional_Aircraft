# -*- coding: utf-8 -*-
"""
Created on Thu May 28 09:48:30 2020

@author: malfl
"""

# input parameters
import numpy as np
import Envelope
# import Class_2_estimation as cl2

g = 9.80665  # [m/s^2]
from math import radians

MTOM = 32846.208787 # [kg]  #maximum takeoff mass calculated in class 2
OEW = 22981.846450 #kg calculated in class 2
MTOW = MTOM * g #N  !!!!!!!!!!!!!make sure this is in newtons!!!!!!!!!!!!!
MLW = 25000 * g  # maximum landing weight [N], to be calculated
AR = 8  # estimate, [-], Aspect ratio
half_sweep = np.cos(radians(27))  # estimate, [degrees], sweep at half chord for main wing

n_max = 2.5  # from envelope, update manually, max loading factor
n_ult = 1.5 * n_max #ultimate loading factor
wingloading = 4375.84  # estimate, N/m^2
powerloading = 0.44  # thrust over weight
S = MTOW /wingloading #m^2, wing area
Tto = powerloading * MTOW #thrust at takeoff in newtons
t_over_c = 0.1  # estimate, [] , maximum thickness over chord ratio for main wing
taper = 0.4  # estimate, []
mach_h = 0.5  # estimate, [] #max Mach at SL
rho = 1.225 * 0.0624279606  # estimate, in lbs/ft3
rho_zero = 0.00237  # fucking americans, this is slug/ft3

lf = 30  #lenght of fuselage m estimate, lil shorter than CRJ as 5 seat rows are used
hf = 2.5  # height of fuselage estimate
A_inlet = 1.17  # m2, engine inlet area
ln = 0.8129  # m 1/4 of CRJ engine length, length of nacelle
b = (S * AR)**0.5 #wingspan [m]
t_r = 1.0  # maximum thickness at root [m] #bullshit estimation
widthf = 4.24  # m max fuselage width
S_fgs = widthf * np.pi * lf * 0.9  # fuselage gross shell area
lh = 15  # very random estimate, distance between wing and tail aerodynamic centers

Kgr = 1.08  # constant for the gear, torenbeek parameter
V_pax = 282.391  # m^3, cabin volume
lpax = 20  # estimate, meters, cabin length
Npax = 75 #number of passengers
N_fdc = 2  # probably, pilots
N_cc = 2  # probably, cabin crew
P_c = 74682.5  # Pa, design cabin pressure
Sff = 7.6  # freight floor area estimate
N_eng = 2  # number of engines
N_t = 2  # number of fuel tanks
rho_hydrogen = 70  # g/l
K_fsp = 0.820  # kg/l, jet A, jet fuel density
H_to_ker_ratio = 0  # fuck hydrogen atm yo, hydrogen to kerosene ratio
Sh = 20.75  # m2 crj700 shizzle yo, horizontal tail area
half_chord_sweep_hor = np.radians(20)  # deg, sweep at half chord of horizontal tail
half_chord_sweep_vert = np.radians(35)  # deg, sweep at half chord of vertical tail

Sv = 13.36  # m2 crj700 shizzle yo
bv = 7.57  # m good ol' CRJ700
bh = 8.54  # estimation, m
zh = bv * 0.95

#Lift/drag ratios
LD_c = 15
LD_c2 = 17
LD_loiter = 17


#Flight envelope
V_C=Envelope.V_C #KNOTS            #Velocities from flight envelope, ask George
V_S=Envelope.V_S #KNOTS
V_S2=Envelope.V_S2 #KNOTS
V_dive=Envelope.V_D #KNOTS
V_A=Envelope.V_A #KNOTS
V_B=Envelope.V_B #KNOTS

nlim=Envelope.nlimpos

V_C2 = 0.8 * V_C
V_loiter = 0.6 * V_C


#Class 1 weight estimation
cj_ck = 1.6 * 10 ** (-5)  # kerosene cj               #Parameters about class 1 weight estimation, ask Jari
cj_c = cj_ck * 0.349 * H_to_ker_ratio + cj_ck * (1 - H_to_ker_ratio)

cj_ck2 = cj_ck * 0.9
cj_c2 = cj_ck2 * 0.349 * H_to_ker_ratio + cj_ck2 * (1 - H_to_ker_ratio)

cj_kloiter = cj_ck * 0.7
cj_loiter = cj_kloiter * 0.349 * H_to_ker_ratio + cj_kloiter * (1 - H_to_ker_ratio)
t_loiter = 2700  # s, as in 45 minutes




W_pax = 93  # total weight per passenger, includes luggage kg
W_cargo = 1000  # kg #Extra cargo weight
n_crew = N_fdc + N_cc
W_payload = Npax * W_pax + W_cargo
Design_range = 2000  # [km]
hydrogen_cost = 2.4  # US DOLLARS per KG

H_to_ker_ratio = 1

# Flight performance
rho0 = 1.225  # kg/m^3
CD0 = 0.01277  # [-], to be refined (roskam) DONE IN MIDTERM, TALK TO JORN
CD0_togd = 0.01277 + .015 + .02  # [-], to be refined as this comes from roskam statistics
CD0_landGD = CD0 + .02 + .065  # [-], to be refined as this comes from roskam statistics
e = 0.85  # [-], Oswald effiency factor, to be refined as this comes from roskam statistics
CLmax_land = 2.25  # TBD
CLmax_clean = 1.8  # TBD
CLmax_to = 2.1  # TBD
mu = 0.04  # runway friction coefficient at take-off, to be reconsidered
mu_br = 0.3  # braking coefficient during landing, to be reconsidered
h_sc = 50 * 0.3048  # screen height equal to 50 ft [m]
gamma_cl = np.radians(7)  # climb angle right after rotation, to be refined [rad]
gamma_ap = np.radians(3)  # approach angle (glide slope) [rad]
Trev = 50000  # [N], maximum thrust reverse force applied during braking
c_t = 0.0002  # [1/s] specific fuel consumption, to be refined
H = 120E6  # Heating value of hydrogen, refine if we fly on kerosene and hydrogen simultenously, or 141.7E6 (higher value of hydrogen)
rho_c = 0.4135  # [kg/m^3], cruise density (this is the one for 10 km cruise altitude)

# parameters for Carbon Footprint
Range_CRJ = 2593  # design range
Pax_CRJ = 78  # Number of passengers
Fuel_use_CRJ = 4740  # Fuel mass at design range
Cruise_alt_max_CRJ = 12497  # Max operating altitude
Cruise_alt = 10  # Max operating altitude in km

# H2 NOx emission: Depends on engine characteristics
A = 14  # Correlation constant for emission index based on Jet-A fuel (advanced LDI tech as reference)
eq = 0.4  # equivalence ratio (fuel/air // fuel/air stoichiometric)
fa_st = 1. / 34.33  # stoichiometric fuel/air ratio for H2
fa = eq * fa_st  # actual fuel/air ratio
P3 = 0.7  # fuel injector inlet pressure MPA
T3 = 800  # fuel injector inlet temperature 600 K approach, 700 K cruise, 800K take-off
dPP = 5  # dP/P fuel injector air flow pressure drop ratio
# kg NOx/ kg fuel
NOx_H2 = A * P3 ** 0.594 * np.exp(T3 / 350) * fa ** 1.6876 * (100 * dPP) ** -0.56 / 1000

# GWP
# Global Warming Potential (equivalent emission ratio to CO2 on 100 year scale (CO2-eq))
# For CO2, H20, Nox
# Altitudes 0 to 15 km
GWP_alt = np.array([[1., 0., -7.1],
                    [1., 0., -7.1],
                    [1., 0., -7.1],
                    [1., 0., -4.3],
                    [1., 0., -1.5],
                    [1., 0., 6.5],
                    [1., 0., 14.5],
                    [1., 0., 37.5],
                    [1., 0., 60.5],
                    [1., 0, 64.7],
                    [1., 0.34, 57.7],
                    [1., 0.43, 46.5],
                    [1., 0.53, 25.6],
                    [1., 0.62, 4.6],
                    [1., 0.72, 0.6]])
GWP = GWP_alt[Cruise_alt]  # Specify the altitude in km
