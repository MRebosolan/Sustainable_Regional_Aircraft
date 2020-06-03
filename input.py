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
def atmosphere_calculator(h):
    T_grad = -0.0065
    T = 288.15 + T_grad * h
    P = 101325 * (T / 288.15) ** (-9.81 / (T_grad * 287))
    rho = P / (T * 287)
    a = (1.4 * 287 * T) ** 0.5
    return (T, P, rho, a)


MTOM = 32846.208787 # [kg]  #maximum takeoff mass calculated in class 2
print('please use MTOM from class 2')
OEW = 22981.846450 #kg calculated in class 2
print('please use OEW from class 2')
MTOW = MTOM * g #N  !!!!!!!!!!!!!make sure this is in newtons!!!!!!!!!!!!!
MLW = 25000 * g  # maximum landing weight [N], to be calculated
AR = 8  # estimate, [-], Aspect ratio
half_sweep = np.cos(radians(27))  # estimate, [degrees], sweep at half chord for main wing
LE_sweep = np.radians(25)         #Leading egde wing sweepTBD, assumed backward let Rick know if it turns out to be forward sweep, calculate with Adsee formula 

n_max = 2.5  # from envelope, update manually, max loading factor
n_ult = 1.5 * n_max #ultimate loading factor
wingloading = 4375.84  # estimate, N/m^2
powerloading = 0.44  # thrust over weight
S = MTOW /wingloading #m^2, wing area
print('please use S from class 2')
Tto = powerloading * MTOW #thrust at takeoff in newtons
print('please use Tto from class 2')
t_over_c = 0.1  # estimate, [] , maximum thickness over chord ratio for main wing

mach_h = 0.5  # estimate, [] #max Mach at SL
rho = 1.225   # estimate, in kg/m3
# rho_zero = rho  # kg/m3


Cr = 5.5                                 #Wing root chord [m]
Ct = 1.2                               #Wing tip chord [m]

taper = Ct / Cr                         #wing taper ratio [-]

Dfus = 2.6, "if you encounter an error here, make your program dependent on a different variable"       #Fuselage diameter [m]



Cla_aileron = 6.4                     #1/rad, sectional lift curve slope at wing section where aileron is located, determine by datcom method or airfoil simulation
Cd0_aileron = 0.0039                     #zero drag coefficient [-] at wing section where aileron is located, determine by airfoil simulation


lf = 30  #lenght of fuselage m estimate, lil shorter than CRJ as 5 seat rows are used
hf = 2.5  # height of fuselage estimate
A_inlet = 1.17  # m2, engine inlet area
ln = 0.8129  # m 1/4 of CRJ engine length, length of nacelle
b = (S * AR)**0.5 #wingspan [m]
print('please use b from class 2')
t_r = 1.0  # maximum thickness at root [m] #bullshit estimation
widthf = 2.8  # m max fuselage width
A_fuselage = np.pi*widthf*hf
ellipse_fuselage = 2*np.pi * (((widthf/2)**2 + (hf/2)**2)/2)**0.5
S_fgs = ellipse_fuselage * lf * 0.9  # fuselage gross shell area, APPROXIMATION
lh = 15  # very random estimate, distance between wing and tail aerodynamic centers
lv = 16 #very random estimate, distance between wing and vertical tail aerodynamic centers
x_ac = 12           #x location of wing aerodynamic center measured from the nose of the aircraft, TBD
x_apu = 20            #cg location of the apu measured from the nose of the aircraft [m], TBD


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
V_C=Envelope.V_C #KNOTS   cruise speed         #Velocities from flight envelope, ask George
V_S=Envelope.V_S #KNOTS stall speed
V_S2=Envelope.V_S2 #KNOTS  stall speed for negative cl max
V_dive=Envelope.V_D #KNOTS  dive speed
V_A=Envelope.V_A #KNOTS  maximum gust intensity speed
V_B=Envelope.V_B #KNOTS   flaps deflected gust speed



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
x_first_pax = 7.5   #x-location measured from the nose [m] where first passenger row is located
seat_pitch = 30     #Seat pitch [inch]!!!!!!!!!!!!!!!!
n_rows = 15         #Number of passenger rows [-] (=n_pax/n_seatsabreast)
W_cargo = 1000  # kg #Extra cargo weight
cargo_fwd_fraction = 1/3 #estimate, amount of cargo in fwd hold
cargo_aft_fraction = 2./3 #estimate, amount of cargo in aft hold
n_crew = N_fdc + N_cc
W_payload = Npax * W_pax + W_cargo
Design_range = 2000  # [km]
hydrogen_cost = 2.4  # US DOLLARS per KG

H_to_ker_ratio = 1

# Flight performance
rho0 = rho  # kg/m^3
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
# rho_c = 0.4135  # [kg/m^3], cruise density (this is the one for 10 km cruise altitude)
v_approach = 66 # m/s, RICK FIX THIS
mach_app = v_approach/340.3 # RICK FIX THIS
V_to = 1.05 * ((MTOW/S)*(2/1.225)*(1/CLmax_to))**0.5 #takeoff speed



Cruise_alt = 10  # Max operating altitude in km

T, P, rho_c, a = atmosphere_calculator(Cruise_alt*1000)


V_C_TAS = V_C * 0.514444444 / ((rho_c/rho)**0.5) #m/s cruise speed, according to flight envelope, 
V_C_estimate = 230 #m/s, design parameter
mach_cruise = V_C_estimate/a

# parameters for Carbon Footprint
Range_CRJ = 2593  # design range
Pax_CRJ = 78  # Number of passengers
Fuel_use_CRJ = 4740  # Fuel mass at design range
Cruise_alt_max_CRJ = 12497  # Max operating altitude


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

# Aerodynamics:
x_start_Cr = 12                 #x location where root chord starts, measured from the nose of the aircraft [m], TBD


MAC =  2 / 3 * Cr * ((1 + taper + taper**2) / (1 + taper)) #length of mean aerodynamic chord, formula taken from Adsee II
y_MAC = b / 6 * ((1 + 2 * taper) / (1 + taper))             #spanwise location of mean aerodynamic chord
x_lemac_rootchord = y_MAC * np.tan(LE_sweep)                          #x position of mac at leading edge [m], measured from the start of the root choord!!!!
x_LEMAC_nose = x_start_Cr + x_lemac_rootchord

# Aerodynamics for scissor plot:

cl0 = 0.153333 #preliminary estimate 
cm0 = -0.018 #preliminary estimate


tail_speedratio = 1**0.5 # sead, T tail

