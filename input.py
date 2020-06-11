import numpy as np
import Envelope
from math import radians
# import Class_2_estimation as cl2

#------------------------------------------------------------------------------------------------------------------
rho = 1.225                       # estimate, in kg/m3
rho0 = rho                        # kg/m^3
Cruise_alt = 10                   # Max operating altitude in km
g = 9.80665                       # [m/s^2]
Design_range = 2000               # [km]
#------------------------------------------------------------------------------------------------------------------

mach_h = 0.5                      # estimate, [] #max Mach at sea-level

              
def atmosphere_calculator(h):
    T_grad = -0.0065
    T = 288.15 + T_grad * h
    P = 101325 * (T / 288.15) ** (-9.81 / (T_grad * 287))
    rho = P / (T * 287)
    a = (1.4 * 287 * T) ** 0.5
    return (T, P, rho, a)
T, P, rho_c, a = atmosphere_calculator(Cruise_alt*1000)



#Masses and weights
#------------------------------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------------------------------
MTOM = 28083                      # [kg]  #maximum takeoff mass calculated in class 2
print('please use MTOM from class 2')
OEW = 18332                # kg calculated in class 2
print('please use OEW from class 2')
MTOW = MTOM * g                   # [N] Make sure it is in Newtons!
MLW = MTOW                   # maximum landing weight [N], same as MTOW because we burn low fuel



#Wing and Empennage parameters
#------------------------------------------------------------------------------------------------------------------
t_over_c = 0.14                  # estimate, [] , maximum thickness over chord ratio for main wing
AR = 8                            # estimate, [-], Aspect ratio

half_sweep = 0.3985698252407915         # estimate, [degrees], sweep at half chord for main wing
LE_sweep = 0.5051251603789663         # Leading egde wing sweepTBD, assumed backward let Rick know if it turns out to be forward sweep, calculate with Adsee formula
quarter_sweep = 0.45322773323204535
wingloading = 4375.84             # [N/m2] estimate

wingloading = 4375.84             # [N/m2] estimate 
AR_h = 4                          #Aspect ratio of the horizontal tail [-], TBD
AR_v = 1.7                        #AR vtail; 1< AR_v<2
taper_v = 0.6                     #taper ratio vertical tail
e_tail = 0.85                     #oswald efficiency factor of tail, TBD
#------------------------------------------------------------------------------------------------------------------
# half_sweep = 0.663         # estimate, [degrees], sweep at half chord for main wing
# LE_sweep =0.7485         # Leading egde wing sweepTBD,
# quarter_sweep = 0.707


S = MTOW /wingloading             # [m2] wing area
print('please use S from class 2')

Cr = 4.292998                     # Wing root chord [m]
Ct = 1.328058                     # Wing tip chord [m]
taper = Ct / Cr                   # wing taper ratio [-]
b = (S * AR)**0.5                 # wingspan [m]

print('please use b from class 2')

Sv = 13.36                        # [m2] CRJ700 | Obtain realistic value from Vtail area sizing
bv = 7.57                         # [m] vertical tail span CRJ700
zh = bv * 0.95                    # Height of horizontal stabilizer measured from the bottom of the vertical tail [m]
Sh_over_S = 0.3725
Sh = Sh_over_S * S                        # m2 crj700 shizzle yo, horizontal tail area
half_chord_sweep_hor = np.radians(20)   # deg, sweep at half chord of horizontal tail

bh = (AR_h*Sh)   **0.5                     # [m] Horizontal tail span 
half_chord_sweep_vert = np.radians(35)  # deg, sweep at half chord of vertical tail







#Fuselage, cabin and loading parameters

#------------------------------------------------------------------------------------------------------------------

hf = 3.486                        # height of fuselage estimate
widthf = 3.486                      # m max fuselage width
lpax = 17.76                        # estimate, meters, cabin length, based on a seat pitch of 32 inch
seat_pitch = 32                   # Seat pitch [inch]
Npax = 75                           # number of passengers
pax_abreast = 3+2                   #seating abreast configuration
N_fdc = 2                         # number of pilots [probably]
N_cc = 2                          # number of cabin crew [probably]
P_c = 74682.5                     # Pa, design cabin pressure, [N/m^2]
P_cruise = 22632.1                # Pa, outside pressure at cruise
Sff = 7.6                         # freight floor area estimate
W_pax = 93                        # total weight per passenger, includes luggage kg
seatlength = 21 * 0.0254 *0.5      #cg of first row of people
x_first_pax = 6.9+ seatlength            # x-location measured from the nose [m] where first passenger row is located
n_rows = 15                       # Number of passenger rows [-] (=n_pax/n_seatsabreast)
W_cargo = 1000                    # kg #Extra cargo weight
n_crew = N_fdc + N_cc             # Total amount of crew
W_payload = Npax * W_pax + W_cargo
cockpit_length=3.6                #Bit smaller than A220
A_fuselage = np.pi*widthf*0.5*hf*0.5 # Area of the fuselage in m^2
ellipse_fuselage = 2*np.pi * (((widthf/2)**2 + (hf/2)**2)/2)**0.5  #circumference of the fuselage [m]
#------------------------------------------------------------------------------------------------------------------
lf = 28                           # length of fuselage m estimate, lil shorter than CRJ as 5 seat rows are used, ADD TANK CYLINDER LENGTH
S_fgs = ellipse_fuselage * lf * 0.9  # fuselage gross shell area, APPROXIMATION

V_pax = A_fuselage*lpax                 # m^3, cabin volume, can be estimated with a small calculation
cargo_fwd_fraction = 1/3          # estimate, amount of cargo in fwd hold
cargo_aft_fraction = 2./3         # estimate, amount of cargo in aft hold
x_cg_fwd_cargo = 6                # cg of forward cargo compartment, measured from the aircraft nose [m]
x_cg_aft_cargo = 22               # cg of aft cargo compartment, measured from the aircraft nose [m]


#Propulsion related parameters
#------------------------------------------------------------------------------------------------------------------
Kgr = 1.08                        # constant for the gear, torenbeek parameter
N_eng = 2                         # number of engines
rho_hydrogen = 70                 # g/l
K_fsp = 0.820                     # kg/l, jet A, jet fuel density
N_t = 0                           # number of kerosene fuel tanks
H_to_ker_ratio = 1                # hydrogen to kerosene ratio

#------------------------------------------------------------------------------------------------------------------

powerloading = 0.44               # thrust over weight from loading diagram, update based on flight performance take off length
Tto = powerloading * MTOW         #thrust at takeoff in newtons
print('please use Tto from class 2')
A_inlet = 1.58                    # m2, engine inlet area
ln = 0.8129                       # m 1/4 of CRJ engine length, length of nacelle
bn = 1.2                          #maximum width of engine [m]
z_engine = -bn/2 - 0.5             # vertical placement of engine w.r.t. wing root


#Parameters related to emissions and cost
#------------------------------------------------------------------------------------------------------------------
hydrogen_cost = 2.4               # US$ per kg
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
GWP = GWP_alt[Cruise_alt]          # Specify the altitude in km
#------------------------------------------------------------------------------------------------------------------



# H2 NOx emission: Depends on engine characteristics
#------------------------------------------------------------------------------------------------------------------
# parameters for Carbon Footprint
Range_CRJ = 2593                   # design range
Pax_CRJ = 78                       # Number of passengers
Fuel_use_CRJ = 4740                # Fuel mass at design range
Cruise_alt_max_CRJ = 12497         # Max operating altitude
#------------------------------------------------------------------------------------------------------------------
A = 14                             # Correlation constant for emission index based on Jet-A fuel (advanced LDI tech as reference)
eq = 0.4                           # equivalence ratio (fuel/air // fuel/air stoichiometric)
fa_st = 1. / 34.33                 # stoichiometric fuel/air ratio for H2
fa = eq * fa_st                    # actual fuel/air ratio
P3 = 0.7                           # fuel injector inlet pressure MPA
T3 = 800                           # fuel injector inlet temperature 600 K approach, 700 K cruise, 800K take-off
dPP = 5                            # dP/P fuel injector air flow pressure drop ratio
# kg NOx/ kg fuel
NOx_H2 = A * P3 ** 0.594 * np.exp(T3 / 350) * fa ** 1.6876 * (100 * dPP) ** -0.56 / 1000




#Parameters related to flight performance
#------------------------------------------------------------------------------------------------------------------
gamma_ap = np.radians(3)           # approach angle (glide slope) [rad]
gamma_cl = np.radians(7)           # climb angle right after rotation, to be refined [rad]
H = 120E6                          # Heating value of hydrogen, or 141.7E6 (higher value of hydrogen)
rho_c = 0.4135                     # [kg/m^3], cruise density (this is the one for 10 km cruise altitude)                   # TBD
mu = 0.04                          # runway friction coefficient at take-off, to be reconsidered
mu_br = 0.3                        # braking coefficient during landing, to be reconsidered
h_sc = 50 * 0.3048                 # screen height equal to 50 ft [m]
e = 0.85                           # [-], Oswald effiency factor, to be refined as this comes from roskam statistics
#------------------------------------------------------------------------------------------------------------------
CD0 = 0.01277                      # [-], to be refined (roskam) DONE IN MIDTERM, TALK TO JORN
CD0_togd = 0.01277 + .015 + .02    # [-], to be refined as this comes from roskam statistics
CD0_landGD = CD0 + .02 + .065      # [-], to be refined as this comes from roskam statistics
Trev = 50000                       # [N], maximum thrust reverse force applied during braking
c_t = 0.0002                       # [1/s] specific fuel consumption, to be refined

v_approach = 81.13                # m/s, update from flight performance
mach_app = v_approach/340.3        # 




#Parameters regarding aerodynamics
#------------------------------------------------------------------------------------------------------------------
#Airfoil
#CleanConfiguration
Clmax_clean   = 2.38
Cl0_clean     = 0.5
alpha0_clean  = -3.936
Cldes_clean   = 0.555368
Cd0_clean 	  = 0.007
Cla_clean 	  = 6.28
Re_clean	  = 19500000

#Take-off Configuration
Clmax_TO  = 2.38
Re_TO	  = 15071059

#Landing Configuration
Clmax_Lnd  = 2.38
Re_Lnd     = 15071059

#Wing 
#Clean Configuration
CLmax_clean  = 1.8
CL0_clean    = 0.405
Alpha0_clean = -3.936
CLdes_clean  = 0.44888
CD0_clean    = 0 
CLa_clean    = 0.09858

#Take-off Configuration             
CLmax_to  = 2.1
CL0_to    = 0.763
Alpha0_to = -7.61995                

#Landing Configuration
CLmax_land  = 2.25
CL0_land    = 0.946
Alpha0_land = -9.4619

#------------------------------------------------------------------------------------------------------------------

V_to = 1.05 * ((MTOW/S)*(2/1.225)*(1/CLmax_to))**0.5 #takeoff speed
t_r = t_over_c * Cr                          # maximum thickness at root [m] #bullshit estimation
Cla_aileron = 6.48                 #1/rad, sectional lift curve slope at wing section where aileron is located, determine by datcom method or airfoil simulation
Cd0_aileron = 0.007               #zero drag coefficient [-] at wing section where aileron is located, determine by airfoil simulation
#------------------------------------------------------------------------------------------------------------------
x_start_Cr = 10.7                    # x-location where root chord starts, measured from the nose of the aircraft [m], TBD
MAC =  2 / 3 * Cr * ((1 + taper + taper**2) / (1 + taper)) #length of mean aerodynamic chord, formula taken from Adsee II
y_MAC = b / 6 * ((1 + 2 * taper) / (1 + taper))            #spanwise location of mean aerodynamic chord
x_lemac_rootchord = y_MAC * np.tan(LE_sweep)               #x position of mac at leading edge [m], measured from the start of the root choord!!!!
x_LEMAC_nose = x_start_Cr + x_lemac_rootchord

Cla_aileron = 6.48                  #1/rad, sectional lift curve slope at wing section where aileron is located, determine by datcom method or airfoil simulation
Cd0_aileron = 0.007               #zero drag coefficient [-] at wing section where aileron is located, determine by airfoil simulation


# Aerodynamics for scissor plot:
#------------------------------------------------------------------------------------------------------------------
tail_speedratio = 1**0.5           # SEAD, T tail
cl0 = 0.488                   # preliminary estimate, TBD from airfoil analysis
CL0takeoff = 0.8                #zero angle lift at takeoff condition
CL0land = 0.95
cm0 = -0.119                    # preliminary estimate, TBD from airfoil analysis
zero_lift_angle = np.radians(-3.941)# degrees, PRELIMINARY estimate, TBD from airfoil analysis
cl_htail_max  = -0.8                #estimate coming from sead: maximum lift coefficient of tail
horizontal_margin = 0.15
#------------------------------------------------------------------------------------------------------------------

z_position_wing = hf - 0.6         # m, PRELIMINARY, still requires thought, for downwash calc
z_position_horizontal = zh + hf    # where tail is positioned, for downwash calc

z_cg = 0.5*hf  

taper_h = 0.5                   #approximation of taper ratio of horizontal tail
sweep_LE_H = np.arctan(np.tan(half_chord_sweep_hor) - 4 / AR_h * ((0 - 25) / 100 * (1 - taper_h) / (1 + taper_h)))

y_MAC_h = bh / 6 * ((1 + 2 * taper_h) / (1 + taper_h))            #spanwise location of mean aerodynamic chord
x_lemac_rootchord_h = y_MAC_h * np.tan(sweep_LE_H)               #x position of mac at leading edge [m], measured from the start of the root choord!!!!
c_root_h = 2*Sh / ((1 + taper_h)*bh)
c_tip_h = taper_h * c_root_h
c_mac_h = 2/3 * c_root_h * (1 + taper_h + taper_h**2)/(1 + taper_h)
x_rootchord_h = lf -c_root_h - 0.4



#Parameters regarding Class I      # Parameters about class 1 weight estimation, ask Jari
#------------------------------------------------------------------------------------------------------------------
cj_ck = 1.6 * 10 ** (-5)           # kerosene cj of a high bypass engine  
cj_c=(H_to_ker_ratio-1)*0.349*cj_ck/((H_to_ker_ratio-1)*0.349-H_to_ker_ratio)+0.349*(cj_ck-(H_to_ker_ratio-1)*0.349*cj_ck/((H_to_ker_ratio-1)*0.349-H_to_ker_ratio))


cj_ck2=cj_ck*0.9
cj_c2=(H_to_ker_ratio-1)*0.349*cj_ck2/((H_to_ker_ratio-1)*0.349-H_to_ker_ratio)+0.349*(cj_ck2-(H_to_ker_ratio-1)*0.349*cj_ck2/((H_to_ker_ratio-1)*0.349-H_to_ker_ratio))

cj_kloiter=cj_ck*0.7
cj_loiter=(H_to_ker_ratio-1)*0.349*cj_kloiter/((H_to_ker_ratio-1)*0.349-H_to_ker_ratio)+0.349*(cj_kloiter-(H_to_ker_ratio-1)*0.349*cj_kloiter/((H_to_ker_ratio-1)*0.349-H_to_ker_ratio))
t_loiter = 2700                    # s, as in 45 minutes
#------------------------------------------------------------------------------------------------------------------


#Parameters regarding flight envelope #Velocities from flight envelope, ask George
#------------------------------------------------------------------------------------------------------------------
V_C=Envelope.V_C                   # KNOTS  cruise speed       
V_S=Envelope.V_S                   # KNOTS  stall speed
V_S2=Envelope.V_S2                 # KNOTS  stall speed for negative cl max
V_dive=Envelope.V_D                # KNOTS  dive speed
V_A=Envelope.V_A                   # KNOTS  maximum gust intensity speed
V_B=Envelope.V_B                   # KNOTS  flaps deflected gust speed
nlim=Envelope.nlimpos              #limit load factor [-]
nmax = Envelope.nmax               #maximum load factor [-]
n_max = nmax                        # from envelope, update manually, max loading factor
n_ult = 1.5 * n_max                # ultimate loading factor
V_C2 = 0.8 * V_C                   #second cruise speed [m/s]
V_loiter = 0.6 * V_C               #loiter speed [m/s]
V_C_TAS = V_C * 0.514444444 / ((rho_c/rho)**0.5)  # m/s cruise speed, according to flight envelope, 
V_C_estimate = V_C_TAS             # m/s, design parameter
mach_cruise = V_C_estimate/a       #Cruise mach
#------------------------------------------------------------------------------------------------------------------




# Landing gear sizing
#------------------------------------------------------------------------------------------------------------------
theta = 15                         #Clearance angle [deg]
#------------------------------------------------------------------------------------------------------------------

d_wheel_main_lg = 1             #Main landing gear diameter [m]
d_wheel_nose_lg = 0.5

strut_length_main_lg = 1        #[m]
strut_length_nose_lg = 0.75
#Lift/drag ratios
#------------------------------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------------------------------
LD_c = 15                          # L/D cruise, determine after wing sizing and drag analysis
LD_c2 = 17                         # L/D cruise2, determine after wing sizing and drag analysis
LD_loiter = 17                     # L/D Loiter, determine after wing sizing and drag analysis


#Relevant distances for aerodynamic calculations, tail sizing, landing gear sizing
#------------------------------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------------------------------


lh = 15                            # very random estimate, distance between wing and tail QUARTER CHORD
lv = 16                            # very random estimate, distance between wing and vertical tail aerodynamic centers
x_apu = lf-0.8                        # cg location of the apu measured from the nose of the aircraft [m], TBD

x_engine = 13                      # cg location of engines, measured from the nose of the aircraft [m], TBD
x_nacelle = 13                     # cg location of engine nacelles, measured from the nose of the aircraft [m], TBD

x_engine_start = - 1.5             #m, RANDOM, begin of engine measured from the lemac, negative means ... [m] closer to the nose

#Drag parameters
k = 0.152 * 10**-5         # Surface factor for skin friction coefficient, for polished sheet metal (need to reconsider if composites are used)
IF_wing   = 1.0         # Interference factors
IF_tailv  = 1.0
IF_tailh  = 1.04
IF_fus    = 1.0
IF_nacelle = 1.0

