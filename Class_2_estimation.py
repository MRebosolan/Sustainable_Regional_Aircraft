"""
Created with love by Jorn and Matteo

inputs: 
    -Class 1 weight estimation (OEW, MTOW, fuel weight, payload weight, zero fuel weight)
    - any system weight that is calculated in more detail!
    - Aspect Ratio
    - half chord sweep
    - ultimate load factor
    - thickness over chord
    - taper ratio
    - design mach number at sea level
    SL rho
    Length of fuselage
    height of fuselage
    width of fuselage
    max thickness at wing root
    fuselage gross shell area
    length between 1/4 mac of wing and tail
    gear constant Kgr
    passenger cabin volume
    length of passenger cabin
    number of pax
    number of pilots
    number of cabin crew
    cabin pressure
    freight floor area
    specific fuel consumption
    vertical tail area
    vertical tail span
    horizontal tail area
    half chord sweep of both tails
    height of placement of horizontal tail wrt vertical tail
    range
    number of engines
    number of fuel tanks
    wingloading
    powerloading
    hydrogen to kerosene ratio

OUTPUTS:
    - weights of all subsystems
    - Improved OEW
    - a nice graph of class 1 and 2 converging
    - a table in latex with subsystem weights
    - a new take-off thrust
    - a new wing size
    
PURPOSE: 
    - make a class 2 weight estimation, and let it converge with the class 1 weight estimation, while providing updated aircraft parameters
"""





import APU_weight_estimation as apu
import Cargo_handling_weight_estimation as cargo
import electrical_system_weight_estimation as electrical
import Engine_weight_estimation as engine
import flight_controls_weight_estimation as flightcontrols
import fuel_system_weight_estimation as fuelsystem
import furnishing_weight_estimation as furnishing
import fuselage_weight_estimation as fuselage
import instrumentation_weight_estimation as instrumentation
import LG_weight_estimation as LG
import nacelle_weight_estimation as nacelles
import Oxygen_system_weight_estimation as oxygen
import paint_weight_estimation as paint
import power_controls_weight_estimation as powercontrols
import Pressurization_system_weight_estimation as pressurization
import Tail_weight_estimation as tail
import Wing_weight_estimation as wing
from Class_1_estimation import CLASS1WEIGHTHYBRID
import input


import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def to_pounds(kg):
    return kg * 2.20462262

def to_kg (lbs):
    return lbs/2.20462262

def toft (m):
    return m /0.3048

def tom (ft):
    return ft*0.3048

def tosqft(m2):
    return m2/(0.3048*0.3048)

def tom2(sqft):
    return sqft *0.3048*0.3048

#------------- INPUT PARAMETERS ----------------#

AR = input.AR
half_sweep = input.half_sweep
n_max =  input.n_max
n_ult = 1.5* n_max

# S = tosqft(input.S)
t_over_c = input.t_over_c
taper = input.taper	
mach_h = input.mach_h
rho = input.rho *0.0624279606
rho_zero = input.rho /515.378818
V_dive = input.V_dive #is already in KEAS
lf = toft(input.lf)
hf = toft(input.hf)
A_inlet = tosqft(input.A_inlet)
ln = toft(input.ln)


# b = toft(input.b)
t_r= toft(input.t_r) #max thickness at root
widthf = toft(input.widthf) #max fuselage width
S_fgs = tosqft(input.S_fgs) #fuselage gross shell area
lh = toft(input.lh) #ft
# T_TO = to_pounds(input.T_TO/9.81)

Kgr = input.Kgr # constant
V_pax = input.V_pax / (0.3048**3) #ft3
lpax = toft(input.lpax)  #ft
Npax =input.Npax #
N_fdc =input.N_fdc #
N_cc =input.N_cc #
P_c = input.P_c * 0.02089 #should be in psf
Sff = tosqft(input.Sff)

K_fsp = input.K_fsp * 8.34540445
#tails:
Sv = tosqft(input.Sv)
half_chord_sweep_hor = input.half_chord_sweep_hor
half_chord_sweep_vert = input.half_chord_sweep_vert
bv = toft(input.bv)
Sh = tosqft(input.Sh)
zh = toft(input.zh)

range = input.Design_range *0.539956803
N_eng = input.N_eng
N_t = input.N_t
wingloading = input.wingloading #still metric
powerloading = input.powerloading

W_payload = to_pounds(input.W_payload)


ratio = input.H_to_ker_ratio

#setting initial values to zero
OEWINPUT = 1
OEW_class1_kg = 2
OEW_plot_class1 = []
OEW_plot_class2 = []



while abs((OEW_class1_kg - OEWINPUT)*100/OEWINPUT)>= 0.01:
    class1 = CLASS1WEIGHTHYBRID(ratio,OEWINPUT)
    MTOW_kg = class1[0]
    S_metric = MTOW_kg*9.81 /wingloading
    b_metric = (S_metric * AR)**0.5
    T_TO_newton = MTOW_kg *9.81 *powerloading
    T_TO = to_pounds(T_TO_newton/9.81)
    T_dry_SL = 0.5* T_TO
    b = toft(b_metric)
    S = tosqft(S_metric)
    
    
    OEW_class1_kg = class1[1]
    OEW_plot_class1.append(OEW_class1_kg)
    M_zfw_kg = class1[4]
    M_fuel_kg = class1[5]+class1[6]
    W_fuel = to_pounds(M_fuel_kg)
    
    OEW_class1 = to_pounds(OEW_class1_kg)
    MTOW = to_pounds(MTOW_kg)
    W_zfw = to_pounds(M_zfw_kg)


#--------- STRUCTURAL WEIGHT --------------#

    W_wing_gd = wing.gd_wing(MTOW, AR, half_sweep, n_ult, S, t_over_c, taper, mach_h)
    W_fuselage_GD = fuselage.W_fuselage_gd (rho_zero, V_dive, MTOW, lf, hf)

    
    W_wing = wing.W_wing(W_zfw, b, half_sweep, n_ult, S, t_r)
    W_empennage = tail.vert_tail_weight(Sv, V_dive, half_chord_sweep_vert, bv, Sh, zh) + tail.hor_tail_weight(Sh, V_dive, half_chord_sweep_hor)
    W_fuselage = fuselage.W_fuselage_torenbeek(V_dive, lh, widthf, hf, S_fgs)
    W_nacelles = nacelles.W_nacelle_torenbeek(T_TO)
    W_landing_gear_nose = LG.LG_weight(Kgr, MTOW)[1]
    W_landing_gear_main= LG.LG_weight(Kgr, MTOW)[0]
    W_landing_gear = W_landing_gear_nose + W_landing_gear_main
    
    #--------- EQUIPMENT WEIGHT ------------#
    

    instrumentation_weight = (instrumentation.instrumentation_torenbeek(MTOW)+ instrumentation.instrumentation_2(OEW_class1, range)+ instrumentation.instrumentation_gd(MTOW,N_fdc))/3
    flight_control_weight = flightcontrols.flight_controls(MTOW)                            #Verified for LE devices + spoilers
    electrical_system_weight = electrical.electrical_torenbeek(V_pax)                       #Verified
    airconditioning_pressurization_weight = pressurization.pressure_system_weight(lpax)     #Verified
    oxygen_system_weight = oxygen.oxygen_system_weight(Npax)                                #Verified
    APU_weight = apu.APU_weight_estimation(MTOW)                                            #Verified
    furnishing_weight = furnishing.furnishing_gd(N_fdc, Npax, N_cc, MTOW, P_c)              #Verified
    cargo_equipment_weight = cargo.cargo_handling_weight(Sff)                               #Verified
    paint_weight = paint.paint(MTOW)                                                        #Verified
    
    #------------ POWER PLANT WEIGHT ------------#
    
    W_engines = engine.engine_weight(T_dry_SL, N_eng)                                       #Verified
    W_fuel_system_kerosene = 0
    
    W_fuel_system_hydrogen = to_pounds(class1[9])
    W_power_controls = powercontrols.total(lf, b, W_engines, pneumatic = False)
    
    
    
    
    #------------- SUMMATION OF COMPONENTS -------------#
    
    W_struct = W_wing + W_empennage + W_fuselage + W_nacelles + W_landing_gear
    
    W_powerplant = W_engines + W_fuel_system_kerosene + W_fuel_system_hydrogen + W_power_controls
    
    W_equipment = APU_weight + cargo_equipment_weight + furnishing_weight + instrumentation_weight + oxygen_system_weight + paint_weight + airconditioning_pressurization_weight + flight_control_weight + electrical_system_weight
    
    OEW_class2 = W_struct + W_powerplant + W_equipment
    

    OEWINPUT = to_kg(OEW_class2)
    OEW_plot_class2.append(OEWINPUT)

plt.close()
plt.figure()
plt.plot(np.arange(0, len(OEW_plot_class1)), OEW_plot_class1, label = 'class 1')
plt.plot(np.arange(0, len(OEW_plot_class2)), OEW_plot_class2, label = 'class 2')
plt.xlabel("iterations")
plt.ylabel("OEW in kg")
plt.legend()
plt.show()

df = pd.DataFrame({'Component': ['MTOW','OEW'],
'SRA': [MTOW, OEW_class2],
'F28': [65000, 31219],
'737-200':[115500,60210]})

zfw_fuel=[{'Component': 'Zero fuel weight', 'SRA': W_zfw, 'F28':31219+14380, '737-200': 60210+32790},
         {'Component': 'Max fuel weight', 'SRA': W_fuel, 'F28':17331, '737-200': 34781},
         {'Component': 'Max payload', 'SRA': W_payload, 'F28':14380, '737-200': 34790}]
df = df.append(zfw_fuel, ignore_index = True, sort = False)

wng = [{'Component': 'Wing group', 'SRA': W_wing, 'F28':7330, '737-200': 10613},
       {'Component': 'Empennage', 'SRA': W_empennage, 'F28':1632 , '737-200':2718},
       {'Component': 'Fuselage', 'SRA': W_fuselage, 'F28':7043, '737-200':12108},
       {'Component': 'Nacelle', 'SRA': W_nacelles, 'F28':834, '737-200':1392},
       {'Component': 'Landing gear', 'SRA': W_landing_gear, 'F28':2759, '737-200':4354},
       {'Component': 'Main LG', 'SRA': W_landing_gear_main},
       {'Component': 'Nose LG', 'SRA': W_landing_gear_nose},
       {'Component': 'Total structural', 'SRA': W_struct, 'F28':19598, '737-200':31185}]
df = df.append(wng, ignore_index = True, sort = False)


power =[{'Component': 'Engines', 'SRA': W_engines, 'F28':4495, '737-200':6217},
       {'Component': 'Exhaust', 'F28':127, '737-200':1007},
       {'Component': 'Kerosene system', 'SRA': W_fuel_system_kerosene, 'F28':545, '737-200':575},
       {'Component': 'Hydrogen tanks', 'SRA': W_fuel_system_hydrogen},
       {'Component': 'Power controls', 'SRA': W_power_controls, 'F28':215, '737-200':378},
       {'Component': 'Total propulsion', 'SRA': W_powerplant, 'F28':5382, '737-200':8177}]
df = df.append(power, ignore_index = True, sort = False)

equipment = [{'Component': 'Electrical systems', 'SRA': electrical_system_weight, 'F28':1892, '737-200':1066+956},
             {'Component': 'Instruments', 'SRA': instrumentation_weight, 'F28':302, '737-200':625},
             {'Component': 'Flight controls', 'SRA': flight_control_weight, 'F28':1387+364, '737-200':2348+873},
             {'Component': 'APU', 'SRA': APU_weight, 'F28':346, '737-200':836},
             {'Component': 'Air conditioning', 'SRA': airconditioning_pressurization_weight + oxygen_system_weight, 'F28':1074, '737-200':1416},
             {'Component': 'Furnishing', 'SRA': furnishing_weight, 'F28':4030, '737-200':6643},
             {'Component': 'Cargo handling', 'SRA': cargo_equipment_weight},
             {'Component': 'Miscellanous/paint', 'SRA': paint_weight, '737-200':124},
             {'Component': 'Total fixed equipment', 'SRA': W_equipment, 'F28':9395, '737-200':14887}]
df = df.append(equipment, ignore_index = True, sort = False)

df['SRA fraction'] = (df['SRA']/MTOW).round(3)
df['F28 fraction'] = (df['F28']/df['F28'][0]).round(3)
df['737 fraction'] = (df['737-200']/df['737-200'][0]).round(3)
df['SRA'] = to_kg(df['SRA']).round(1)
df['F28'] = to_kg(df['F28']).round(1)
df['737-200'] = to_kg(df['737-200']).round(1)
df = df.set_index('Component')

x = df.loc['MTOW','SRA']


aircraftpar = pd.DataFrame()
wing = [{'Component': 'Wing Area', 'SRA': S, 'F28':tom2(1), '737-200':tom2(1)},
               ]
aircraftpar = aircraftpar.append(wing, ignore_index = True, sort = False)

df['SRA']['MTOW'] = input.MTOM
# latex = df.to_latex(index = True, caption = "Summary table for Class-II estimation", label = 'tab:class2table', na_rep = 'None')
# print("Uncomment the caption for the final version")

# print(df)
# # print(aircraftpar)



# try:
#     file = open('C://Users//jornv//Google Drive//DSE upload//Class2dataframe.txt', 'w')
#     file.write(latex)
#     file.close()
# except:
#     print()
#     #print('you cannot update files, ask jorn if necessary')
    
    

try:
    file = open('C://Users//jornv//Google Drive//DSE upload//Class2dataframe.txt', 'w')
    file.write(latex)
    file.close()
except:
    print()
    #print('you cannot update files, ask jorn if necessary')

print ('yahoo',to_kg(LG.LG_weight(Kgr, MTOW)[0]))
print ('the nose',to_kg(LG.LG_weight(Kgr, MTOW)[1]))

S = input.S
b = input.b
MTOM = to_kg(MTOW)
MTOW = 9.81*MTOM
OEM = OEWINPUT
Tto = T_TO_newton
print (W_empennage*0.45359237)
print (W_fuselage*0.45359237)
# print (S,b**2/8)
