# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 09:33:37 2020

@author: Gebruiker
"""

# loading diagram generator

import matplotlib.pyplot as plt
import numpy as np
import Class_2_estimation as cl2
import input
from loading_diagram_generator import wing_cg, cg_OEW_wrt_lemac, loadingcg, passenger_loading

#Raw inputs
MTOW = cl2.MTOM                 #kg
OEW = cl2.OEM                   #kg
MZF = cl2.M_zfw_kg              #kg
Npax = input.Npax               #Number of passengers [-]
w_person = input.W_pax          #Weight of each passenger + luggage [kg]
l_f = input.lf                  #Fuselage length [m]
cargo = input.W_cargo           #Total cargo weight [kg]
fuel_weight = cl2.M_fuel_kg     #Total fuel weight taken onboard, including reserves [kg]
MAC = input.MAC                 #Length of MAC [m]
y_MAC = input.y_MAC             #spanwise location of mean aerodynamic chord, measured from the centerline of the aircraft[m]
x_lemac_Cr = input.x_lemac_rootchord #x location of leading edge mac measured from root chord [m]
seat_start = input.x_first_pax  #x-location measured from the nose where first passenger row is located
pitch = input.seat_pitch        #seat pitch [inch]
rows = input.n_rows             #number of passegner rows [-]
lh = input.lh                   #distance between wing and horizontal tail aerodynamic centers
lv = input.lv                   #distance between wing and vertical tail aerodynamic centers
x_apu = input.x_apu             #cg location of the apu measured from the nose of the aircraft [m]

Cr = input.Cr                   #wing root chord length [m]
Ct = input.Ct                   #wing tip chord length [m]
b = cl2.b                       #wing span [m]
sweep = input.LE_sweep          #Leading edge wing sweep (If we use forward sweep, please let Rick know)
x_cargo_fwd = input.x_cg_fwd_cargo #front cargo cg measured from nose [m]
x_cargo_aft = input.x_cg_aft_cargo #aft cargo cg measured from nose [m]
pax_abreast = input.pax_abreast

x_start_Cr = np.arange(0.2 * l_f, 0.8 * l_f, 1)   #x-location measured from the nose where root chord starts
x_lemac = [i + x_lemac_Cr for i in x_start_Cr]
x_engine = x_lemac       #Assume engine cg is at lemac
x_nacelle = x_engine     #assume nacelle cg is at engine cg
x_ac = [i + 0.25 * MAC for i in x_lemac]            #assume ac at quarter chord point
#Small calculations with raw inputs
pax_cabin = Npax * w_person
fwd_cargo_max = cargo * input.cargo_fwd_fraction
aft_cargo_max = cargo * input.cargo_aft_fraction

seatloc = []
rows = Npax/pax_abreast
for j in range(int(rows)):
    row = seat_start + j * pitch * 0.0254  # convert to meters
    seatloc.append(row)




#Calculate x_cg & OEW
w_engine = cl2.df['SRA']['Engines']  # kg
w_nacelle = cl2.df['SRA']['Nacelle']  # kg  
w_empennage = cl2.df['SRA']['Empennage']    #kg
w_wing = cl2.df['SRA']['Wing group'] #kg 
w_apu = cl2.df['SRA']['APU']    #kg
w_tank = 500.
x_tank = 20.
print("change w_tank and x_tank to variables used in other files once decided on a fuel tank configuration")
x_fuel = x_tank                 #fuel cg measured from nose, assumed same as tank cg as most likely the tank will be symmetrical
w_lg_main = cl2.df['SRA']['Main LG']    #kg
w_lg_front = cl2.df['SRA']['Nose LG']    #kg


x_empennage = [i + (lh + lv) / 2 for i in x_ac] #Assume cg of empennage is in the middle of the aerodynamic center of horizontal and vertical tail, measured from the nose
x_lg_front = 3     #cg location of front landing gear [m], measured from the nose, assumed to be 3 m (used for calculating cg at oew, not to be changed per se)
x_lg_main = [i + 2 * Cr / 3 for i in x_start_Cr]     #cg location of main landing gear [m], assumed 2/3 root chord length further than start of root chord (used for calculating cg at oew, not to be changed per se)
print("In calculation of cg @ OEW, take into account the exact tank placement and cg location once agreed on a specific configuration")

def cg_excursion_wing_shift():
    plt.close()
    cg_fwd_excursion_lst = []
    cg_aft_excursion_lst = []
    for i in range(len(x_start_Cr)):
        x_cg_wing_nose, x_cg_wing_mac = wing_cg(sweep, b, Cr, Ct, MAC, x_lemac_Cr, x_lemac[i])
        cg_oew_wrt_lemac, cg_oew_nose = cg_OEW_wrt_lemac(x_engine[i], w_engine, x_nacelle[i], w_nacelle, x_empennage[i], w_empennage, x_apu, w_apu, x_tank, w_tank, x_cg_wing_nose, w_wing, x_lg_front, w_lg_front, x_lg_main[i], w_lg_main, OEW, x_lemac[i], MAC)
        onlyfwdcargo = loadingcg(OEW, cg_oew_nose, fwd_cargo_max, x_cargo_fwd)
        onlyaftcargo = loadingcg(OEW, cg_oew_nose, aft_cargo_max, x_cargo_aft)
        bothcargo = loadingcg(onlyfwdcargo[1], onlyfwdcargo[0], aft_cargo_max, x_cargo_aft)
        window = passenger_loading(bothcargo[1], bothcargo[0], multiplication=2)
        window_back = passenger_loading(bothcargo[1], bothcargo[0], multiplication=2, seatloc=seatloc[::-1])
        middle = passenger_loading(window[1][-1], window[0][-1], multiplication=2)
        middle_back = passenger_loading(window[1][-1], window[0][-1], multiplication=2, seatloc=seatloc[::-1])
        aisle = passenger_loading(middle[1][-1], middle[0][-1])
        aisle_back = passenger_loading(middle[1][-1], middle[0][-1], seatloc=seatloc[::-1])
        fully_loaded = loadingcg(aisle[1][-1], aisle[0][1], fuel_weight, x_fuel)
        cg_excursion = np.array([[onlyfwdcargo[0]], [onlyaftcargo[0]], [bothcargo[0]], [window[0]], window_back[0], 
                             middle[0], middle_back[0], aisle[0], aisle_back[0], fully_loaded[0]]) 
        cgmin_lst = []
        cgmax_lst = []
        for j in range(len(cg_excursion)):
            cgmin = np.min(cg_excursion[j])
            cgmax = np.max(cg_excursion[j])
            cgmin_lst.append(cgmin)
            cgmax_lst.append(cgmax)
        
        cg_fwd = (np.min(cgmin_lst) - x_lemac[i]) / MAC * 0.98      #subtract 2% margin, assuming most forward cg is after lemac
        cg_aft = (np.max(cgmax_lst) - x_lemac[i]) / MAC * 1.02      #add 2% margin
        cg_fwd_excursion_lst.append(cg_fwd)
        cg_aft_excursion_lst.append(cg_aft)
        
    plt.plot(cg_fwd_excursion_lst, [i / l_f for i in x_lemac], cg_aft_excursion_lst, [i / l_f for i in x_lemac])
    plt.xlabel('x_cg / MAC [-]')
    plt.ylabel('x_lemac / l_fus [-]')
    plt.show()
    
    return 
cg_excursion_wing_shift()




