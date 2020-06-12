# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 09:33:37 2020

@author: Rick
"""

# loading diagram generator

import matplotlib.pyplot as plt
import numpy as np
import Class_2_estimation as cl2
import input
from loading_diagram_generator import wing_cg, cg_OEW_wrt_lemac, loadingcg, passenger_loading
from cabindesign import cabin_design
#Several cabin and fuel config parameters

t_cyl,m_cyl, tm_cyl, d_cyl,l_cyl,t_tail,m_tail, tm_tail, d_tail,l_tail\
           ,t_top,m_top,tm_top,d_top,l_top,t_pod,m_pod,tm_pod,d_pod,l_pod,totalcabinlength,V_tank_cyl, V_tank_tail, V_tank_top,V_tank_pod,\
           tm_tanksystem,CGtank,CGfuelfull,CGcomb,totdrag,fuselage_weight,CDzerofus,FFbody,Cfturb,fuselage_area,CDzeropods,fusdrag,poddrag,tailcone_length=cabin_design(0.5,1,26,0)
           
           
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
x_apu = input.x_apu             #cg location of the apu measured from the nose of the aircraft [m]

Cr = input.Cr                   #wing root chord length [m]
Ct = input.Ct                   #wing tip chord length [m]
b = cl2.b                       #wing span [m]
sweep = input.LE_sweep          #Leading edge wing sweep (If we use forward sweep, please let Rick know)
x_cargo_fwd = input.x_cg_fwd_cargo #front cargo cg measured from nose [m]
x_cargo_aft = input.x_cg_aft_cargo #aft cargo cg measured from nose [m]
pax_abreast = input.pax_abreast

x_start_Cr = np.arange(0.2 * l_f, 0.8 * l_f, 0.1)   #x-location measured from the nose where root chord starts
x_lemac = [i + x_lemac_Cr for i in x_start_Cr]
x_engine = x_lemac       #Assume engine cg is at lemac
print('Make x_eninge not wing location dependent if we end up having fuselage mounted engines')
x_nacelle = x_engine     #assume nacelle cg is at engine cg
x_ac = [i + 0.25 * MAC for i in x_lemac]            #assume ac at quarter chord point
lh_fix = input.lh                   #distance between wing and horizontal tail aerodynamic centers
lv_fix = input.lv                   #distance between wing and vertical tail aerodynamic centers
xh = input.x_lemac_rootchord_h + input.x_rootchord_h + 0.25*input.c_mac_h
lh = [xh - i - 0.25*MAC for i in x_start_Cr]
lv = [lv_fix - (i - input.x_start_Cr) for i in x_start_Cr]


#Small calculations with raw inputs
pax_cabin = Npax * w_person
fwd_cargo_max = cargo * input.cargo_fwd_fraction
aft_cargo_max = cargo * input.cargo_aft_fraction


seatloc = []
rows = Npax/pax_abreast
for j in range(int(rows)):
    if j > 7:
        emergency_exit = 10 * 0.0254
    else:
        emergency_exit = 0
    row = seat_start + j * pitch * 0.0254 + emergency_exit  # convert to meters
    seatloc.append(row)





#Calculate x_cg & OEW
w_engine = cl2.df['SRA']['Engines']  # kg
w_nacelle = cl2.df['SRA']['Nacelle']  # kg  
w_empennage = cl2.df['SRA']['Empennage']    #kg
w_wing = cl2.df['SRA']['Wing group'] #kg 
w_apu = cl2.df['SRA']['APU']    #kg


# w_tank = cl2.df['SRA']['Hydrogen tanks']
x_pod_tank = np.array(x_lemac)
x_cyl_tank=totalcabinlength+l_cyl/2+input.cockpit_length
x_tail_tank=totalcabinlength+l_cyl+l_tail/2+input.cockpit_length
w_pod_tank=tm_pod
w_tail_tank=tm_tail
w_cyl_tank=tm_cyl

w_tank = w_pod_tank+w_tail_tank+w_cyl_tank
x_tank = (w_pod_tank * x_pod_tank + w_tail_tank * x_tail_tank + w_cyl_tank * x_cyl_tank)/(w_tank)

w_pod_fuel=V_tank_pod*input.rho_hydrogen
w_tail_fuel=V_tank_tail*input.rho_hydrogen
w_cyl_fuel=V_tank_cyl*input.rho_hydrogen
# x_fuel = x_tank                 #fuel cg measured from nose, assumed same as tank cg as most likely the tank will be symmetrical

x_fuel_fuselage = (w_tail_fuel *x_tail_tank + w_cyl_fuel * x_cyl_tank)/(w_tail_fuel + w_cyl_fuel)
w_fuel_fuselage = w_tail_fuel + w_cyl_fuel

w_lg_main = cl2.df['SRA']['Main LG']    #kg
w_lg_front = cl2.df['SRA']['Nose LG']    #kg


x_empennage = [x_ac[i] + (lh[i] + lv[i]) / 2 for i in range(len(x_start_Cr))] #Assume cg of empennage is in the middle of the aerodynamic center of horizontal and vertical tail, measured from the nose
x_lg_front = input.x_lg_front    #cg location of front landing gear [m], measured from the nose, assumed to be 3 m (used for calculating cg at oew, not to be changed per se)
x_lg_main = [i + 4.7 for i in x_start_Cr]     #cg location of main landing gear [m], assumed 2/3 root chord length further than start of root chord (used for calculating cg at oew, not to be changed per se)
print("In calculation of cg @ OEW, take into account the exact tank placement and cg location once agreed on a specific configuration")

def cg_excursion_wing_shift():
    plt.close()
    cg_fwd_excursion_lst = []
    cg_aft_excursion_lst = []
    cg_loaded_lst = []
    for i in range(len(x_start_Cr)):
        x_cg_wing_nose, x_cg_wing_mac = wing_cg(sweep, b, Cr, Ct, MAC, x_lemac_Cr, x_lemac[i])
        cg_oew_wrt_lemac, cg_oew_nose = cg_OEW_wrt_lemac(x_engine[i], w_engine, x_nacelle[i], w_nacelle, x_empennage[i], w_empennage, x_apu, w_apu, x_tank[i], w_tank, x_cg_wing_nose, w_wing, x_lg_front, w_lg_front, x_lg_main[i], w_lg_main, OEW, x_lemac[i], MAC)
       
        onlyfuel = loadingcg(OEW, cg_oew_nose, 0* w_fuel_fuselage, x_fuel_fuselage)
        onlyfwdcargo = loadingcg(OEW, cg_oew_nose, fwd_cargo_max, x_cargo_fwd)
        onlyaftcargo = loadingcg(OEW, cg_oew_nose, aft_cargo_max, x_cargo_aft)
        bothcargo = loadingcg(onlyfwdcargo[1], onlyfwdcargo[0], aft_cargo_max, x_cargo_aft)
        window = passenger_loading(bothcargo[1], bothcargo[0], multiplication=2)
        window_back = passenger_loading(bothcargo[1], bothcargo[0], multiplication=2, seatloc=seatloc[::-1])
        middle = passenger_loading(window[1][-1], window[0][-1], multiplication=2)
        middle_back = passenger_loading(window[1][-1], window[0][-1], multiplication=2, seatloc=seatloc[::-1])
        aisle = passenger_loading(middle[1][-1], middle[0][-1])
        aisle_back = passenger_loading(middle[1][-1], middle[0][-1], seatloc=seatloc[::-1])
        # fully_loaded = loadingcg(aisle[1][-1], aisle[0][1], fuel_weight, x_fuel)
        
        onlyfuselagefuel = loadingcg(aisle[1][-1], aisle[0][-1], w_fuel_fuselage, x_fuel_fuselage)
        onlypodfuel = loadingcg(aisle[1][-1], aisle[0][-1], w_pod_fuel, x_pod_tank[i])

        bothfuel = loadingcg(onlyfuselagefuel[1], onlyfuselagefuel[0], w_pod_fuel, x_pod_tank[i])
        # plt.plot(100 * (np.array([aisle[0][-1], onlyfuselagefuel[0], bothfuel[0]]) - x_lemac) / MAC,
        #                   [MZF, onlyfuselagefuel[1], bothfuel[1]], marker='^', color='cyan', label = 'Hydrogen')
        
        # bothfuel2 = loadingcg(onlypodfuel[1], onlypodfuel[0], w_fuel_fuselage, x_fuel_fuselage)
        # plt.plot(100 * (np.array([aisle[0][-1], onlypodfuel[0], bothfuel2[0]]) - x_lemac) / MAC,
        #                   [MZF, onlypodfuel[1], bothfuel2[1]], marker='^', color='brown', label = 'Hydrogen fwd first')
        
        
        cg_excursion = np.array([ onlyfuselagefuel[0], onlypodfuel[0], bothfuel[0]]) 

        cg_loaded_lst.append(bothfuel[0])
        cgmin_lst = []
        cgmax_lst = []
        for j in range(len(cg_excursion)):
            cgmin = np.min(cg_excursion[j])
            cgmax = np.max(cg_excursion[j])
            cgmin_lst.append(cgmin)
            cgmax_lst.append(cgmax)
        
        cg_fwd = (np.min(cgmin_lst) - x_lemac[i]) / MAC * 0.97      #subtract 2% margin, assuming most forward cg is after lemac
        cg_aft = (np.max(cgmax_lst) - x_lemac[i]) / MAC * 1.03      #add 2% margin
        cg_fwd_excursion_lst.append(cg_fwd)
        cg_aft_excursion_lst.append(cg_aft)
        
    plt.plot(cg_fwd_excursion_lst, [i / l_f for i in x_lemac], cg_aft_excursion_lst, [i / l_f for i in x_lemac])
    plt.xlabel('x_cg / MAC [-]')
    plt.ylabel('x_lemac / l_fus [-]')
    plt.show()
    

    
    return cg_fwd_excursion_lst, cg_aft_excursion_lst, cg_loaded_lst
cg_fwd_excursion_lst, cg_aft_excursion_lst, cg_loaded_lst = cg_excursion_wing_shift()
print (min(cg_fwd_excursion_lst),max(cg_aft_excursion_lst))





