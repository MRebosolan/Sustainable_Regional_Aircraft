# loading diagram generator

import matplotlib.pyplot as plt
import numpy as np
import Class_2_estimation as cl2
import input
from cabindesign import cabin_design()

#Several cabin and fuel config parameters

t_cyl,m_cyl, tm_cyl, d_cyl,l_cyl,t_tail,m_tail, tm_tail, d_tail,l_tail\
           ,t_top,m_top,tm_top,d_top,l_top,t_pod,m_pod,tm_pod,d_pod,l_pod,totalcabinlength,V_tank_cyl, V_tank_tail, V_tank_top,V_tank_pod,\
           tm_tanksystem,CGtank,CGfuelfull,CGcomb,totdrag,fuselage_weight,CDzerofus,FFbody,Cfturb,fuselage_area,CDzeropods,fusdrag,poddrag,empennage_length=cabin_design(0,0.25,25,0)

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
x_lemac = input.x_LEMAC_nose    #x lemac measured from the nose of the aircraft
x_start_Cr = input.x_start_Cr   #x-location measured from the nose where root chord starts
seat_start = input.x_first_pax  #x-location measured from the nose where first passenger row is located
pitch = input.seat_pitch        #seat pitch [inch]
rows = input.n_rows             #number of passegner rows [-]
lh = input.lh                   #distance between wing and horizontal tail aerodynamic centers
lv = input.lv                   #distance between wing and vertical tail aerodynamic centers
x_ac = x_lemac + MAC / 4               #x location of wing aerodynamic center measured from the nose of the aircraft
x_apu = input.x_apu             #cg location of the apu measured from the nose of the aircraft [m]
x_engine = input.x_engine       #cg location of engines, measured from the nose of the aircraft [m]
x_nacelle = input.x_nacelle     #cg location of engine nacelles, measured from the nose of the aircraft [m]
Cr = input.Cr                   #wing root chord length [m]
Ct = input.Ct                   #wing tip chord length [m]
b = cl2.b                       #wing span [m]
sweep = input.LE_sweep          #Leading edge wing sweep (If we use forward sweep, please let Rick know)
x_cargo_fwd = input.x_cg_fwd_cargo #front cargo cg measured from nose [m]
x_cargo_aft = input.x_cg_aft_cargo #aft cargo cg measured from nose [m]
pax_abreast = input.pax_abreast

#Small calculations with raw inputs
pax_cabin = Npax * w_person
fwd_cargo_max = cargo * input.cargo_fwd_fraction
aft_cargo_max = cargo * input.cargo_aft_fraction


seatloc = []
rows = Npax/pax_abreast
for j in range(int(rows)):
    row = seat_start + j * pitch * 0.0254  # convert to meters
    seatloc.append(row)


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
w_tank = cl2.df['SRA']['Hydrogen tanks']
x_pod_tank = 
x_cyl_tank=totalcabinlength+l_cyl/2+input.cockpit_length
x_tail_tank=totalcabinlength+l_cyl+l_tail/2+input.cockpit_length

print("change w_tank and x_tank to variables used in other files once decided on a fuel tank configuration")
x_fuel = x_tank                 #fuel cg measured from nose, assumed same as tank cg as most likely the tank will be symmetrical
w_lg_main = cl2.df['SRA']['Main LG']    #kg
w_lg_front = cl2.df['SRA']['Nose LG']    #kg
w_fuselage = cl2.df['SRA']['Fuselage']
x_fuselage = l_f/2 
w_powercontrols = cl2.df['SRA']['Power controls']
x_powercontrols = x_engine
w_electrical = cl2.df['SRA']['Electrical systems']
x_electrical = input.x_LEMAC_nose
w_instruments = cl2.df['SRA']['Instruments']
x_instruments = input.x_first_pax /3
w_flightcontrols = cl2.df['SRA']['Flight controls']
x_flightcontrols = x_start_Cr+Cr
w_airconditioning = cl2.df['SRA']['Air conditioning']
x_airconditioning = l_f/2
w_furnishing = cl2.df['SRA']['Furnishing']
x_furnishing = l_f/2
w_cargohandling = cl2.df['SRA']['Cargo handling']
x_cargohandling = x_cargo_aft *input.cargo_aft_fraction + x_cargo_fwd *input.cargo_fwd_fraction 
w_paint = cl2.df['SRA']['Miscellanous/paint']
x_paint = l_f/2


x_empennage = x_ac + (lh + lv) / 2 #Assume cg of empennage is in the middle of the aerodynamic center of horizontal and vertical tail, measured from the nose
x_lg_front = 3     #cg location of front landing gear [m], measured from the nose, assumed to be 3 m (used for calculating cg at oew, not to be changed per se)
x_lg_main = x_start_Cr + 2 * Cr / 3      #cg location of main landing gear [m], assumed 2/3 root chord length further than start of root chord (used for calculating cg at oew, not to be changed per se)
print("In calculation of cg @ OEW, take into account the exact tank placement and cg location once agreed on a specific configuration")


def wing_cg(sweep, b, Cr, Ct, MAC, x_lemac_Cr, x_lemac):
    """
    This function computes the wing cg w.r.t. the nose of the aircraft and the lemac
    """
    c = b / (2 * np.cos(sweep))
    a = np.sqrt(c**2 - (b / 2)**2)
    x_1 = 2 * a / 3
    A1 = 0.5 * b / 2 * a
    A2 = b / 2 * (Cr - a)
    x_2 = a + (Cr - a) / 2
    A3 = 0.5 * b / 2 * (Ct - (Cr - a))
    x_3 = a +  Ct / 3
    Atot = A1 + A2 + A3
    x_cg_wrt_xlemac_Cr = (A1 * (x_1 - x_lemac_Cr) + A2 * (x_2 - x_lemac_Cr) + A3 * (x_3 - x_lemac_Cr)) / Atot
    x_cg_nose = x_start_Cr + + x_lemac_Cr + x_cg_wrt_xlemac_Cr
    x_cg_mac= x_cg_nose - x_lemac 
    return x_cg_nose, x_cg_mac

x_cg_wing_nose, x_cg_wing_mac = wing_cg(sweep, b, Cr, Ct, MAC, x_lemac_Cr, x_lemac)
print(x_cg_wing_nose, x_cg_wing_mac)
#vary x_start_Cr

print((w_engine +  w_nacelle +  w_empennage +  w_apu + w_tank +  w_wing +  w_lg_front +  w_lg_main + w_fuselage  + w_powercontrols \
                   + w_electrical  + w_instruments + w_flightcontrols  + w_airconditioning + w_furnishing + w_cargohandling +w_paint) / OEW, 'OEW')


def cg_OEW_wrt_lemac(x_engine, w_engine, x_nacelle, w_nacelle, x_empennage, w_empennage, x_apu, w_apu, x_tank, w_tank, x_cg_wing_nose, w_wing, x_lg_front, w_lg_front, x_lg_main, w_lg_main, OEW, x_lemac, MAC):
    cg_oew_nose = (x_engine * w_engine + x_nacelle * w_nacelle + x_empennage * w_empennage + x_apu * w_apu + x_tank * w_tank + x_cg_wing_nose * w_wing + x_lg_front * w_lg_front + x_lg_main * w_lg_main + w_fuselage * x_fuselage + w_powercontrols* x_powercontrols \
                   + w_electrical * x_electrical + w_instruments*x_instruments + w_flightcontrols * x_flightcontrols + w_airconditioning*x_airconditioning + w_furnishing*x_furnishing + w_cargohandling*x_cargohandling + w_paint*x_paint) / OEW
    cg_oew_wrt_lemac = (cg_oew_nose - x_lemac) / MAC
    return cg_oew_wrt_lemac, cg_oew_nose




cg_oew_wrt_lemac, cg_oew_nose =  cg_OEW_wrt_lemac(x_engine, w_engine, x_nacelle, w_nacelle, x_empennage, w_empennage, x_apu, w_apu, x_tank, w_tank, x_cg_wing_nose, w_wing, x_lg_front, w_lg_front, x_lg_main, w_lg_main, OEW, x_lemac, MAC)  
print('C.G. @ OEW = ', cg_oew_wrt_lemac, 'MAC')   
    

# new cg calculation, returns new cg and new weight
def loadingcg(w_old, cg_old, w_item, cg_item):
    x_cg_new = (w_old * cg_old + w_item * cg_item) / (w_old + w_item)
    return x_cg_new, w_old + w_item



def maccie(x, x_lemac, MAC):
    return 100 * (x - x_lemac) / MAC



############# PLOTTING BELOW -----------------------------


# passenger loading calculation with new cg, outputs lists with cgs and weights
def passenger_loading(current_weight, current_cg, multiplication=1, seatloc=seatloc):
    cglist = [current_cg]
    weight_list = [current_weight]

    for loc in seatloc:
        cg_and_weight = loadingcg(current_weight, current_cg, w_person * multiplication, loc)
        current_cg = cg_and_weight[0]
        cglist.append(current_cg)
        current_weight = cg_and_weight[1]
        weight_list.append(current_weight)

    return cglist, weight_list


def loading():

    plt.close()
    plt.figure()
    onlyfwdcargo = loadingcg(OEW, cg_oew_nose, fwd_cargo_max, x_cargo_fwd)
    bothcargo = loadingcg(onlyfwdcargo[1], onlyfwdcargo[0], aft_cargo_max, x_cargo_aft)
    cargo1 = plt.plot(100 * (np.array([cg_oew_nose, onlyfwdcargo[0], bothcargo[0]]) - x_lemac) / MAC,
                      [OEW, onlyfwdcargo[1], bothcargo[1]], label='Cargo', marker='x', color='brown')

    onlyaftcargo = loadingcg(OEW, cg_oew_nose, aft_cargo_max, x_cargo_aft)
    bothcargo2 = loadingcg(onlyaftcargo[1], onlyaftcargo[0], fwd_cargo_max, x_cargo_fwd)
    cargo2 = plt.plot(100 * (np.array([cg_oew_nose, onlyaftcargo[0], bothcargo2[0]]) - x_lemac) / MAC,
                      [OEW, onlyaftcargo[1], bothcargo2[1]], marker='x', color='cyan')

    window = passenger_loading(bothcargo[1], bothcargo[0], multiplication=2)
    window_back = passenger_loading(bothcargo[1], bothcargo[0], multiplication=2, seatloc=seatloc[::-1])
    window1 = plt.plot(100 * (np.array(window[0]) - x_lemac) / MAC, window[1], label='Window seats', marker='.',
                       color='blue')
    plt.plot(100 * (np.array(window_back[0]) - x_lemac) / MAC, window_back[1], marker='.', color='orange')

    middle = passenger_loading(window[1][-1], window[0][-1], multiplication=2)
    middle_back = passenger_loading(window[1][-1], window[0][-1], multiplication=2, seatloc=seatloc[::-1])
    plt.plot(100 * (np.array(middle[0]) - x_lemac) / MAC, middle[1], label='Middle seats', marker='1', color='black')
    plt.plot(100 * (np.array(middle_back[0]) - x_lemac) / MAC, middle_back[1], marker='1', color='yellow')

    aisle = passenger_loading(middle[1][-1], middle[0][-1])
    aisle_back = passenger_loading(middle[1][-1], middle[0][-1], seatloc=seatloc[::-1])
    plt.plot(100 * (np.array(aisle[0]) - x_lemac) / MAC, aisle[1], label='Aisle seats', marker='2', color='red')
    plt.plot(100 * (np.array(aisle_back[0]) - x_lemac) / MAC, aisle_back[1], marker='2', color='green')

    fully_loaded = loadingcg(aisle[1][-1], aisle[0][1], fuel_weight, x_fuel)
    plt.plot(100 * (np.array([aisle[0][-1], fully_loaded[0]]) - x_lemac) / MAC, [MZF, fully_loaded[1]], marker='^',
             color='magenta', label='Fuel')
    
    plt.legend()
    plt.grid()
    plt.ylabel('mass [kg]')
    plt.xlabel('xcg [% of MAC]')
    plt.show()
    cg_excursion = np.array([[onlyfwdcargo[0]], [onlyaftcargo[0]], [bothcargo[0]], [window[0]], window_back[0], 
                             middle[0], middle_back[0], aisle[0], aisle_back[0], fully_loaded[0]]) 
    cgmin_lst = []
    cgmax_lst = []
    for i in range(len(cg_excursion)):
        cgmin = np.min(cg_excursion[i])
        cgmax = np.max(cg_excursion[i])
        cgmin_lst.append(cgmin)
        cgmax_lst.append(cgmax)
     
        
    cg_fwd = (np.min(cgmin_lst) - x_lemac) / MAC * 0.98      #subtract 2% margin, assuming most forward cg is after lemac
    cg_aft = (np.max(cgmax_lst) - x_lemac) / MAC *1.02      #add 2% margin
    print("most forward cg should be positive")
   
    
    return cg_fwd, cg_aft


cg_fwd, cg_aft = loading()  
print('CG range:', cg_fwd, '-', cg_aft, 'MAC')
