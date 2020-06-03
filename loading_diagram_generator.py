# loading diagram generator

import matplotlib.pyplot as plt
import numpy as np
import Class_2_estimation as cl2
import input


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
x_ac = input.x_ac               #x location of wing aerodynamic center measured from the nose of the aircraft
x_apu = input.x_apu             #cg location of the apu measured from the nose of the aircraft [m]
x_engine = input.x_engine       #cg location of engines, measured from the nose of the aircraft [m]
x_nacelle = input.x_nacelle     #cg location of engine nacelles, measured from the nose of the aircraft [m]
Cr = input.Cr                   #wing root chord length [m]
Ct = input.Ct                   #wing tip chord length [m]
b = cl2.b                       #wing span [m]
sweep = input.LE_sweep          #Leading edge wing sweep (If we use forward sweep, please let Rick know)
x_cargo_fwd = input.x_cg_fwd_cargo #front cargo cg measured from nose [m]
x_cargo_aft = input.x_cg_aft_cargo #aft cargo cg measured from nose [m]


#Small calculations with raw inputs
pax_cabin = Npax * w_person
fwd_cargo_max = cargo * input.cargo_fwd_fraction
aft_cargo_max = cargo * input.cargo_aft_fraction

seatloc = []
for j in range(14):
    row = seat_start + j * pitch * 0.0254  # convert to meters
    seatloc.append(row)




#Calculate x_cg & OEW
w_engine = cl2.df['SRA']['Engines']  # kg
w_nacelle = cl2.df['SRA']['Nacelle']  # kg  
w_empennage = cl2.df['SRA']['Empennage']    #kg
w_wing = cl2.df['SRA']['Wing group'] #kg 
w_apu = cl2.df['SRA']['APU']    #kg
w_tank = 500
x_tank = 20
print("change w_tank and x_tank to variables used in other files once decided on a fuel tank configuration")
x_fuel = x_tank                 #fuel cg measured from nose, assumed same as tank cg as most likely the tank will be symmetrical
w_lg_main = cl2.df['SRA']['Main LG']    #kg
w_lg_front = cl2.df['SRA']['Nose LG']    #kg

x_engine = 15       #x_location of c.g. of engines measured from the nose [m]
x_nacelle = 15      #x_location of c.g. of engine nacelles measured from the nose [m]
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
    x_cg_nose = x_start_Cr + x_cg_wrt_xlemac_Cr
    x_cg_mac= x_cg_nose - x_lemac 
    return x_cg_nose, x_cg_mac

x_cg_wing_nose, x_cg_wing_mac = wing_cg(sweep, b, Cr, Ct, MAC, x_lemac_Cr, x_lemac)

#vary x_start_Cr



def cg_OEW_wrt_lemac(x_engine, w_engine, x_nacelle, w_nacelle, x_empennage, w_empennage, x_apu, w_apu, x_tank, w_tank, x_cg_wing_nose, w_wing, x_lg_front, w_lg_front, x_lg_main, w_lg_main, OEW, x_lemac, MAC):
    cg_oew_nose = (x_engine * w_engine + x_nacelle * w_nacelle + x_empennage * w_empennage + x_apu * w_apu + x_tank * w_tank + x_cg_wing_nose * w_wing + x_lg_front * w_lg_front + x_lg_main * w_lg_main) / OEW
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


# plot potato diagram
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
    window2 = plt.plot(100 * (np.array(window_back[0]) - x_lemac) / MAC, window_back[1], marker='.', color='orange')

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
    extreme_cg = maccie(onlyfwdcargo[0], x_lemac, MAC), maccie(onlyaftcargo[0], x_lemac, MAC)
    return fully_loaded[0], extreme_cg


cgMTOW, extreme_cg = loading()  # extreme cg is in % of mac, first fwd, then aft, fully loaded is cg















# tail_armh = xlemac_horizontal + 0.25 * MAC_horizontal - cgMTOW
# tail_armv = xlemac_vertical + 0.25 * MAC_vertical - cgMTOW
# C_h = horizontal_area * tail_armh / (wing_area * MAC)
# C_v = vertical_area * tail_armv / (wing_area * span)

# clalpha = (1 - -1) / 17.5
# cmalpha = (0.0329 + 0.0093) / (7.25 + 13.75)
# cl0 = -1.25  # degrees
# cm0 = -0.0184
# area_fuselage = diameter * 14 * (10 / 21)
# mach = 0.78
# v_approach = 129  # kts
# beta = (1 - 0.85 * (mach ** 2)) ** 0.5
# n = 0.95
# stabilizer_sweep = np.radians(27)
# span_tail = 12.3
# AR_tail = span_tail ** 2 / horizontal_area
# AR = span ** 2 / wing_area
# sweep = np.radians(24.5)

# clalpha_datcom = 2 * np.pi * AR / (2 + ((4 + ((AR * beta) ** 2) * (1 + (np.tan(sweep) ** 2) / beta ** 2)) ** 0.5))

# clalpha_tail = 2 * np.pi * AR_tail / (
#             2 + ((4 + ((AR_tail * beta / n) ** 2) * (1 + (np.tan(stabilizer_sweep) ** 2) / beta ** 2)) ** 0.5))
# clalpha_tail_degrees = np.radians(clalpha_tail)

# clalpha_tailless = clalpha * (1 + 2.15 * (diameter / span)) * (wing_area - area_fuselage) / wing_area + (
#             np.pi * diameter ** 2) / (2 * wing_area)

# # tail-wing speed ratio
# v_hv = 0

# # data_weight = { "Variable": ['MTOW', 'OEW', 'Payload', 'Fuel weight','Pax & cabin luggage', 'Forward cargo', 'Aft cargo'],
# #                 'CS100': [MTOW, OEW, PL, fuel_weight, pax_cabin, fwd_cargo_max, aft_cargo_max],
# #                 'Percentage of CS100 MTOW':[100, OEW_percent, PL_percent, fuel_weight_percent, pax_cabin_percent, 100*fwd_cargo_max/MTOW, 100*aft_cargo_max/MTOW]
# #                     }

# # table_weight = pd.DataFrame(data_weight, columns = ['Variable', 'Value', 'Percentage of MTOW']).to_latex(caption = 'caption', label = 'label', column_format= '|l|l|l|', multirow = True, float_format="%.2f", index = False)

# # slash = "\\"
# # table_weight = table_weight.replace(slash+slash, slash+slash + ' '+slash+"hline" ).replace("toprule",'hline').replace(".00",'')


# # data_cg = {"CG": ['Forward CG','Aft CG',  'OEW', 'Fuel', 'Forward hold', 'Aft hold', 'Main landing gear'],
# #            "% of MAC": [ extreme_cg[0],extreme_cg[1], cgoew_mac, maccie(cgfuel), maccie(fwdcargo_cg), maccie(aftcargo_cg), maccie(cgmain)],
# #            'Condition':["Aft cargo loaded", 'Fwd cargo loaded', 'Empty', '', '', '', '']}

# # table_cg = pd.DataFrame(data_cg, columns = ["CG", "% of MAC", 'Condition']).to_latex(caption = 'caption', label = 'label', column_format= '|l|l|l|', multirow = True, float_format="%.2f", index = False)

# # slash = "\\"
# # table_cg = table_cg.replace(slash+slash, slash+slash + ' '+slash+"hline" ).replace("toprule",'hline').replace(".00",'')

# # scissor_required = [mach, v_approach, v_hv, clalpha_tailless, 'C_{l'+slash+'alpha_{t-h}}_{wing}','C_{l'+slash+'alpha_{t-h}}_{fuselage}', 'wing downwash gradient','AC_{t-h_{cruise}}',  'AC_{t-h_{cruise_{wing}}}', 'AC_{t-h_{cruise_{fuselage}}}','AC_{t-h_{approach}}', 'AC_{t-h_{approach_{wing}}}', 'AC_{t-h_{approach_{fuselage}}}', 'AC_{t-h_{nacelles}}', 'C_{l_{t}max}', 'C_{m0}' , 'C_{m0_{wing}}', 'C_{m0_{flaps}}' , 'C_{m0_{fuselage}}', 'C_{m0_{nacelles}}']
# # data_scissor = {'Variable':['Cruise speed', 'Approach speed', 'Tail-wing speed ratio', 'C_{l'+slash+'alpha_{t}}', 'C_{l'+slash+'alpha_{t-h}}_{wing}','C_{l'+slash+'alpha_{t-h}}_{fuselage}', 'wing downwash gradient','AC_{t-h_{cruise}}',  'AC_{t-h_{cruise_{wing}}}', 'AC_{t-h_{cruise_{fuselage}}}','AC_{t-h_{approach}}', 'AC_{t-h_{approach_{wing}}}', 'AC_{t-h_{approach_{fuselage}}}', 'AC_{t-h_{nacelles}}', 'C_{l_{t}max}', 'C_{m0}' , 'C_{m0_{wing}}', 'C_{m0_{flaps}}' , 'C_{m0_{fuselage}}', 'C_{m0_{nacelles}}'],
# #                 'New':['Cruise speed', 'Approach speed', 'Tail-wing speed ratio', 'C_{l'+slash+'alpha_{t}}', 'C_{l'+slash+'alpha_{t-h}}_{wing}','C_{l'+slash+'alpha_{t-h}}_{fuselage}', 'wing downwash gradient','AC_{t-h_{cruise}}',  'AC_{t-h_{cruise_{wing}}}', 'AC_{t-h_{cruise_{fuselage}}}','AC_{t-h_{approach}}', 'AC_{t-h_{approach_{wing}}}', 'AC_{t-h_{approach_{fuselage}}}', 'AC_{t-h_{nacelles}}', 'C_{l_{t}max}', 'C_{m0}' , 'C_{m0_{wing}}', 'C_{m0_{flaps}}' , 'C_{m0_{fuselage}}', 'C_{m0_{nacelles}}'],
# #                 'Units':['Cruise speed', 'Approach speed', 'Tail-wing speed ratio', 'C_{l'+slash+'alpha_{t}}', 'C_{l'+slash+'alpha_{t-h}}_{wing}','C_{l'+slash+'alpha_{t-h}}_{fuselage}', 'wing downwash gradient','AC_{t-h_{cruise}}',  'AC_{t-h_{cruise_{wing}}}', 'AC_{t-h_{cruise_{fuselage}}}','AC_{t-h_{approach}}', 'AC_{t-h_{approach_{wing}}}', 'AC_{t-h_{approach_{fuselage}}}', 'AC_{t-h_{nacelles}}', 'C_{l_{t}max}', 'C_{m0}' , 'C_{m0_{wing}}', 'C_{m0_{flaps}}' , 'C_{m0_{fuselage}}', 'C_{m0_{nacelles}}']}

# # table_scissor = pd.DataFrame(data_scissor, columns = ["Variable", "New", 'Units']).to_latex(escape = False, caption = 'caption', label = 'label', column_format= '|l|l|l|', multirow = True, float_format="%.2f", index = False)

# # slash = "\\"
# # table_scissor= table_scissor.replace(slash+slash, slash+slash + ' '+slash+"hline" ).replace("toprule",'hline').replace(".00",'')


# # data_wing = {'Variable':["MAC", "XLEMAC","MAC_{horizontal}", "XLEMAC_{horizontal}","MAC_{vertical}", "XLEMAC_{vertical}", "Horizontal tail arm","Vertical tail arm", "S_{horizontal}", "S_{horizontal}/S_{wing}", "S_{vertical}", "S_{vertical}/S_{wing}", "C_h", "C_v", "X_{nose gear}","X_{main gear}", "X_{main gear_{MAC}}", "b", "b_{horizontal}"],
# #              'CS100':[MAC, xlemac,MAC_horizontal, xlemac_horizontal,MAC_vertical, xlemac_vertical, tail_armh,tail_armv, horizontal_area, horizontal_area/wing_area, vertical_area, vertical_area/wing_area, C_h, C_v, cgnose, cgmain, maccie(cgmain), span, span_tail],
# #              'New':[MAC, xlemac,MAC_horizontal, xlemac_horizontal,MAC_vertical, xlemac_vertical, tail_armh,tail_armv, horizontal_area, horizontal_area/wing_area, vertical_area, vertical_area/wing_area, C_h, C_v, cgnose, cgmain, maccie(cgmain), span, span_tail]}


# # table_wing = pd.DataFrame(data_wing, columns = ["Variable", "CS100", 'New']).to_latex(escape = False, caption = 'caption', label = 'label', column_format= '|l|l|l|', multirow = True, float_format="%.2f", index = False)

# # slash = "\\"
# # table_wing= table_wing.replace(slash+slash, slash+slash + ' '+slash+"hline" ).replace("toprule",'hline').replace(".00",'')

# # print(table_wing)

# # cars = {'Brand': ['Honda Civic','Toyota Corolla','Ford Focus','Audi A4'],
# #         'Price': [22000,25000,27000,35000]
# #         }

# # df = pd.DataFrame(cars, columns = ['Brand', 'Price'])


# # ltx = df.to_latex(caption = 'caption', label = 'label', column_format= '|l|l|l|', multirow = True)
# # slash = "\\"
# # print(ltx.replace(slash+slash, slash+slash + ' '+slash+"hline" ))


# # hundredpercent = OEW_percent + pax_cabin_percent+fuel_weight_percent+cargo_percent

# # percentages = [OEW_percent, pax_cabin_percent,fwd_cargo_percent, aft_cargo_percent, fuel_weight_percent]
# # labels = 'OEW', 'Passengers', 'Fwd Cargo', 'Aft Cargo','Fuel'

# #     #make py chart
# # def pychart():
# #     fig1, ax1 = plt.subplots()
# #     ax1.pie(percentages,  labels=labels, autopct='%1.1f%%', startangle=90)
# #     ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
# #     plt.show()

# # wing_area_ft = 112.3*10.7639104
# # MTOW_lbs = MTOW * 2.2046226218488
# # vertical_area_ft = vertical_area *10.7639104
# # horizontal_area_ft = horizontal_area *10.7639104

# # l_tail = 4.9+2.5
# # l_nose = 34.9/12.5 *2
# # s_wetted = np.pi * diameter * (34.1-l_tail -l_nose ) + (l_nose+l_tail)*diameter*0.5*np.pi


# # weight calculation
# # w_wing = 10*wing_area_ft*0.45359237
# # w_horizontal = 5.5* horizontal_area_ft*0.45359237
# # w_vertical = 5.5 *vertical_area_ft*0.45359237
# # w_fuselage = 5* s_wetted*0.45359237 *10.7639104
# # w_noselandinggear = 0.043*0.15*MTOW
# # w_mainlandinggear = 0.043*0.85*MTOW
# # w_engineinstalled = 1.3*w_engine
# # all_else = 0.17*MTOW

# # calculated_OEW = w_wing+w_horizontal+w_vertical+w_noselandinggear+w_mainlandinggear+w_engineinstalled+all_else +w_fuselage
# # weight_multiplication_factor = OEW/calculated_OEW

# # w_wing = 10*wing_area_ft*0.45359237 * weight_multiplication_factor
# # w_horizontal = 5.5* horizontal_area_ft*0.45359237* weight_multiplication_factor
# # w_vertical = 5.5 *vertical_area_ft*0.45359237* weight_multiplication_factor
# # w_fuselage = 5* s_wetted*0.45359237 *10.7639104* weight_multiplication_factor
# # w_noselandinggear = 0.043*0.15*MTOW* weight_multiplication_factor
# # w_mainlandinggear = 0.043*0.85*MTOW* weight_multiplication_factor
# # w_engineinstalled = 1.3*w_engine* weight_multiplication_factor
# # all_else = 0.17*MTOW* weight_multiplication_factor

# # cg calculations
# # cgwing = 0.4*MAC + xlemac
# # cghorizontal = 0.4*MAC_horizontal+xlemac_horizontal
# # cgvertical = 0.4*MAC_vertical+xlemac_vertical
# # cgfuselage = 0.4*fuselage_lenght
# # cgnose = 3.5
# # cgengine = 10+ 4*10/21
# # cgallelse = 0.4*fuselage_lenght


# # xcgoew = (w_wing*cgwing + w_horizontal * cghorizontal + w_vertical*cgvertical + w_noselandinggear*cgnose + w_mainlandinggear*cgmain + w_engineinstalled*cgengine + all_else*cgallelse + w_fuselage*cgfuselage)/OEW

# # cgoew_mac = 100*(xcgoew-xlemac)/MAC
# # print(cgoew_mac)