# loading diagram generator

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import Class_2_estimation as cl2
import input

MTOW = cl2.MTOW_kg
OEW = cl2.OEWINPUT
OEW_percent = 100. * OEW / MTOW

Npax = input.Npax
w_person = input.W_pax
pax_cabin = Npax * w_person

max_PL = input.W_payload
# maxpl_fuel = max_PL

cargo = input.W_cargo
fwd_cargo_max = cargo / 3
aft_cargo_max = cargo * 2 / 3

cargo_percent = 100 * cargo / MTOW
fwd_cargo_percent = fwd_cargo_max * 100 / MTOW
aft_cargo_percent = aft_cargo_max * 100 / MTOW

PL = pax_cabin + cargo
PL_percent = 100 * PL / MTOW
fuel_weight = MTOW - OEW - PL
fuel_weight_percent = fuel_weight * 100 / MTOW

M_zfw = OEW + PL

# wing, stabilizer, fuselage and engine parameters
w_engine = cl2.df['SRA'][8]  # kg
wing_area = cl2.S  # m^2
span = input.b
horizontal_area = input.bh
vertical_area = input.bv
diameter = input.widthf
fuselage_lenght = input.lf  # m
cgmain = 0.6 * fuselage_lenght
cgnose = 1

# determining MAC's
xlemac = 0.4 * fuselage_lenght  # m, datum is front of nose
MAC = 3 * 10 / 21  # m, length of MAC
xlemac_horizontal = 63.5 * 10 / 21
MAC_horizontal = 5.3 * 10 / 21
xlemac_vertical = 73 * 34.9 / 87.4
MAC_vertical = 9.5 * 34.9 / 87.4


def maccie(x):
    return 100 * (x - xlemac) / MAC


# cg OEW
xcgoew = 0.4 * fuselage_lenght  # meter
cgoew_mac = 24.95440889952892  # percent
cgoew_fuselage = xcgoew * 100 / fuselage_lenght

# cargo
fwdcargo_begin = 0.1 * fuselage_lenght
fwdcargo_end = fwdcargo_begin + 3
aftcargo_begin = 0.6 * fuselage_lenght
aftcargo_end = aftcargo_begin + 4

fwdcargo_cg = (fwdcargo_begin + fwdcargo_end) / 2
aftcargo_cg = (aftcargo_begin + aftcargo_end) / 2

# fuel cg
cgfuel = xlemac + MAC * 0.4

# passenger cg
seat_start = 4.8 + 0.8 + 0.81 - 1
seatloc = [seat_start]
row = seat_start
pitch = 29  # inch
rows = Npax / 5
for i in range(int(rows - 1)):
    row = row + pitch * 0.0254  # convert to meters
    seatloc.append(row)


# new cg calculation, returns new cg and new weight
def loadingcg(w_old, cg_old, w_item, cg_item):
    x_cg_new = (w_old * cg_old + w_item * cg_item) / (w_old + w_item)
    return x_cg_new, w_old + w_item


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
    onlyfwdcargo = loadingcg(OEW, xcgoew, fwd_cargo_max, fwdcargo_cg)
    bothcargo = loadingcg(onlyfwdcargo[1], onlyfwdcargo[0], aft_cargo_max, aftcargo_cg)
    cargo1 = plt.plot(100 * (np.array([xcgoew, onlyfwdcargo[0], bothcargo[0]]) - xlemac) / MAC,
                      [OEW, onlyfwdcargo[1], bothcargo[1]], label='Cargo', marker='x', color='brown')

    onlyaftcargo = loadingcg(OEW, xcgoew, aft_cargo_max, aftcargo_cg)
    bothcargo2 = loadingcg(onlyaftcargo[1], onlyaftcargo[0], fwd_cargo_max, fwdcargo_cg)
    cargo2 = plt.plot(100 * (np.array([xcgoew, onlyaftcargo[0], bothcargo2[0]]) - xlemac) / MAC,
                      [OEW, onlyaftcargo[1], bothcargo2[1]], marker='x', color='cyan')

    window = passenger_loading(bothcargo[1], bothcargo[0], multiplication=2)
    window_back = passenger_loading(bothcargo[1], bothcargo[0], multiplication=2, seatloc=seatloc[::-1])
    window1 = plt.plot(100 * (np.array(window[0]) - xlemac) / MAC, window[1], label='Window seats', marker='.',
                       color='blue')
    window2 = plt.plot(100 * (np.array(window_back[0]) - xlemac) / MAC, window_back[1], marker='.', color='orange')

    middle = passenger_loading(window[1][-1], window[0][-1], multiplication=2)
    middle_back = passenger_loading(window[1][-1], window[0][-1], multiplication=2, seatloc=seatloc[::-1])
    plt.plot(100 * (np.array(middle[0]) - xlemac) / MAC, middle[1], label='Middle seats', marker='1', color='black')
    plt.plot(100 * (np.array(middle_back[0]) - xlemac) / MAC, middle_back[1], marker='1', color='yellow')

    aisle = passenger_loading(middle[1][-1], middle[0][-1])
    aisle_back = passenger_loading(middle[1][-1], middle[0][-1], seatloc=seatloc[::-1])
    plt.plot(100 * (np.array(aisle[0]) - xlemac) / MAC, aisle[1], label='Aisle seats', marker='2', color='red')
    plt.plot(100 * (np.array(aisle_back[0]) - xlemac) / MAC, aisle_back[1], marker='2', color='green')

    fully_loaded = loadingcg(aisle[1][-1], aisle[0][1], fuel_weight, cgfuel)
    plt.plot(100 * (np.array([aisle[0][-1], fully_loaded[0]]) - xlemac) / MAC, [M_zfw, fully_loaded[1]], marker='^',
             color='magenta', label='Fuel')

    plt.legend()
    plt.grid()
    plt.ylabel('mass [kg]')
    plt.xlabel('xcg [% of MAC]')
    plt.show()
    extreme_cg = maccie(onlyfwdcargo[0]), maccie(onlyaftcargo[0])
    return fully_loaded[0], extreme_cg


cgMTOW, extreme_cg = loading()  # extreme cg is in % of mac, first fwd, then aft, fully loaded is cg

tail_armh = xlemac_horizontal + 0.25 * MAC_horizontal - cgMTOW
tail_armv = xlemac_vertical + 0.25 * MAC_vertical - cgMTOW
C_h = horizontal_area * tail_armh / (wing_area * MAC)
C_v = vertical_area * tail_armv / (wing_area * span)

clalpha = (1 - -1) / 17.5
cmalpha = (0.0329 + 0.0093) / (7.25 + 13.75)
cl0 = -1.25  # degrees
cm0 = -0.0184
area_fuselage = diameter * 14 * (10 / 21)
mach = 0.78
v_approach = 129  # kts
beta = (1 - 0.85 * (mach ** 2)) ** 0.5
n = 0.95
stabilizer_sweep = np.radians(27)
span_tail = 12.3
AR_tail = span_tail ** 2 / horizontal_area
AR = span ** 2 / wing_area
sweep = np.radians(24.5)

clalpha_datcom = 2 * np.pi * AR / (2 + ((4 + ((AR * beta) ** 2) * (1 + (np.tan(sweep) ** 2) / beta ** 2)) ** 0.5))

clalpha_tail = 2 * np.pi * AR_tail / (
            2 + ((4 + ((AR_tail * beta / n) ** 2) * (1 + (np.tan(stabilizer_sweep) ** 2) / beta ** 2)) ** 0.5))
clalpha_tail_degrees = np.radians(clalpha_tail)

clalpha_tailless = clalpha * (1 + 2.15 * (diameter / span)) * (wing_area - area_fuselage) / wing_area + (
            np.pi * diameter ** 2) / (2 * wing_area)

# tail-wing speed ratio
v_hv = 0

# data_weight = { "Variable": ['MTOW', 'OEW', 'Payload', 'Fuel weight','Pax & cabin luggage', 'Forward cargo', 'Aft cargo'],
#                 'CS100': [MTOW, OEW, PL, fuel_weight, pax_cabin, fwd_cargo_max, aft_cargo_max],
#                 'Percentage of CS100 MTOW':[100, OEW_percent, PL_percent, fuel_weight_percent, pax_cabin_percent, 100*fwd_cargo_max/MTOW, 100*aft_cargo_max/MTOW]
#                     }

# table_weight = pd.DataFrame(data_weight, columns = ['Variable', 'Value', 'Percentage of MTOW']).to_latex(caption = 'caption', label = 'label', column_format= '|l|l|l|', multirow = True, float_format="%.2f", index = False)

# slash = "\\"
# table_weight = table_weight.replace(slash+slash, slash+slash + ' '+slash+"hline" ).replace("toprule",'hline').replace(".00",'')


# data_cg = {"CG": ['Forward CG','Aft CG',  'OEW', 'Fuel', 'Forward hold', 'Aft hold', 'Main landing gear'],
#            "% of MAC": [ extreme_cg[0],extreme_cg[1], cgoew_mac, maccie(cgfuel), maccie(fwdcargo_cg), maccie(aftcargo_cg), maccie(cgmain)],
#            'Condition':["Aft cargo loaded", 'Fwd cargo loaded', 'Empty', '', '', '', '']}

# table_cg = pd.DataFrame(data_cg, columns = ["CG", "% of MAC", 'Condition']).to_latex(caption = 'caption', label = 'label', column_format= '|l|l|l|', multirow = True, float_format="%.2f", index = False)

# slash = "\\"
# table_cg = table_cg.replace(slash+slash, slash+slash + ' '+slash+"hline" ).replace("toprule",'hline').replace(".00",'')

# scissor_required = [mach, v_approach, v_hv, clalpha_tailless, 'C_{l'+slash+'alpha_{t-h}}_{wing}','C_{l'+slash+'alpha_{t-h}}_{fuselage}', 'wing downwash gradient','AC_{t-h_{cruise}}',  'AC_{t-h_{cruise_{wing}}}', 'AC_{t-h_{cruise_{fuselage}}}','AC_{t-h_{approach}}', 'AC_{t-h_{approach_{wing}}}', 'AC_{t-h_{approach_{fuselage}}}', 'AC_{t-h_{nacelles}}', 'C_{l_{t}max}', 'C_{m0}' , 'C_{m0_{wing}}', 'C_{m0_{flaps}}' , 'C_{m0_{fuselage}}', 'C_{m0_{nacelles}}']
# data_scissor = {'Variable':['Cruise speed', 'Approach speed', 'Tail-wing speed ratio', 'C_{l'+slash+'alpha_{t}}', 'C_{l'+slash+'alpha_{t-h}}_{wing}','C_{l'+slash+'alpha_{t-h}}_{fuselage}', 'wing downwash gradient','AC_{t-h_{cruise}}',  'AC_{t-h_{cruise_{wing}}}', 'AC_{t-h_{cruise_{fuselage}}}','AC_{t-h_{approach}}', 'AC_{t-h_{approach_{wing}}}', 'AC_{t-h_{approach_{fuselage}}}', 'AC_{t-h_{nacelles}}', 'C_{l_{t}max}', 'C_{m0}' , 'C_{m0_{wing}}', 'C_{m0_{flaps}}' , 'C_{m0_{fuselage}}', 'C_{m0_{nacelles}}'],
#                 'New':['Cruise speed', 'Approach speed', 'Tail-wing speed ratio', 'C_{l'+slash+'alpha_{t}}', 'C_{l'+slash+'alpha_{t-h}}_{wing}','C_{l'+slash+'alpha_{t-h}}_{fuselage}', 'wing downwash gradient','AC_{t-h_{cruise}}',  'AC_{t-h_{cruise_{wing}}}', 'AC_{t-h_{cruise_{fuselage}}}','AC_{t-h_{approach}}', 'AC_{t-h_{approach_{wing}}}', 'AC_{t-h_{approach_{fuselage}}}', 'AC_{t-h_{nacelles}}', 'C_{l_{t}max}', 'C_{m0}' , 'C_{m0_{wing}}', 'C_{m0_{flaps}}' , 'C_{m0_{fuselage}}', 'C_{m0_{nacelles}}'],
#                 'Units':['Cruise speed', 'Approach speed', 'Tail-wing speed ratio', 'C_{l'+slash+'alpha_{t}}', 'C_{l'+slash+'alpha_{t-h}}_{wing}','C_{l'+slash+'alpha_{t-h}}_{fuselage}', 'wing downwash gradient','AC_{t-h_{cruise}}',  'AC_{t-h_{cruise_{wing}}}', 'AC_{t-h_{cruise_{fuselage}}}','AC_{t-h_{approach}}', 'AC_{t-h_{approach_{wing}}}', 'AC_{t-h_{approach_{fuselage}}}', 'AC_{t-h_{nacelles}}', 'C_{l_{t}max}', 'C_{m0}' , 'C_{m0_{wing}}', 'C_{m0_{flaps}}' , 'C_{m0_{fuselage}}', 'C_{m0_{nacelles}}']}

# table_scissor = pd.DataFrame(data_scissor, columns = ["Variable", "New", 'Units']).to_latex(escape = False, caption = 'caption', label = 'label', column_format= '|l|l|l|', multirow = True, float_format="%.2f", index = False)

# slash = "\\"
# table_scissor= table_scissor.replace(slash+slash, slash+slash + ' '+slash+"hline" ).replace("toprule",'hline').replace(".00",'')


# data_wing = {'Variable':["MAC", "XLEMAC","MAC_{horizontal}", "XLEMAC_{horizontal}","MAC_{vertical}", "XLEMAC_{vertical}", "Horizontal tail arm","Vertical tail arm", "S_{horizontal}", "S_{horizontal}/S_{wing}", "S_{vertical}", "S_{vertical}/S_{wing}", "C_h", "C_v", "X_{nose gear}","X_{main gear}", "X_{main gear_{MAC}}", "b", "b_{horizontal}"],
#              'CS100':[MAC, xlemac,MAC_horizontal, xlemac_horizontal,MAC_vertical, xlemac_vertical, tail_armh,tail_armv, horizontal_area, horizontal_area/wing_area, vertical_area, vertical_area/wing_area, C_h, C_v, cgnose, cgmain, maccie(cgmain), span, span_tail],
#              'New':[MAC, xlemac,MAC_horizontal, xlemac_horizontal,MAC_vertical, xlemac_vertical, tail_armh,tail_armv, horizontal_area, horizontal_area/wing_area, vertical_area, vertical_area/wing_area, C_h, C_v, cgnose, cgmain, maccie(cgmain), span, span_tail]}


# table_wing = pd.DataFrame(data_wing, columns = ["Variable", "CS100", 'New']).to_latex(escape = False, caption = 'caption', label = 'label', column_format= '|l|l|l|', multirow = True, float_format="%.2f", index = False)

# slash = "\\"
# table_wing= table_wing.replace(slash+slash, slash+slash + ' '+slash+"hline" ).replace("toprule",'hline').replace(".00",'')

# print(table_wing)

# cars = {'Brand': ['Honda Civic','Toyota Corolla','Ford Focus','Audi A4'],
#         'Price': [22000,25000,27000,35000]
#         }

# df = pd.DataFrame(cars, columns = ['Brand', 'Price'])


# ltx = df.to_latex(caption = 'caption', label = 'label', column_format= '|l|l|l|', multirow = True)
# slash = "\\"
# print(ltx.replace(slash+slash, slash+slash + ' '+slash+"hline" ))


# hundredpercent = OEW_percent + pax_cabin_percent+fuel_weight_percent+cargo_percent

# percentages = [OEW_percent, pax_cabin_percent,fwd_cargo_percent, aft_cargo_percent, fuel_weight_percent]
# labels = 'OEW', 'Passengers', 'Fwd Cargo', 'Aft Cargo','Fuel'

#     #make py chart
# def pychart():
#     fig1, ax1 = plt.subplots()
#     ax1.pie(percentages,  labels=labels, autopct='%1.1f%%', startangle=90)
#     ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
#     plt.show()

# wing_area_ft = 112.3*10.7639104
# MTOW_lbs = MTOW * 2.2046226218488
# vertical_area_ft = vertical_area *10.7639104
# horizontal_area_ft = horizontal_area *10.7639104

# l_tail = 4.9+2.5
# l_nose = 34.9/12.5 *2
# s_wetted = np.pi * diameter * (34.1-l_tail -l_nose ) + (l_nose+l_tail)*diameter*0.5*np.pi


# weight calculation
# w_wing = 10*wing_area_ft*0.45359237
# w_horizontal = 5.5* horizontal_area_ft*0.45359237
# w_vertical = 5.5 *vertical_area_ft*0.45359237
# w_fuselage = 5* s_wetted*0.45359237 *10.7639104
# w_noselandinggear = 0.043*0.15*MTOW
# w_mainlandinggear = 0.043*0.85*MTOW
# w_engineinstalled = 1.3*w_engine
# all_else = 0.17*MTOW

# calculated_OEW = w_wing+w_horizontal+w_vertical+w_noselandinggear+w_mainlandinggear+w_engineinstalled+all_else +w_fuselage
# weight_multiplication_factor = OEW/calculated_OEW

# w_wing = 10*wing_area_ft*0.45359237 * weight_multiplication_factor
# w_horizontal = 5.5* horizontal_area_ft*0.45359237* weight_multiplication_factor
# w_vertical = 5.5 *vertical_area_ft*0.45359237* weight_multiplication_factor
# w_fuselage = 5* s_wetted*0.45359237 *10.7639104* weight_multiplication_factor
# w_noselandinggear = 0.043*0.15*MTOW* weight_multiplication_factor
# w_mainlandinggear = 0.043*0.85*MTOW* weight_multiplication_factor
# w_engineinstalled = 1.3*w_engine* weight_multiplication_factor
# all_else = 0.17*MTOW* weight_multiplication_factor

# cg calculations
# cgwing = 0.4*MAC + xlemac
# cghorizontal = 0.4*MAC_horizontal+xlemac_horizontal
# cgvertical = 0.4*MAC_vertical+xlemac_vertical
# cgfuselage = 0.4*fuselage_lenght
# cgnose = 3.5
# cgengine = 10+ 4*10/21
# cgallelse = 0.4*fuselage_lenght


# xcgoew = (w_wing*cgwing + w_horizontal * cghorizontal + w_vertical*cgvertical + w_noselandinggear*cgnose + w_mainlandinggear*cgmain + w_engineinstalled*cgengine + all_else*cgallelse + w_fuselage*cgfuselage)/OEW

# cgoew_mac = 100*(xcgoew-xlemac)/MAC
# print(cgoew_mac)