import APU_weight_estimation as apu
import Cargo_handling_weight_estimation as cargo
import electrical_system_weight_estimation as electrical
import Engine_weight_estimation as engine
import Envelope as envelope
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
from class_1_estimation import CLASS1WEIGHTHYBRID
import input


def to_pounds(kg):
    return kg * 2.20462262

def to_kg (lbs):
    return lbs/2.20462262
#------------- INPUT PARAMETERS ----------------#

AR = input.AR
half_sweep = input.half_sweep
n_max =  input.n_max
n_ult = 1.5* n_max

S = input.S
t_over_c = input.t_over_c
taper = input.taper	
mach_h = input.mach_h
rho = 1.225 #estimate, still convert to imperial!!!
V_dive = input.V_dive
lf = input.lf
hf =input.hf
A_inlet = input.A_inlet
ln = input.ln
p2 = input.p2
# W_zfw = input.W_zfw
b = input.b
t_r= input.t_r
wf = input.wf #max fuselage width
S_fgs = input.S_fgs #fuselage gross shell area
lh =input.lh
T_TO = input.T_TO

Kgr = input.Kgr
V_pax = input.V_pax
lpax =input.lpax
Npax =input.Npax
N_fdc =input.N_fdc
N_cc =input.N_cc
P_c = input.P_c
Sff = input.Sff
T_dry_SL = T_TO

N_eng = input.N_eng
N_t = input.N_t
K_fsp = input.K_fsp
W_fuel = input.W_fuel


ratio = input.hydrogenratio

iterate = 0
OEWINPUT = 0

while abs((OEW_class1 - OEW_class2)*100/OEW_class2)>= 0.5 and iterate < 5000:
    class1 = CLASS1WEIGHTHYBRID(ratio,OEWINPUT)
    MTOW_kg = class1[0]
    OEW_class1_kg = class1[1]
    M_zfw_kg = class1[4]
    # M_fuel_kg = class1[2]
    
    OEW_class1 = to_pounds(OEW_class1_kg)
    MTOW = to_pounds(MTOW_kg)

    
#--------- STRUCTURAL WEIGHT --------------#

W_wing_gd = wing.gd_wing(MTOW, AR, half_sweep, n_ult, S, t_over_c, taper, mach_h)
W_fuselage_GD = fuselage.W_fuselage_gd (rho, V_dive, MTOW, lf, hf)
W_nacelle_GD = nacelles.W_nacelle_gd (A_inlet, ln, p2)

W_wing = wing.W_wing(W_zfw, b, half_sweep, n_ult, S, t_r)
W_empennage = tail.vert_tail_weight()+ tail.hor_tail_weight()
W_fuselage = fuselage.W_fuselage_torenbeek(V_dive, lh, wf, hf, S_fgs)
W_nacelles = nacelles.W_nacelle_torenbeek(T_TO)
W_landing_gear = LG.LG_weight(Kgr, MTOW, Ag, Bg, Cg, Dg)


#--------- EQUIPMENT WEIGHT ------------#

flight_control_weight = flightcontrols.flight_controls(MTOW)
electrical_system_weight = electrical.electrical_torenbeek(V_pax)
instrumentation_weight = instrumentation.instrumentation_torenbeek(MTOW)
airconditioning_pressurization_weight = pressurization.pressure_system_weight(lpax)
oxygen_system_weight = oxygen.oxygen_system_weight(Npax)
APU_weight = apu.APU_weight_estimation(MTOW)
furnishing_weight = furnishing.furnishing_gd(N_fdc, Npax, N_cc, MTOW, P_c)
cargo_equipment_weight = cargo.cargo_handling_weight(Sff)
paint_weight = paint.paint(MTOW)

#------------ POWER PLANT WEIGHT ------------#

W_engines = engine.engine_weight(T_dry_SL, N_eng)
# W_fuel_system = fuelsystem.W_fuelsystem (N_t, K_fsp, W_fuel)
W_fuel_system = to_pounds(class1[-1])
W_power_controls = powercontrols.total(lf, b, W_engines, pneumatic = True)




#------------- SUMMATION OF COMPONENTS -------------#

W_struct = W_wing + W_empennage + W_fuselage + W_nacelles + W_landing_gear

W_powerplant = W_engines + W_fuel_system + W_power_controls

W_equipment = APU_weight + cargo_equipment_weight + furnishing_weight + instrumentation_weight + oxygen_system_weight + paint_weight + airconditioning_pressurization_weight + flight_control_weight + electrical_system_weight

OEW_class2 = W_struct + W_powerplant + W_equipment