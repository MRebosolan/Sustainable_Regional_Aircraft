import numpy as np

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


#starting values

MTOW = 60000
AR = 8
half_sweep = np.radians(27)
n_max = 3
n_ult = 3*1.5
S = 80 * 0.3048**2
t_over_c = 0.1 
taper = 0.4 
mach_h = 0.5 

rho = 1.225 
V_dive = 300 
lf = 70
hf = 7
wf = 7
S_fgs = np.pi * lf *(hf+wf)/2

T_TO = 20000





W_wing_gd = wing.gd_wing(MTOW, AR, half_sweep, n_ult, S, t_over_c, taper, mach_h)
W_fuselage_GD = fuselage.W_fuselage_gd (rho, V_dive, MTOW, lf, hf)
# W_nacelle_GD = nacelles.W_nacelle_gd (A_inlet, ln, p2)

W_wing = wing.W_wing(W_zfw, b, half_sweep, n_ult, S, t_r)
W_empennage = tail.vert_tail_weight()+ tail.hor_tail_weight()
W_fuselage = fuselage.W_fuselage_torenbeek(V_dive, lh, wf, hf, S_fgs) 
W_nacelles = nacelles.W_nacelle_torenbeek(T_TO)
landing_gear = LG.LG_weight(Kgr, MTOW, Ag, Bg, Cg, Dg)


W_engines = engine.engine_weight(T_TO, num_engines)

W_fuel_system = fuelsystem.W_fuelsystem (N_t, K_fsp, W_f)

W_power_controls = powercontrols.total(lf, b, W_engines, pneumatic = True)






W_struct = W_wing + W_empennage + W_fuselage + W_nacelles + landing_gear

W_powerplant = W_engines + W_fuel_system + W_power_controls 

W_equipment = APU + cargo_equipment + furnishing + instrumentation + oxygen_system + paint + pressurization + flight_control + electrical









# OEW =

# structural =
# wing_group
# tail_group
# body_group
# landing_group
# nacelle_group

#  #still missing

# propulsion_group =
# engine_installation #matteo

# fuel_system
# power_controls

# equipment=

flight_control_weight = flightcontrols.flight_controls(MTOW)
electrical_system_weight = electrical.electrical_torenbeek(V_pax)
instrumentation_weight = instrumentation.instrumentation_torenbeek(MTOW)
airconditioning_pressurization_weight = pressurization.pressure_system_weight(lpax)
oxygen_system_weight = oxygen.oxygen_system_weight(Npax)
APU_weight = apu.APU_weight_estimation(MTOW)
furnishing_weight = furnishing.furnishing_gd(N_fdc, N_pax, N_cc, MTOW, P_c)
baggage_weight = cargo.cargo_handling_weight(Sff)
paint_weight = paint.paint(MTOW)






