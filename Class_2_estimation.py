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



W_wing_gd = wing.gd_wing(MTOW, AR, half_sweep, n_ult, S, t_over_c, taper, mach_h)
W_fuselage_GD = fuselage.W_fuselage_gd (rho, V_dive, MTOW, lf, hf)
W_nacelle_GD = nacelle.W_nacelle_gd (A_inlet, ln, p2)

W_wing = wing.W_wing(W_zfw, b, half_sweep, n_ult, S, t_r)
W_empennage = tail.vert_tail_weight()+ tail.hor_tail_weight()
W_fuselage = fuselage.W_fuselage_torenbeek(V_d, lh, wf, hf, S_fgs) 
W_nacelles = nacelles.W_nacelle_torenbeek(T_TO)
landing_gear = LG.LG_weight(Kgr, Wto, Ag, Bg, Cg, Dg)


W_engines = engine.engine_weight(dry_thrust_SL, num_engines)

W_fuel_system = fuelsystem.W_fuelsystem (N_t, K_fsp, W_f)

W_power_controls = powercontrols.total(lf, b, W_e, pneumatic = True)






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

# flight_control
# hydraulics #see appendix A
# electrical system
# instrumentation #done
# airconditioning_pressurization# done
# oxygen_system #done
# APU # matteo
# furnishing #jorn
# baggage_cargo
# operational_items #included in furnishing
# ballast #tbd
# paint

# misc


# # furnishing =
# # flight_deck
# # cabin_accomodations
# # emergency_equipment


# # operational_items =
# # crew
# # supplies
# # water
# # safety_equipment
# # trapped_fuel
# # cargo_equipment

