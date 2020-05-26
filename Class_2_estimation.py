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

















W_struct = wing + empennage + fuselage + nacelles + landing_gear

W_powerplant = engines + fuel_system + power_controls 

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

