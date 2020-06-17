import input
import matplotlib.pyplot as plt
import Class_2_estimation as C2

C2.OEM = C2.OEW_class2
C2.MTOW = C2.MTOW/input.g

labels = 'Wing', 'Empennage', 'Fuselage', 'Propulsion system ','Landing gear', 'Electrical system','Instruments and controls', 'Other'
sizes = [C2.W_wing/C2.OEM, C2.W_empennage/C2.OEM, C2.W_fuselage/C2.OEM, C2.W_powerplant/C2.OEM,(C2.flight_control_weight+C2.instrumentation_weight)/C2.OEM,C2.W_landing_gear/C2.OEM,C2.electrical_system_weight/C2.OEM,(C2.OEM-(C2.flight_control_weight+C2.instrumentation_weight)-C2.electrical_system_weight-C2.W_landing_gear-C2.W_powerplant-C2.W_fuselage-C2.W_wing-C2.W_empennage)/C2.OEM]
labels2 = 'Payload', 'Fuel','OEW'
sizes2 = [C2.W_payload/C2.MTOW,C2.W_fuel/C2.MTOW,C2.OEM/C2.MTOW]


fig1, ax1 = plt.subplots()
ax1.pie(sizes, labels=labels, autopct='%1.1f%%', startangle=-15)
ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.

fig2,ax2 = plt.subplots()
ax2.pie(sizes2,labels=labels2, autopct='%1.1f%%', startangle=110)
ax2.axis('equal')
plt.show()
