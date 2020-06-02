import numpy as np
import input as inp
import matplotlib.pyplot as plt
import Envelope

"""
inputs
    CLmax
    Altitude
    Equivalent Airspeed
    Wing Loading
    Class I Weight Estimation Outputs

outputs 
    Wing Area
    Wing Sweep
    Empennage Area 
    Dihedral Angle
    
Description 
    This script calculates the main aerodynamic values that determine the wing and tail configuration of the aircraft.

"""

# ---------------------------- Import Parameters

M_cruise = 0.8
S = inp.S
AR = inp.AR
MTOW = inp.MTOW

if M_cruise >= 0.7:
    sweep_c4 = np.arccos(0.75*(0.935/(0.03 + M_cruise)))
else:
    sweep_c4 = np.arccos(1)

V_C=Envelope.V_C #Cruise Speed
V_D=Envelope.V_D #Dive Speed
V_S=Envelope.V_S #Stall Speed
V_A=Envelope.V_A #Max Gust Speed

# --------------------------- Equations 
    
taper = 0.2 * (2 - sweep_c4)

b = np.sqrt(S*AR)
c_root = 2*S / ((1 + taper)*b)
c_tip = taper * c_root

c_mac = 2/3 * c_root * (1 + taper + taper**2)/(1 + taper)
y_mac = b/2 * 1/3 * (1 + 2*taper)/(1 + taper)

p_cruise = 101325 * (1 - 0.0065*inp.Cruise_alt*1000/288)**(9.80665/(287*0.0065))
q = 0.5 * 1.4 * p_cruise * M_cruise**2
CL_cruise = MTOW/(q*S)
sweep_c2 = np.arctan(np.tan(sweep_c4) - 4/AR * ((50-25)/100 * (1 - taper)/(1 + taper))) #* 180/np.pi
t_c = min((np.cos(sweep_c2)**3 * (0.935 - (M_cruise + 0.03) * np.cos(sweep_c2)) - 0.115 * CL_cruise**1.5) \
      / (np.cos(sweep_c2)**2), 0.18)

Dihedral = 3
s = sweep_c4*180/np.pi

while s > 0:
    s -= 10
    Dihedral -= 1
Dihedral -= 1






