import numpy as np
import input as inp
import matplotlib.pyplot as plt

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

M_cruise = 0.6
S = inp.S
AR = inp.AR

if M_cruise >= 0.7:
    sweep_c4 = np.acos(0.935/(0.03 + M_cruise))
else:
    sweep_c4 = np.acos(1)

taper = 0.2 * (2 - sweep_c4)

b = np.sqrt(S*AR)
c_root = 2*S / ((1 + taper)*b)
c_tip = taper * c_root