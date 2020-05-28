# fuselage
from math import *


def dynamic_pressure(rho, V):
    q = 0.5*rho* V*V
    q_psf = q * 0.020885
    return q, q_psf

def W_fuselage_gd (rho, V_dive, MTOW, lf, hf): #V_dive in m/s
    q_d = dynamic_pressure(rho, V_dive)[1]
    return 10.43 * ((q_d/100)**0.283) * ((MTOW/1000)**0.95)*((lf/hf)**0.71) # q is design dive dynamic pressure, MTOW in lbs, lf and hf in feet

def W_fuselage_torenbeek(V_d, lh, wf, hf, S_fgs): #v_dive in keas, lh = distance from wing root quarter chord to horizontal tail root quarter chord
    #sfgs = fuselage gross shell area in square feet
    kf = 1.08 # for a pressurized fuselage
    return 0.021*kf * sqrt((V_d*lh /(wf+hf))) * S_fgs**1.2




