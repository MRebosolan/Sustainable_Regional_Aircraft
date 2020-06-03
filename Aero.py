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
    Aerodynamic Requirements
    Wing Area
    
outputs     
    Wing Sweep
    CL cruise
    Empennage Area 
    Dihedral Angle
    Airfoil Geometry
    
description 
    This script calculates the main aerodynamic values that determine the wing and tail configuration of the aircraft.
"""

# ---------------------------- Import Parameters

M_cruise = 0.75
S = inp.S
AR = inp.AR
MTOW = inp.MTOW

V_C = Envelope.V_C  # Cruise Speed
V_D = Envelope.V_D  # Dive Speed
V_S = Envelope.V_S  # Stall Speed
V_A = Envelope.V_A  # Max Gust Speed

def wing_geometry(M_cruise, S, AR, MTOW, V_C):

    if M_cruise >= 0.7:
        sweep_c4 = np.arccos(0.75*(0.935/(0.03 + M_cruise)))
    else:
        sweep_c4 = np.arccos(1)

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

    dihedral = 3
    s = sweep_c4*180/np.pi

    while s > 0:
        s -= 10
        dihedral -= 1
    dihedral -= 1

    sweep_cLE = np.arctan(np.tan(sweep_c4) - 4 / AR * ((0 - 25) / 100 * (1 - taper) / (1 + taper)))

    x_wing = [0, b/2, b/2, 0, 0]
    y_wing = [0, -b/2*np.tan(sweep_cLE), (-b/2*np.tan(sweep_cLE) - c_tip), - c_root, 0]

    plt.plot(x_wing, y_wing)
    plt.plot()
    plt.show()

    AR_check = 17.7 * (2 - taper) * np.exp(- 0.043 * sweep_c4/np.pi*180)
    print(AR_check)

    WS_cr_start = 0.9843800695598843 * MTOW / S

    WS_cr_end = 0.9629656887889539 * MTOW / S

    CL_des = 1.1/q * (0.5*(WS_cr_start + WS_cr_end))
    Cl_des = CL_des / np.cos(sweep_c4)**2
    print("Cl design =", CL_des, Cl_des)

    T_alt = 288 * (1 - 0.0065*inp.Cruise_alt*1000/288)
    visc_k = 9.2e-6

    Re = V_C * 0.514444 * c_mac / visc_k
    print(Re)
    # With CL_max = 1.8 we could take airfoil NACA 63(3)-618 (supercritical with 0.18 t/c)

    # CLmax take-off: 2.1 , Clmax landing: 2.25
    #     # target for take-off: Delta CLmax = 0.3
    #     # target for landing: Delta CLmax = 0.45

    dCLmax_land = 0.45
    dCLmax_to   = 0.3

    print(dCLmax_land)

    return sweep_c4, taper, c_root, c_tip, c_mac, y_mac, t_c, dihedral, Cl_des, dCLmax_land, dCLmax_to

wing_geometry(M_cruise, S, AR, MTOW, V_C)


def airfoilplot(datfile):
    f=open('airfoil2','r')
    lines=f.readlines(f.split('/n '))
#    xy=lines.split('/n')
#    xcoord=[]
#    ycoord1=[]
#    ycoord2=[]
#    for i in range(0,len(xy)):
#        xcoord.append(float(xy[i][0]))
#        ycoord1.append(float(xy[i][1]))
#    return xcoord, ycoord1, ycoord2
#
##def airfoilplot('airfoil1.txt'):
#f=open('airfoil1.txt','r')
#lines=f.readlines()
#result=[]
#for x in lines:
#    result.append(x.split('\n'))
#f.close()
#    return result

#    f=open('datfile','r')
#    lines=f.readlines()
#    xy=lines.split('/n')
#    xcoord=[]
#    ycoord1=[]
#    ycoord2=[]
#    for i in range(0,len(xy)):
#        xcoord.append(float(xy[i][0]))
#        ycoord1.append(float(xy[i][1]))
#    return xcoord, ycoord1, ycoord2




