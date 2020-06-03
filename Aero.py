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
widthf = inp.widthf

V_C = Envelope.V_C  # Cruise Speed
V_D = Envelope.V_D  # Dive Speed
V_S = Envelope.V_S  # Stall Speed
V_A = Envelope.V_A  # Max Gust Speed

def wing_geometry(M_cruise, S, AR, MTOW, V_C, widthf):

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
    x_fus = [widthf/2, widthf/2]
    y_fus = [0, -6]

    geom = [x_wing, y_wing, x_fus, y_fus]

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

    dClmax_land = 1.3
    hinge_c     = 80 #percent
    sweep_hinge = np.arctan(np.tan(sweep_c4) - 4/AR * ((hinge_c-25)/100 * (1 - taper)/(1 + taper)))

    SwfS = dCLmax_land/ (0.9 * dClmax_land * np.cos(sweep_hinge))

    x1 = widthf/2

    a = -2 * (c_root - c_tip)/b

    D = (2 * (a * x1 + c_root))**2 + 4 * a * SwfS * S

    x2 = max((-2 * (a * x1 + c_root) + np.sqrt(D))/(2*a),  (-2 * (a * x1 + c_root) - np.sqrt(D))/(2*a))
    print(x2)

    wing = [sweep_c4, taper, c_root, c_tip, c_mac, y_mac, t_c, dihedral,
            Cl_des, dCLmax_land, dCLmax_to]

    return wing, geom


wing, geom = wing_geometry(M_cruise, S, AR, MTOW, V_C, widthf)

hld(dCLmax_land, dCLmax_to, sweep_c4, taper)

# -------- Line Intersection Point--------------------

def line_intersect(Ax1, Ay1, Ax2, Ay2, Bx1, By1, Bx2, By2):
    """ returns a (x, y) tuple or None if there is no intersection """
    d = (By2 - By1) * (Ax2 - Ax1) - (Bx2 - Bx1) * (Ay2 - Ay1)
    if d:
        uA = ((Bx2 - Bx1) * (Ay1 - By1) - (By2 - By1) * (Ax1 - Bx1)) / d
        uB = ((Ax2 - Ax1) * (Ay1 - By1) - (Ay2 - Ay1) * (Ax1 - Bx1)) / d
    else:
        return
    if not(0 <= uA <= 1 and 0 <= uB <= 1):
        return
    x = Ax1 + uA * (Ax2 - Ax1)
    y = Ay1 + uA * (Ay2 - Ay1)

    return x, y

cross1 = line_intersect(x_fus[0],yfus[0],xfus[1],y_fus[1],0.0,0.0)

#----------------------------- .txt File Airfoil Coordinates

f=open('airfoil2.txt','r')
lines=f.readlines()
xcoord1=[]
xcoord2=[]
ycoord1=[]
ycoord2=[]
for i in lines[:26]:
    xcoord1.append(float(i.split('     ')[0]))
    ycoord1.append(float(i.split('     ')[1].strip('\n')))
for i in lines[26:]:
    xcoord2.append(float(i.split('     ')[0]))
    ycoord2.append(float(i.split('     ')[1].strip('\n')))

xcoord2.insert(0,0.0)
ycoord2.insert(0,0.0)

print(lines)
print(xcoord1)
print(ycoord1)
print(ycoord2)

plt.figure(0)
plt.plot(geom[0], geom[1], geom[2], geom[3])
plt.plot()
plt.show()


#----------------------------- Plotting

plt.figure(1)
plt.grid(True,which="major",color="#999999")
plt.grid(True,which="minor",color="#DDDDDD",ls="--")
plt.minorticks_on()
plt.plot(xcoord1,ycoord1,color='r')
plt.plot(xcoord2,ycoord2,color='r')
plt.xlim(0,1)
plt.ylim(-0.3,0.3)
plt.text(0.0,0.0,'LE')
plt.text(1.0,0.0,'TE')

plt.show()





