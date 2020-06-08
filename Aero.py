import numpy as np
import input as inp
import matplotlib.pyplot as plt
import Envelope
import Class_2_estimation as CL2
#import Aileron_sizing


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
S = CL2.S
AR = inp.AR
MTOW = inp.MTOW
widthf = inp.widthf

V_C = inp.V_C  # Cruise Speed knots
V_D = inp.V_dive  # Dive Speed knots
V_S = inp.V_S  # Stall Speed knots
V_A = inp.V_A  # Max Gust Speed knots
v_approach = inp.v_approach                  # approach speed m/s
V_C_TAS = inp.V_C_TAS    # True air speed cruise m/s

b1 = 8
b2 = 11.2

# ---------------------------- Line Intersection Point

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

# --------------------------------- Chord length at a spanwise point x

def chord_length(c_root, c_tip, x, b):
    l_chord = -2 * ((c_root - c_tip)/b) * x + c_root
    return l_chord

# --------------------------------- Wing Geometry

def wing_geometry(M_cruise, S, AR, MTOW, V_C, widthf, V_S, v_approach, V_C_TAS):

    # Wing sweep, also consider M_crit?
    if M_cruise >= 0.7:
        sweep_c4 = np.arccos(0.75*(0.935/(0.03 + M_cruise)))
    else:
        sweep_c4 = np.arccos(1)

    sweep_c4 = 40.535802011 * np.pi / 180 # from airfoil selection

    print("Sweep =", sweep_c4 * 180 / np.pi)

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
          / (np.cos(sweep_c2)**2), 0.18)  #Upper limit for wing thickness

    print("t/c = ", t_c)

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
    Cl_des = Cl_des * np.sqrt(1 - M_cruise**2)
    print("Cl design =", CL_des, Cl_des)

    T_alt = 288 * (1 - 0.0065*inp.Cruise_alt*1000/288)
    mu = 1.458e-6 * T_alt**1.5 / (T_alt + 110.4)
    rho = p_cruise / (287 * T_alt)

    Re = (rho * V_C_TAS * c_mac) / mu
    M_cruise2 = V_C_TAS / np.sqrt(1.4 * 287 * T_alt)
    Re2 = (rho * M_cruise * np.sqrt(1.4 * 287 * T_alt)  * c_mac) / mu
    print("V_C=", V_C, V_C_TAS)

    Re_sea = (1.225 * 66 * c_mac) / 1.802e-5

    M_sea = v_approach / np.sqrt(1.4*287*288)
    print("M_sea = ", M_sea, M_cruise2)

    print("Re =", Re, Re2, Re_sea)
    # With CL_max = 1.8
    # CLmax = 2.464
    # CLmax take-off: 2.1 , Clmax landing: 2.25
    #     # target for take-off: Delta CLmax = 0.3
    #     # target for landing: Delta CLmax = 0.45

    k = (0.115 * 180 / np.pi) / (2 * np.pi)

    beta = np.sqrt(1 - M_cruise**2)
    CLalpha = (2*np.pi*AR)/ (2 + np.sqrt(4 + (AR * beta / k) * (1 + (np.tan(sweep_c2)**2/beta**2))))

    alpha0L = -3.667/180 * np.pi

    #CL = CLalpha * (alpha - alpha0L)

    alpha_trim = CL_des / CLalpha + alpha0L

    CLmax_corrected = 0.74 * 2.432
    print("CLmax_corrected =", CLmax_corrected)

    dCLmax_land = 0.45
    dCLmax_to   = 0.3

    dClmax_land = 0.9

    hinge_c     = 75 #percent
    aileron_C   = 75 #percent

    sweep_hinge = np.arctan(np.tan(sweep_c4) - 4/AR * ((hinge_c-25)/100 * (1 - taper)/(1 + taper)))
    print(sweep_hinge)
    SwfS = dCLmax_land/ (0.9 * dClmax_land * np.cos(sweep_hinge))

    Df = widthf/2 * 1

    a = -2 * (c_root - c_tip)/b
    ch = 1 - (hinge_c/100)

    D = (-4 * a * Df + 2 * c_root)**2 + 4 * 2 * a * (-2 * SwfS * S)

    print("D = ", D)

    x2 = max((-(-4 * a * Df + 2 * c_root) + np.sqrt(D))/ (-4 * a), (-(-4 * a * Df + 2 * c_root) - np.sqrt(D)) / (-4 * a))
    print("x2 = ", (x2 + Df))

    Sw_check = ((-2 * a * Df + c_root) + (-2 * a * (Df + x2) + c_root)) * x2 / 2 / S   #verified
    print(SwfS, Sw_check)

    d_alpha_0 = -5 * np.pi / 180

    d_alpha = d_alpha_0 * SwfS * np.cos(sweep_hinge)


    wing = [sweep_c4, taper, c_root, c_tip, c_mac, y_mac, t_c, dihedral,
            Cl_des, dCLmax_land, dCLmax_to]

    print("wing =", wing)

    cross1 = line_intersect(x_fus[0],y_fus[0],x_fus[1],y_fus[1],x_wing[0],y_wing[0],x_wing[1],y_wing[1])

    x_hld = [Df, Df, x2 + Df, x2 + Df]
    y_hld = [(-Df*np.tan(sweep_cLE) - chord_length(c_root, c_tip, Df, b)),
             (-Df * np.tan(sweep_cLE) - (hinge_c/100) * chord_length(c_root, c_tip, Df, b)),
             (-(Df + x2) * np.tan(sweep_cLE) - (hinge_c/100) * chord_length(c_root, c_tip, (Df + x2), b)),
             (-(Df + x2) * np.tan(sweep_cLE) - chord_length(c_root, c_tip, (Df + x2), b))]
    
    x_ail = [b1, b1, b2, b2]
    y_ail = [(-b1*np.tan(sweep_cLE) - chord_length(c_root, c_tip, b1, b)),
             (-b1 * np.tan(sweep_cLE) - (hinge_c/100) * chord_length(c_root, c_tip, b1, b)),
             (-b2 * np.tan(sweep_cLE) - (hinge_c/100) * chord_length(c_root, c_tip, b2, b)),
             (-b2 * np.tan(sweep_cLE) - chord_length(c_root, c_tip, b2, b))]
    
    ail   = [x_ail,y_ail]
    hld   = [x_hld, y_hld]



    return wing, geom,cross1, hld, ail, x2


wing, geom, cross1, hld, ail, x2 = wing_geometry(M_cruise, S, AR, MTOW, V_C, widthf, V_S, v_approach, V_C_TAS)


#----------------------------- .txt File Airfoil Coordinates

#Read from file

#Create empty lists

lines  = [[],[],[],[]]
xcoord1= [[],[],[],[]]
xcoord2= [[],[],[],[]]
ycoord1= [[],[],[],[]]
ycoord2= [[],[],[],[]]
camline= [[],[],[],[]]

#Read from file
f0=open('airfoil1.txt','r')
f1=open('airfoil2.txt','r')
f2=open('airfoil3.txt','r')
f3=open('airfoil4.txt','r')

lines0=f0.readlines()
lines1=f1.readlines()
lines2=f2.readlines()
lines3=f3.readlines()

for i in range(0,103):
    xcoord1[0].append(float(lines0[i].split()[0]))
    xcoord1[1].append(float(lines1[i].split()[0]))
    xcoord1[2].append(float(lines2[i].split()[0]))
    ycoord1[0].append(float(lines0[i].split()[1]))
    ycoord1[1].append(float(lines1[i].split()[1]))
    ycoord1[2].append(float(lines2[i].split()[1]))        
for i in range(0,70):
    xcoord1[3].append(float(lines3[i].split()[0]))
    ycoord1[3].append(float(lines3[i].split()[1]))
    
for i in range(103,205):
    xcoord2[0].append(float(lines0[i].split()[0]))
    xcoord2[1].append(float(lines1[i].split()[0]))
    xcoord2[2].append(float(lines2[i].split()[0]))
    ycoord2[0].append(float(lines0[i].split()[1]))
    ycoord2[1].append(float(lines1[i].split()[1]))
    ycoord2[2].append(float(lines2[i].split()[1]))        
for i in range(70,139):
    xcoord2[3].append(float(lines3[i].split()[0]))
    ycoord2[3].append(float(lines3[i].split()[1]))

for i in range(0,4):
    xcoord2[i].insert(-1,1.0)
    ycoord2[i].insert(-1,0.0)
    xcoord1[i]=xcoord1[i][::-1]
    ycoord1[i]=ycoord1[i][::-1]

##Camber Line
for i in range(0,len(xcoord1[0])):
    camline[0].append((ycoord1[0][i]+ycoord2[0][i])/2)
    camline[1].append((ycoord1[1][i]+ycoord2[1][i])/2)
    camline[2].append((ycoord1[1][i]+ycoord2[1][i])/2)
for i in range(0,len(xcoord1[3])):
    camline[3].append((ycoord1[3][i]+ycoord2[3][i])/2)


#----------------------------- Plotting

#plt.figure(0)
#plt.plot(geom[0], geom[1], geom[2], geom[3], hld[0], hld[1],ail[0],ail[1])
#plt.text(cross1[0],cross1[1],'Fuselage Wall Line')
#plt.grid(True,which="major",color="#999999")
#plt.grid(True,which="minor",color="#DDDDDD",ls="--")
#plt.minorticks_on()
#plt.ylim(-12.0,2.0)
#plt.ylabel('x [m]')
#plt.xlabel('y [m]')
#
plt.figure(1)
plt.grid(True,which="major",color="#999999")
plt.grid(True,which="minor",color="#DDDDDD",ls="--")
plt.minorticks_on()
plt.plot(xcoord1[0],ycoord1[0],color='r')
plt.plot(xcoord2[0],camline[0],'--',color='r')
plt.plot(xcoord2[0],ycoord2[0],color='r')
plt.xlim(0,1)
plt.ylim(-0.3,0.3)
plt.text(0.0,0.0,'LE')
plt.text(1.0,0.0,'TE')
plt.ylabel('y/c [-]')
plt.xlabel('x/c [-]')

#plt.figure(2)
#plt.grid(True,which="major",color="#999999")
#plt.grid(True,which="minor",color="#DDDDDD",ls="--")
#plt.minorticks_on()
#plt.plot(xcoord1[1],ycoord1[1],color='r')
#plt.plot(xcoord2[1],camline[1],'--',color='r')
#plt.plot(xcoord2[1],ycoord2[1],color='r')
#plt.xlim(0,1)
#plt.ylim(-0.3,0.3)
#plt.text(0.0,0.0,'LE')
#plt.text(1.0,0.0,'TE')
#plt.ylabel('y/c [-]')
#plt.xlabel('x/c [-]')
#
#plt.figure(3)
#plt.grid(True,which="major",color="#999999")
#plt.grid(True,which="minor",color="#DDDDDD",ls="--")
#plt.minorticks_on()
#plt.plot(xcoord1[2],ycoord1[2],color='r')
#plt.plot(xcoord2[2],camline[2],'--',color='r')
#plt.plot(xcoord2[2],ycoord2[2],color='r')
#plt.xlim(0,1)
#plt.ylim(-0.3,0.3)
#plt.text(0.0,0.0,'LE')
#plt.text(1.0,0.0,'TE')
#plt.ylabel('y/c [-]')
#plt.xlabel('x/c [-]')
#
#plt.figure(4)
#plt.grid(True,which="major",color="#999999")
#plt.grid(True,which="minor",color="#DDDDDD",ls="--")
#plt.minorticks_on()
#plt.plot(xcoord1[3],ycoord1[3],color='r')
#plt.plot(xcoord2[3],camline[3],'--',color='r')
#plt.plot(xcoord2[3],ycoord2[3],color='r')
#plt.xlim(0,1)
#plt.ylim(-0.3,0.3)
#plt.text(0.0,0.0,'LE')
#plt.text(1.0,0.0,'TE')
#plt.ylabel('y/c [-]')
#plt.xlabel('x/c [-]')
#
plt.show()





