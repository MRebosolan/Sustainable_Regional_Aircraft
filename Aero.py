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

M_cruise = inp.mach_cruise
S = CL2.S
print(S)
AR = inp.AR
MTOW = inp.MTOW
widthf = inp.widthf

V_C = inp.V_C  # Cruise Speed knots
V_D = inp.V_dive  # Dive Speed knots
V_S = inp.V_S  # Stall Speed knots
V_A = inp.V_A  # Max Gust Speed knots
v_approach = inp.v_approach                  # approach speed m/s
V_C_TAS = inp.V_C_TAS    # True air speed cruise m/s

b1 = 7.6               #aileron inside y position , starts where flap ends
b2 = b1 + 2.707943    # The roll rate requirement is met with a difference of 9.196882743367496e-06 [deg/s]

############ for Drag #############

k = inp.k
CL_minD = inp.Cd0_clean   # Cd0 of the airfoil
IF_wing   = inp.IF_wing
IF_tailv  = inp.IF_tailv
IF_tailh  = inp.IF_tailh
IF_fus    = inp.IF_fus
IF_nacelle = inp.IF_nacelle
cds_nose   = inp.cds_nose
dflap = inp.dflap

#### TAIL INPUTS
# Vertical tail: Tobias,     horizontal tail: jorn
Sh = 1       # horizontal tail area
Sv = 15       # vertical tail area                                        # WILL CHANGE
taperh = inp.taper_h   # taper h tail
taperv = inp.taper_v   # taper v tail
AR_h   = inp.AR_h
AR_v   = inp.AR_v
sweep_c4h = np.arctan(np.tan(inp.half_chord_sweep_hor) - 4 / AR_h * ((25 - 50) / 100 * (1 - taperh) / (1 + taperh)))
sweep_c4v = np.arctan(np.tan(inp.half_chord_sweep_vert) - 4 / AR_v * ((25 - 50) / 100 * (1 - taperv) / (1 + taperv)))
c_MACh = inp.c_mac_h   # MAC length horizontal tail
c_MACv = 3            # MAC length vertical tail                          # THIS IS GONNA CHANGE
x_cm_wing = 0.36       #x/c max thickness
x_cm_tailv = 0.30         #x/c max thickness vertical tail
x_cm_tailh = 0.30         #x/c max thickness horizontal tail               # Check
t_c_wing = 0.14        #t/c wing airfoil
t_c_tailv = 0.12           #t/c vertical tail airfoil
t_c_tailh = 0.12           #t/c horizontal tail airfoil                   # Check

#### FUSELAGE INPUTS
Df = widthf
lf = inp.lf          # length fuselage
print("lf=",lf)
L1 = 0.2 * lf               # nosecone length L1,  guesss
L2 = 0.6 * lf               # cabin length    L2,  guesss
L3 = 0.2 * lf               # tail length     L3,  guesss
Amax_fus = inp.A_fuselage                    # estimate
upsweep = inp.theta * np.pi / 180                       # Clearance angel [deg]
#### ENGINE
l_nacelle = inp.ln     # nacelle length l_nacelle
Amax_nacelle = np.pi / 4 * inp.bn**2   # max area nacelle Amax_nacelle
B = 10                # Bypass ratio, take from inp
Tto = inp.Tto          # take-off thrust

#### LANDING GEAR
d_nose = inp.d_wheel_nose_lg + inp.strut_length_nose_lg    # length nose gear
w_nose = 1                                              # total width nose gear
d_main = inp.d_wheel_main_lg + inp.strut_length_main_lg         # length main gear
w_main = 1                                              # total width main gear
main_amount = 2                                         # Main landing gear amount
S_mlg = d_main * w_main                                 # reference frontal area main landing gear
S_nlg = d_nose * w_nose                                 # reference area nose landing gear
Sa_main = 0.8 * S_mlg                                   # actual frontal area main landing gear

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
    if M_cruise >= 0.7:
        sweep_c4 = np.arccos(0.75*(0.935/(0.03 + M_cruise)))
    else:
        sweep_c4 = np.arccos(1)

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
    Cl_des = CL_des / (np.cos(sweep_c4)**2)
    #Cl_des = Cl_des * np.sqrt(1 - M_cruise**2)
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
    M_to = (1.2 * V_S) / np.sqrt(1.4*287*288)
    Re_to = (1.225 * 1.2 * V_S * 0.514444 * c_mac) / 1.802e-5

    print("M_sea = ", M_sea, M_cruise2)

    print("Re =", Re, Re2, Re_sea, Re_to)
    # With CL_max = 1.8
    # CLmax = 2.464
    # CLmax take-off: 2.1 , Clmax landing: 2.25
    #     # target for take-off: Delta CLmax = 0.3
    #     # target for landing: Delta CLmax = 0.45

    dCLmax_land = 0.2984
    dCLmax_to   = 0.16316

    dClmax_land = 0.9

    hinge_c     = 75 #percent
    aileron_C   = 75 #percent

    sweep_hinge = np.arctan(np.tan(sweep_c4) - 4/AR * ((hinge_c-25)/100 * (1 - taper)/(1 + taper)))
    SwfS = dCLmax_land / (0.9 * dClmax_land * np.cos(sweep_hinge))

    Df = widthf/2 * 1.25 # clearance of 1/8 * fuselage diameter completely arbitrary

    a = -2 * (c_root - c_tip)/b
    ch = 1 - (hinge_c/100)

    D = (-4 * a * Df + 2 * c_root)**2 + 4 * 2 * a * (-2 * SwfS * S)

    print("D = ", D)

    x2 = max((-(-4 * a * Df + 2 * c_root) + np.sqrt(D))/ (-4 * a), (-(-4 * a * Df + 2 * c_root) - np.sqrt(D)) / (-4 * a))
    print("x2 = ", x2)
    print("x3 = ", (x2 + Df))

    Sw_check = ((-2 * a * Df + c_root) + (-2 * a * (Df + x2) + c_root)) * x2 / 2 / S   #verified
    print("SwfS, Sw_check", SwfS, Sw_check)

    Swf_actual = (0.25 * (-2 * a * Df + c_root) + 0.25 * (-2 * a * (Df + x2) + c_root)) * x2 / 2
    print("Swf_actual=", Swf_actual)
    # NEW LIFT CURVE

    eta = 0.95

    beta = np.sqrt(1 - M_cruise ** 2)
    CLalpha = (2 * np.pi * AR) / (2 + np.sqrt(4 + (AR * beta / eta)**2 * (1 + (np.tan(sweep_c2) / beta)**2)))
    d_alpha0l_land = -15 / 180 * np.pi * SwfS * np.cos(sweep_hinge)
    d_alpha0l_to   = -10 / 180 * np.pi * SwfS * np.cos(sweep_hinge)

    print("d_alpha0l_land, d_alpha0l_to=",d_alpha0l_land * 180 / np.pi, d_alpha0l_to * 180 / np.pi)

    alpha0L = -3.936 / 180 * np.pi     # at Re = 19520133

    CL = CLalpha * (0 - alpha0L)
    print("CLalpha=", CLalpha * np.pi / 180)
    print("alpha0l_land, alpha0l_to=", (alpha0L + d_alpha0l_land) * 180 / np.pi, (alpha0L + d_alpha0l_to) * 180 / np.pi)

    alpha_trim = CL_des / CLalpha + alpha0L
    print("alpha_trim=", alpha_trim * 180 / np.pi)

    # CLmax/Clmax vs sweepLE plot: dy > 2.5, -> CLmax/Clmax = 0.82
    # Re takeoff = 15071059, Re landing = 16326981

    CLmax_wing_landing = 0.82 * 2.38
    CLmax_wing_to      = 0.82 * 2.362

    print("CLmax at landing=",CLmax_wing_landing)
    print("CLmax at takeoff=", CLmax_wing_to)

    d_alpha_CLmax = 2.6
    alpha_stall_landing = CLmax_wing_landing/(CLalpha * np.pi / 180) + alpha0L + d_alpha_CLmax
    alpha_stall_to      = CLmax_wing_to/(CLalpha * np.pi / 180) + alpha0L + d_alpha_CLmax
    print("alpha_stall_landing=", alpha_stall_landing)
    print("alpha_stall_takeoff=", alpha_stall_to)

    alpha_range = [range(-10, 17, 1), range(-10, 14, 1), range(-10, 15, 1)]
    CL_clean_list = []
    CL_landing_list = []
    CL_to_list = []

    r = 0.9      # random factor to take into account viscosity (not done by xflr5), such that slope moves towards stall point better
    # clean
    for i in alpha_range[0]:
        CL_clean = r * CLalpha * np.pi / 180 * (i - alpha0L * 180 / np.pi)
        CL_clean_list.append(CL_clean)
    # landing
    for j in alpha_range[1]:
        CL_landing = r * CLalpha * np.pi / 180 * (j - (alpha0L + d_alpha0l_land) * 180 / np.pi)
        CL_landing_list.append(CL_landing)
    # take - of
    for k in alpha_range[2]:
        CL_to = r * CLalpha * np.pi / 180 * (k - (alpha0L + d_alpha0l_to) * 180 / np.pi)
        CL_to_list.append(CL_to)

    CLmax_list = np.array([[1.9516, 2.25, 2.1], [20.32, 16.795, 18.4946]])

    wing = [sweep_c4, sweep_c2, sweep_cLE, taper, c_root, c_tip, c_mac, y_mac, t_c, dihedral,
            Cl_des, dCLmax_land, dCLmax_to, SwfS, Re, CL_des]

    print("wing =", wing)
    print("Sweep =", sweep_cLE * 180 / np.pi)

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

    return wing, geom,cross1, hld, ail, x2, CL_clean_list, CL_landing_list, CL_to_list, alpha_range, CLmax_list


def drag():
    sweep_c4 = wing[0]
    Re = wing[14]
    SwfS = wing[13]
    c_MAC = wing[6]
    taper = wing[3]
    CL_des = wing [15]

    ####################### Zero lift drag estimation

    # wetted area, verified with https://onlinelibrary.wiley.com/doi/pdf/10.1002/9781118568101.app1

    S_wet_wing = 1.07 * 2 * S
    S_wet_tailh = 1.05 * 2 * Sh
    S_wet_tailv = 1.05 * 2 * Sv
    S_wet_fus = (np.pi * Df / 4) * ( 1/(3*L1**2) * ((4 * L1**2 + Df**2/4)**1.5 - Df**3/8) - Df + 4*L2 + 2*np.sqrt(L3**2 + Df**2/4))
    S_wet_nacelle = (25 * (1 + B)**0.2 * (Tto/ 100e3)**0.8) * 2

    ####### skin friction coeff

    # wing
    Re_wing = 44.62 * (c_MAC/k)**1.053 * M_cruise*1.16

    Cf_lam_wing = 1.328 / np.sqrt(Re_wing)
    Cf_tur_wing = 0.455 / ((np.log10(Re_wing) ** 2.58) * (1 + 0.144 * M_cruise ** 2) ** 0.65)

    Cftot_wing = 0.35 * Cf_lam_wing + 0.65 * Cf_tur_wing  # values for smooth metal

    # v tail
    Re_vtail = 44.62 * (c_MACv/k)**1.053 * M_cruise*1.16

    Cf_lam_vtail = 1.328 / np.sqrt(Re_vtail)
    Cf_tur_vtail = 0.455 / ((np.log10(Re_vtail) ** 2.58) * (1 + 0.144 * M_cruise ** 2) ** 0.65)

    Cftot_tailv = 0.35 * Cf_lam_vtail + 0.65 * Cf_tur_vtail  # values for smooth metal

    # h tail
    Re_htail = 44.62 * (c_MACh/k)**1.053 * M_cruise*1.16

    Cf_lam_htail = 1.328/np.sqrt(Re_htail)
    Cf_tur_htail = 0.455/((np.log10(Re_htail)**2.58) * (1 + 0.144 * M_cruise**2)**0.65)

    Cftot_tailh = 0.35 * Cf_lam_htail + 0.65 * Cf_tur_htail  # values for smooth metal

    # fus
    Re_fus = 44.62 * ((L1 + L2 + L3) / k) ** 1.053 * M_cruise * 1.16

    Cf_lam_fus = 1.328 / np.sqrt(Re_fus)
    Cf_tur_fus = 0.455 / ((np.log10(Re_fus) ** 2.58) * (1 + 0.144 * M_cruise ** 2) ** 0.65)

    Cftot_fus  = 0.1 * Cf_lam_fus + 0.9 * Cf_tur_fus                # values for smooth metal

    # Nacelle
    Re_nacelle = 44.62 * ((l_nacelle) / k) ** 1.053 * M_cruise * 1.16

    Cf_lam_nacelle = 1.328 / np.sqrt(Re_nacelle)
    Cf_tur_nacelle = 0.455 / ((np.log10(Re_nacelle) ** 2.58) * (1 + 0.144 * M_cruise ** 2) ** 0.65)

    Cftot_nacelle = 0.1 * Cf_lam_nacelle + 0.9 * Cf_tur_nacelle  # values for smooth metal

    print("Re_specific=", Re_wing, Re_vtail, Re_htail, Re_fus, Re_nacelle)
    print("Cf_lam_sp=", Cf_lam_wing, Cf_lam_vtail, Cf_lam_htail, Cf_lam_fus, Cf_lam_nacelle)
    print("Cf_tur_sp=", Cf_tur_wing, Cf_tur_vtail, Cf_tur_htail, Cf_tur_fus, Cf_tur_nacelle)
    print("Cf = ", Cftot_wing, Cftot_tailv, Cftot_tailh, Cftot_fus, Cftot_nacelle)
    ## Form Factor

    sweep_m_wing = np.arctan(np.tan(sweep_c4) - 4 / AR * ((x_cm_wing * 100 - 25) / 100 * (1 - taper) / (1 + taper)))
    sweep_m_tailv = np.arctan(np.tan(sweep_c4v) - 4 / AR * ((x_cm_tailv * 100 - 25) / 100 * (1 - taperv) / (1 + taperv)))
    sweep_m_tailh = np.arctan(np.tan(sweep_c4h) - 4 / AR * ((x_cm_tailh * 100 - 25) / 100 * (1 - taperh) / (1 + taperh)))

    FF_wing = (1 + 0.6 / x_cm_wing * t_c_wing + 100 * t_c_wing ** 4) * (
                1.34 * M_cruise ** 0.18 * (np.cos(sweep_m_wing)) ** 0.28)
    FF_tailh = (1 + 0.6 / x_cm_tailh * t_c_tailh + 100 * t_c_tailh ** 4) * (
                1.34 * M_cruise ** 0.18 * (np.cos(sweep_m_tailh)) ** 0.28)
    FF_tailv = (1 + 0.6 / x_cm_tailv * t_c_tailv + 100 * t_c_tailv ** 4) * (
                1.34 * M_cruise ** 0.18 * (np.cos(sweep_m_tailv)) ** 0.28)

    f_fus = lf / np.sqrt(4 * Amax_fus / np.pi)
    FF_fus = (1 + 60 / (f_fus ** 3) + f_fus / 400)
    f_nacelle = l_nacelle / np.sqrt(4 * Amax_nacelle / np.pi)
    FF_nacelle = 1 + 0.35 / f_nacelle

    ######## Miscellaneous drag
    # Wave drag
    Mdd = 0.935/np.cos(sweep_c4) - 0.14 /(np.cos(sweep_c4)**2) - CL_des/ (10*(np.cos(sweep_c4)**3))
    print("Mdd=", Mdd)
    if Mdd > M_cruise:
        wavedrag = 0.002 * (1 + 2.5 * (Mdd - M_cruise)/0.05)**(-1)
    else:
        wavedrag = 0.002 * (1 + 2.5 * (M_cruise - Mdd)/0.05)**(2.5)
    # Fuselage base drag  -> look this up
    drag_fusbase = (0.139 + 0.419 * (M_cruise - 0.161)**2) * Amax_fus

    # Drag due to fuselage upsweep (upsweep in rad, Amax is max cross-sectional area)

    dragupsweep = 3.83 * upsweep**2.5 * Amax_fus           # (not delta_cd but D/q  ??)

    # landing gear drag (add this from ADSEE)
                                            #TBD
    S_mlg = d_main * w_main

    S_nlg = d_nose * w_nose

    cds_main = main_amount * 0.04955 * np.exp(5.615 * Sa_main / S_mlg)
    drag_lg = (cds_nose + cds_main) * (S_nlg + main_amount * S_mlg) / S

    # flap drag

    if dflap > 10:
        drag_flap = 0.0144 * SwfS * (dflap - 10)
    else:
        drag_flap = 0

    drag_misc = wavedrag + drag_fusbase + dragupsweep + drag_flap + drag_lg
    leakage   = 1.05                                     # 2-5 % of total CDO



    ############ FINAL ZERO LIFT DRAG

    CD0 = 1 / S * ((S_wet_wing * Cftot_wing * IF_wing * FF_wing)
                    + (S_wet_tailh * Cftot_tailh * IF_tailh * FF_tailh)
                    + (S_wet_tailv * Cftot_tailv * IF_tailv * FF_tailv)
                    + (S_wet_fus * Cftot_fus * IF_fus * FF_fus)
                    + (S_wet_nacelle * Cftot_nacelle * IF_nacelle * FF_nacelle)
                    + drag_misc) * leakage

    print("CD0=",CD0)
    ####################### Lift induced drag
    df1 = 0      # flap deflection - clean
    df2 = 1.047  # flap deflection - Lnd
    df3 = 0.349  # flap deflection - TO

    oswaldclean = 1.78 * (1 - 0.045 * AR**0.68) - 0.64 + 0.0046 * df1
    oswaldTO = 1.78 * (1 - 0.045 * AR**0.68) - 0.64 + 0.0046 * df3
    oswaldLnd = 1.78 * (1 - 0.045 * AR**0.68) - 0.64 + 0.0046 * df2

    d_CD_twist = 0          #0.00004 * (phi_tip - phi_MGC) #effect of twist

    dAR = 0                 #effect of wing tips
    AR_eff = AR + dAR
    print("AReff=", AR_eff)

#    K_ground = (33 * (h/b)**1.5)/ (1 + 33 * (h/b)**1.5) #  ground effect
#
#    ######################### Total drag polar #######################
    clrange = np.arange(-1.0, 2.5, 0.1)
    Draglist = []

    for i in clrange:
        C_D = CD0 + (i - CL_des)**2/(np.pi*AR_eff*oswaldclean)
        Draglist.append(C_D)

    return oswaldclean, oswaldTO, oswaldLnd, Draglist, clrange

wing, geom, cross1, hld, ail, x2, CL_clean_list, CL_landing_list, CL_to_list, alpha_range, CLmax_list = wing_geometry(M_cruise, S, AR, MTOW, V_C, widthf, V_S, v_approach, V_C_TAS)

osclean,osTO,osLnd, Draglist, clrange = drag()

#----------------------------- .txt File Airfoil Coordinates

#Read from file

#Create empty lists

# lines  = [[],[],[],[]]
# xcoord1= [[],[],[],[]]
# xcoord2= [[],[],[],[]]
# ycoord1= [[],[],[],[]]
# ycoord2= [[],[],[],[]]
# camline= [[],[],[],[]]
#
# #Read from file
# f0=open('airfoil1.txt','r')
# f1=open('airfoil2.txt','r')
# f2=open('airfoil3.txt','r')
# f3=open('airfoil4.txt','r')
#
# lines0=f0.readlines()
# lines1=f1.readlines()
# lines2=f2.readlines()
# lines3=f3.readlines()
#
# for i in range(0,103):
#     xcoord1[0].append(float(lines0[i].split()[0]))
#     xcoord1[1].append(float(lines1[i].split()[0]))
#     xcoord1[2].append(float(lines2[i].split()[0]))
#     ycoord1[0].append(float(lines0[i].split()[1]))
#     ycoord1[1].append(float(lines1[i].split()[1]))
#     ycoord1[2].append(float(lines2[i].split()[1]))
# for i in range(0,70):
#     xcoord1[3].append(float(lines3[i].split()[0]))
#     ycoord1[3].append(float(lines3[i].split()[1]))
#
# for i in range(103,205):
#     xcoord2[0].append(float(lines0[i].split()[0]))
#     xcoord2[1].append(float(lines1[i].split()[0]))
#     xcoord2[2].append(float(lines2[i].split()[0]))
#     ycoord2[0].append(float(lines0[i].split()[1]))
#     ycoord2[1].append(float(lines1[i].split()[1]))
#     ycoord2[2].append(float(lines2[i].split()[1]))
# for i in range(70,139):
#     xcoord2[3].append(float(lines3[i].split()[0]))
#     ycoord2[3].append(float(lines3[i].split()[1]))
#
# for i in range(0,4):
#     xcoord2[i].insert(-1,1.0)
#     ycoord2[i].insert(-1,0.0)
#     xcoord1[i]=xcoord1[i][::-1]
#     ycoord1[i]=ycoord1[i][::-1]
#
# ##Camber Line
# for i in range(0,len(xcoord1[0])):
#     camline[0].append((ycoord1[0][i]+ycoord2[0][i])/2)
#     camline[1].append((ycoord1[1][i]+ycoord2[1][i])/2)
#     camline[2].append((ycoord1[1][i]+ycoord2[1][i])/2)
# for i in range(0,len(xcoord1[3])):
#     camline[3].append((ycoord1[3][i]+ycoord2[3][i])/2)
#
#
# #----------------------------- Plotting
#
plt.figure(0)
plt.plot(geom[0], geom[1], geom[2], geom[3], hld[0], hld[1],ail[0],ail[1])
plt.text(cross1[0],cross1[1],'Fuselage Wall Line')
plt.grid(True,which="major",color="#999999")
plt.grid(True,which="minor",color="#DDDDDD",ls="--")
plt.minorticks_on()
plt.ylim(-12.0,2.0)
plt.ylabel('x [m]')
plt.xlabel('y [m]')
#
# plt.figure(1)
# plt.grid(True,which="major",color="#999999")
# plt.grid(True,which="minor",color="#DDDDDD",ls="--")
# plt.minorticks_on()
# plt.plot(xcoord1[0],ycoord1[0],color='r')
# plt.plot(xcoord2[0],camline[0],'--',color='r')
# plt.plot(xcoord2[0],ycoord2[0],color='r')
# plt.xlim(0,1)
# plt.ylim(-0.4,0.4)
# plt.text(0.0,0.0,'LE')
# plt.text(1.0,0.0,'TE')
# plt.ylabel('y/c [-]')
# plt.xlabel('x/c [-]')

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
#plt.text(0.0,0.0,'LE')git pul
#plt.text(1.0,0.0,'TE')
#plt.ylabel('y/c [-]')
#plt.xlabel('x/c [-]')
#



plt.figure(5)
plt.grid("minor")
plt.plot(alpha_range[0], CL_clean_list)
plt.plot(alpha_range[2], CL_to_list)
plt.plot(alpha_range[1], CL_landing_list)
plt.plot(CLmax_list[1][0], CLmax_list[0][0], marker=".", color="blue")
plt.plot(CLmax_list[1][2], CLmax_list[0][2], marker=".", color="orange")
plt.plot(CLmax_list[1][1], CLmax_list[0][1], marker=".", color="green")
plt.ylim(0, 2.5)
plt.xlim(-5, 25)
plt.xlabel(r"$\alpha$ [deg]")
plt.ylabel(r"$C_L$ [-]")
plt.legend(["clean", "take-off", "landing",  "stall clean", "stall take-off", "stall landing"], loc="upper left")

plt.figure(6)
plt.grid()
plt.plot(clrange, Draglist)
plt.show()





