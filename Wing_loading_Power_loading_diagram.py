# -*- coding: utf-8 -*-
"""
Created on Mon May 11 14:35:32 2020

@author: Gebruiker
"""

import matplotlib.pyplot as plt
import numpy as np

# Constants
rho_0 = 1.225  # [kg/m^3]
x_axis = np.arange(1, 6000)  # upper limit to be adjusted
s_land = 1500  # [m]

# Colours
blue = ['aqua', 'darkturquoise', 'steelblue', 'mediumblue', 'darkblue']
green = ['greenyellow', 'lawngreen', 'limegreen', 'forestgreen', 'darkgreen']
black = ['lightgrey', 'silver', 'darkgrey', 'dimgray', 'black']
red = ['lightsalmon', 'coral', 'orangered', 'firebrick', 'red']
purple = ['violet', 'palevioletred', 'deeppink', 'mediumvioletred', 'purple']
orange = ['bisque', 'orange', 'darkorange', 'goldenrod', 'darkgoldenrod']
dark = ['mediumslateblue', 'mediumpurple', 'rebeccapurple', 'blueviolet']
j = 0

# Sizing for stall speed (vertical line for clean configuration stall and landing configuration stall)
CLmax_land = [1.9, 2.2, 2.25, 2.8]  # Assume for each concept [1.9-3.3] conventional
CLmax_clean = [1.5, 1.6, 1.7, 1.8, 1.9]  # Assume for each concept [1.5-1.9] conventional
Vstall = 63  # 56      #[m/s]    #Assume for each concept, based on an assumed approach speed / 1.3


def stallspeedland(CLmax_land):
    return 0.5 * rho_0 * Vstall ** 2 * CLmax_land


def stallspeedclean(CLmax_clean):
    return 0.5 * rho_0 * Vstall ** 2 * CLmax_clean


# Sizing for take off (assume range of CL_max at take-off)
TOP_jet = 125 * 47.89  # Assume if concept has jet engines, based on graph on slide 100 and multiply with 47.89
TOP_prop = 650 * 3.53  # Assume if concept has turboprop engines, based on graph on slide 100 and multiply with 3.53
eta_prop = 0.95  # Assume a propeller effieciency in case your concept has a turboprop engine
rho_h = 1.225  # 0.974                                   #Specify density at take-off altitude and landing   (check for 0.974), take upper limit of regional airport the aircraft will operate from
sigma = rho_h / rho_0


def takeoff_jet(x, CL_to, sigma):
    return x * (1 / CL_to) * (1 / TOP_jet) * (1 / sigma)


def takeoff_prop(x, CL_to, sigma):
    return (CL_to * TOP_prop * sigma) / (x * eta_prop)


CLmax_to = [1.7, 1.9, 2.0, 2.1, 2.2]  # Assume for each conept conventional [1.7-2.1]
for i in CLmax_to:
    CL_to = i / 1.1 ** 2
    plt.plot(x_axis, takeoff_jet(x_axis, CL_to, sigma), label="CLmaxto = " + str(i),
             color=str(blue[j]))  # Disable if your conept has turboprop engines
    # plt.plot(x_axis, takeoff_prop(x_axis, CL_to,sigma), label="takeoff, CLmaxto = "+str(i), color=str(blue[j]))       #Disable if your conept has jet engines
    j = j + 1
j = 0
# Sizing for landing
f = 0.8058  # Calculate based on class 1 weight estimation for your concept, f = landing weight / take-off weight


def landing(CLmax_land):
    return (CLmax_land * rho_h * s_land) / (2 * 0.5847 * f)


for i in CLmax_land:
    plt.axvline(x=landing(i), label="Land, CLmaxland =" + str(i), color=str(green[j]))
    plt.axvline(x=stallspeedland(i), label="Vs, CLmaxland =" + str(i), color=str(dark[j]))
    j = j + 1
j = 0
for i in CLmax_clean:
    plt.axvline(x=stallspeedclean(i), label="Vs, CLmaxclean =" + str(i), color=str(orange[j]))
    j = j + 1

# Sizing for cruise speed
rho_c = 0.41268  # Determine cruise density based on cruise altitude of your concept, 0.41268 at 10 km altitude
V_c = 230  # [m/s]          #Determine cruise speed for your concept, 230 m/s like CRJ700
CD0 = 0.01277  # Assume CDO in clean configuration for your concept, cdo = 0.0145 for jet aircraft?
A = [7, 8, 9, 10]  # Assume range of feasible aspect ratios for your concept
e = 0.85  # Assume clean configuration oswald factor for your conept
PS_c = 0.9  # Assume the power setting or thrust level at cruise condition, in this case 90%, should be between 0 and 1
f_c = 0.8587  # Mass fraction at cruise condition, cruise weight / MTOW, determine for your concept for start of cruise


def cruisespeed_jet(x, A):  # too high values
    return (f_c / PS_c * (rho_0 / rho_c) ** (3 / 4)) * (
                (CD0 * 0.5 * rho_c * V_c ** 2) / (f_c * x) + x * f_c / (np.pi * A * e * 0.5 * rho_c * V_c ** 2))


def cruisespeed_prop(x, A):
    return (PS_c / f_c * eta_prop * (rho_0 / rho_c) ** (3 / 4) * (
    ((CD0 * 0.5 * rho_c * V_c ** 3) / (f_c * x) + f_c * x / (np.pi * A * e * 0.5 * rho_c * V_c)))) ** (-1)


# sizing for climb rate performance

##delta_takeoff_flaps_cdo = [.01-.02],  delta_e=0.05,
## delta_undercarieage_cdo=[.015-.025], delta_e=0
## delta_landing_flaps_cdo =[.055-.075], delta_e=0.10


CL_c = 0.4  # Assume cruise CL for your concept

CD0_toGU = CD0 + .015  # Assume CD0 at take off configuration with gear up, look at Adsee slides for conversion factors for your conept
CD0_toGD = CD0 + .015 + .02  # Assume CD0 at take off configuration with gear down for your conept
CD0_landGD = CD0 + .02 + .065  # Assume CD0 at landing configuration with gear down for your conept

e_toGD = e + 0.05 + 0  # Assume e at take off configuration with gear down for your conept
e_toGU = e + 0.05  # Assume e at take off configuration with gear up for your conept
e_landGD = e + 0 + 0.1

# Assume e at landing configuration with gear down for your conept


c = 10  # [m/s]               #Assume ROC for your concept, maybe this is related to requirements and set equal for all concepts


def climbrate_jet(x, A, CL_roc):
    return c / (np.sqrt(x) * np.sqrt((2 / rho_0) * (1 / CL_roc))) + 4 * CD0_toGU / CL_roc


def climbrate_prop(x, A, CD_toGU):
    return eta_prop / (c + np.sqrt(2 * x / rho_0) / (1.345 * (A * e) ** (3 / 4) / (CD_toGU ** (1 / 4))))


# Sizing for climb gradient performance
c_over_V = 0.024  # Assume climb gradient for your concept, maybe this is related to requirements and set equal for all concepts
n = 2  # Specify the number of engines your concept has
n_max = 3
v_m = 110  # m/s manoeuvring velocity


def climbgradient_prop(x, CL_grad, CD_toGU):
    return eta_prop / (np.sqrt(x) * (c_over_V + CD_toGU / CL_grad) * np.sqrt(2 / (rho_h * CL_grad)))


def climbgradient_jet(A, CD_toGU):
    return n / (n - 1) * (c_over_V + 2 * np.sqrt(CD_toGU / (np.pi * A * e_toGU)))


def manoeuvring_jet(v_m, n_max, A, e, x):
    return CD0 * 0.5 * rho_0 * v_m ** 2 / x + x * n_max ** 2 / (np.pi * A * e * 0.5 * rho_h * v_m ** 2)


def manoeuvring_prop(v_m, n_max, A, e, x, eta_prop):
    return ((CD0 * 0.5 * rho_0 * v_m ** 3 / x + x * n_max ** 2 / (np.pi * A * e * 0.5 * rho_h * v_m)) / eta_prop) ** (
        -1)


j = 0
for i in A:
    CL_grad = np.sqrt(CD0_toGU * np.pi * i * e_toGU)
    CL_roc = np.sqrt(3 * CD0_toGU * np.pi * i * e_toGU)
    CD_toGU = CD0_toGU + (CLmax_to[1] / 1.21) ** 2 / (
                i * np.pi * e_toGU)  # CHOOSE ONE CL_to value from the array put it in CL_to[]
    plt.plot(x_axis, cruisespeed_jet(x_axis, i), label="Vc, A =" + str(i),
             color=str(green[j]))  # Disable if your concept has turboprop engines
    # plt.plot(x_axis,cruisespeed_prop(x_axis,i),label="cruisespeed, A ="+str(i), color=str(green[j]))        #Disable if your concept has jet engines
    plt.plot(x_axis, climbrate_jet(x_axis, i, CL_roc), label="c, A=" + str(i),
             color=str(black[j]))  # Disable if your concept has turboprop engines
    # plt.plot(x_axis,climbrate_prop(x_axis,i, CD_toGU),label="climbrate, A="+str(i),color=str(black[j]))										#Disable if your concept hasjet engines
    plt.axhline(y=climbgradient_jet(i, CD_toGU), label="c/V, A=" + str(i),
                color=str(red[j]))  # Disable if your concept has turboprop engines
    # plt.plot(x_axis, climbgradient_prop(x_axis,CL_grad, CD_toGU),label="climbgradient, A="+str(i),color=str(red[j]))             			#Disable if your concept has jet engines
    plt.plot(x_axis, manoeuvring_jet(v_m, n_max, i, e, x_axis), label="manoeuvre, A=" + str(i), color=str(purple[j]))
    # plt.plot(x_axis, manoeuvring_prop(v_m,n_max,i,e,x_axis,eta_prop),label="manoeuvring performance, A="+str(i),color=str(purple[j]))
    j = j + 1

plt.xlim(0, 6000)
plt.ylim(0, 1)
plt.xlabel('W/S [N/m$^2$]')
plt.ylabel('T/W [-]')
# plt.legend(loc='center right', bbox_to_anchor=(1, 0.5))
plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05),
           fancybox=True, shadow=True, ncol=8)
plt.show()

# other
Vto = 1.1 * Vstall
Vsland = np.sqrt(s_land / 0.584)
Vapproach = 1.3 * Vsland