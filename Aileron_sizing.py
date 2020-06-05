# -*- coding: utf-8 -*-
"""
Created on Fri May 29 09:44:28 2020

@author: Rick van Rooijen
This file can be used to size the ailerons. By varying the input parameters that are not fixed
i.e. varying the aileron geometry, you can size the aileron until the roll rate requirement is met

input: see input parameters below under the tab "aileron sizing parameters, not to be changed"
output: aileron geometry (spanwise position, roll rate achieved, rolldamping coefficient, dcl/d(delta a), aileron chord/wing chord
"""
import numpy as np
import input
import Class_2_estimation as Cl2

#Aileron sizing input parameters, not to be changed
Preq = np.radians(45/1.45)              #required roll rate from Class II aircraft regulations
Cr = input.Cr                                  #Wing root chord
Ct = input.Ct                               #Wing tip chord
Dfus = input.widthf                                #Fuselage diameter
S = Cl2.S                                 #Wing area
b = Cl2.b                                  #Wing span
Cla = input.Cla_aileron                     #1/rad, sectional lift curve slope at wing section where aileron is located, determine by datcom method or airfoil simulation
Cd0 = input.Cd0_aileron                             #zero drag coefficient at wing section where aileron is located, determine by airfoil simulation

#Aileron sizing input parameters: adjust the geometry of the aileron untill desired roll rate is achieved
tau = 0.5                               #Aileron effectiveness based on control-surfce-to-lifting-surface ratio. Read off graph
delta_a_up_max = np.radians(25)             #Change to obtain required aileron geometry, compare to reference aircraft for realistic values
b1 = 8                                      #aileron inside y position
b2 = 11                                     #aileron outside y position
V = 60                                      #Speed at which roll rate is computed, stall speed is most sizing propably
#Minor calculations with input parameters
Sref = S - Cr * Dfus                    #Wing reference area
delta_a_do_max = delta_a_up_max * 0.75      #In order to prevent adverse yaw, upgoing deflection should be 4/3 times larger than downward max deflection
delta_a = (delta_a_up_max + delta_a_do_max) / 2    #average maximum aileron deflection

def roll_rate(Clda, Clp, delta_a, V, b):
    P = - Clda / Clp * delta_a * 2 * V / b
    return P

def Clda(Cr, Ct, b, Cla, tau, Sref, b1, b2):
    Clda = 2 * Cla * tau / (Sref * b) * (0.5 * Cr * (b2**2 - b1**2) + 2 * (Cr - Ct) / (3 * b) * (b1**3 - b2**3))
    return Clda


def Clp(Cla, Cd0, Sref, b, Cr, Ct):
    Clp = -4 * (Cla + Cd0) / (Sref * b**2) * (Cr * b**3 / 24 - b**3 * (Cr - Ct) / 32)
    return Clp

Clda = Clda(Cr, Ct, b, Cla, tau, Sref, b1, b2)
Clp = Clp(Cla, Cd0, Sref, b, Cr, Ct)
P = np.degrees(roll_rate(Clda, Clp, delta_a, V, b))
print(Clp)
print("Clda = ", Clda)
print("Clp = ", Clp)
print("roll_rate [deg/s]=", P)

if P < 0 or P < np.degrees(Preq):
    print("The aileron geometry parameters DO NOT meet the roll rate requirement")

else:
    print("The roll rate requirement is met with a difference of", P-np.degrees(Preq), "[deg/s]")


