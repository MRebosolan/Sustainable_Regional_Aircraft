import numpy as np
import Class_2_estimation as Cl2
import input

"""
Responsible person: Tobias

This code requires as inputs:
    S          | Wing area
    b          | Wing span
    x_v        | Distance of the aerodynamic center of vtail to aircraft cg
    Vbar_v     | Obtained from statistics
    Sr_over_Sv | Ratio of rudder area over vertical tail area: obtained from statistics

This code gives as outputs:
    S_v        | Vertical tail area (Roskam book II)
    S_r        | Rudder area (Roskam book II)
    

Only value that still has to be determined is x_v
x_c: Distance from AC of Vtail to c.g of the aircraft
Vbar_v (Vertical tail volume) and Sr_over_Sv (Rudder area over Vtail area 
have both been determined from statistics, the excel file can be found on 
the google drive Final report -> stability and control
"""

Vbar_v = 0.086  #Vertical tail volume
Sr_over_Sv = 0.303    #Ratio of rudder area over vtail area

S =   Cl2.S     #input.S  [m]
b =   Cl2.b     #input.b  [m]
x_v =  10.5     #input.x_v  #Distance Aerodynamic centre Vtail to c.g.
                #x_v can be calculated the location of the htail aerodynamic center is know wrt either the aircraft c.g. OR the wing AC, provided the distance of the wing AC is know wrt the nose

#Calculatetes preliminary area of the vertical tail surface and rudder
def S_v(Vbar_v,S,b,x_v):
    S_v = Vbar_v*S*b/x_v
    S_r = Sr_over_Sv*S_v
    return S_v,S_r
S_v,S_r = S_v(Vbar_v,S,b,x_v)
print ('The vertical tail area equals', np.round(S_v,4), 'm^2. The rudder area equals' , np.round(S_r,3),'m^2.')