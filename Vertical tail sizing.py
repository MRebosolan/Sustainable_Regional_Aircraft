import numpy as np
import Class_2_estimation as Cl2
import input

"""
Responsible person: Tobias  | For now: Start from line 47

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



S =   Cl2.S     #input.S  [m]
b =   Cl2.b     #input.b  [m]
x_v =  10.5     #input.x_v  #Distance Aerodynamic centre Vtail to c.g.
                #x_v can be calculated the location of the htail aerodynamic center is know wrt either the aircraft c.g. OR the wing AC, provided the distance of the wing AC is know wrt the nose
T_OEI = 0.5*Cl2.Tto   # Thrust at 1 engine inoperative, 
rho_c = input.rho_c
g = input.g
MTOW = g*Cl2.MTOM
CLmax = input.CLmax_land 
y_engine = 0.35*Cl2.b/2      #y distance of the engine from the c.g.
rho_to = input.rho0
Vbar_v = 0.086  #Vertical tail volume
Sr_over_Sv = 0.303    #Ratio of rudder area over vtail area
Clv_max_r = 1.2  # Maximum lift coefficient vtail at max rudder deflection
Vmin = np.sqrt(MTOW * 2 /(S * rho_to * CLmax))                  #or change to Vs, depends on whether the wing is sized for Vs or for CLmax
Vlof = 1.05 * Vmin
V_c = 0.514444*input.V_C

import numpy as np
import Class_2_estimation as Cl2
import input

rho_c = input.rho_c
g = input.g
MTOW = g*Cl2.MTOM
CLmax = input.CLmax_land 
S = 80
b = 25
x_v =  10.5     #input.x_v  #Distance Aerodynamic centre Vtail to c.g.
T_OEI = 80E3                #x_v can be calculated the location of the htail aerodynamic center is know wrt either the aircraft c.g. OR the wing AC, provided the distance of the wing AC is know wrt the nose
y_engine = 0.35*Cl2.b/2       #y distance of the engine from the c.g.
rho_to = input.rho0
Vbar_v = 0.086  #Vertical tail volume
Sr_over_Sv = 0.303    #Ratio of rudder area over vtail area
Clv_max_r = 1.2  # Maximum lift coefficient vtail at max rudder deflection
Vmin = np.sqrt(MTOW * 2 /(S * rho_to * CLmax))                  #or change to Vs, depends on whether the wing is sized for Vs or for CLmax
Vlof = 1.05 * Vmin
V_c = 0.514444*input.V_C
vtail_sweep = np.radians(30) # Radians

#Calculatetes preliminary area of the vertical tail surface and rudder
def S_v(Vbar_v,S,b,x_v):
    S_v = Vbar_v*S*b/x_v
    S_r = Sr_over_Sv*S_v
    return S_v,S_r
S_v,S_r = S_v(Vbar_v,S,b,x_v)

def rudder_sizing(T_OEI=T_OEI,Clv_max_r=Clv_max_r,y_engine=y_engine,rho_to=rho_to,rho_c=rho_c,Vlof=Vlof,V_c=V_c,vtail_sweep=vtail_sweep):
    T_moment_OEI = T_OEI*y_engine
    Req_Sv = []
    Req_y_engine = []
    Rho = [rho_to,rho_c]
    V = [Vlof/np.cos(vtail_sweep),V_c/np.cos(vtail_sweep)]
    for i in range(len(Rho)):
        Req_Sv1 = T_moment_OEI/(0.5*Rho[i]*V[i]**2*Clv_max_r*y_engine)
        Req_y_engine1 = T_moment_OEI/(0.5*Rho[i]*V[i]**2*Clv_max_r*S_v)
        Req_Sv.append(Req_Sv1)
        Req_y_engine.append(Req_y_engine1)
    return Req_Sv, Req_y_engine
Req_Sv,Req_y_engine = rudder_sizing()

print ('The vertical tail area (Class I) estimate is', np.round(S_v,4), 'm^2. The rudder area equals' , np.round(S_r,3),'m^2.')
print ()
print ('The required vtail area to counteract the OEI moment at take-off is',Req_Sv[0],'[m^2]')
print ('The The maximum allowable engine distance from the c.g. at take-off is',Req_y_engine[0],'[m]')
print ()
print ('The required vtail area to counteract the OEI moment at cruise is',Req_Sv[1],'[m]')
print ('The The maximum allowable engine distance from the c.g. at take-off is',Req_y_engine[1],'[m]')

# Research engine placement regulations (Stability related)
# Slide 447-449 ADSEE I 0.35b/2