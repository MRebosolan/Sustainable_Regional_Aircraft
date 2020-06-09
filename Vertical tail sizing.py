import numpy as np
import Class_2_estimation as Cl2
import input
import scissor_plot_wing_shift as sc_shift
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
S =   Cl2.S                # input.S  [m]
b =   Cl2.b                # input.b  [m]
T_OEI = 0.5*Cl2.Tto        # Thrust at 1 engine inoperative, 
rho_c = input.rho_c
g = input.g
MTOW = g*Cl2.MTOM
CLmax = input.CLmax_land 
y_engine = 0.35*Cl2.b/2    # y distance of the engine from the c.g.
rho0 = input.rho0
rho_to = rho0
Vbar_v = 0.086             # Vertical tail volume
Sr_over_Sv = 0.303         # Ratio of rudder area over vtail area
Clv_max_r = 1.2            # Maximum lift coefficient vtail at max rudder deflection
Vmin = np.sqrt(MTOW * 2 /(S * rho_to * CLmax))                  #or change to Vs, depends on whether the wing is sized for Vs or for CLmax
Vlof = 1.05 * Vmin
V_c = 0.514444*input.V_C
x_v =  sc_shift.x_ac_v_nose-sc_shift.cg_loaded_nose   # input.x_v  #Distance Aerodynamic centre Vtail to c.g.
vtail_sweep = input.half_chord_sweep_vert # Radians
AR_vtail = input.AR_v
taper_v = input.taper_v


import numpy as np
import Class_2_estimation as Cl2
import input

rho_c = input.rho_c
g = input.g
MTOW = g*Cl2.MTOM
CLmax = input.CLmax_land 
S = Cl2.S
b = Cl2.b
x_v =  10.5                 #input.x_v  #Distance Aerodynamic centre Vtail to c.g.
T_OEI = 0.5*Cl2.Tto                #x_v can be calculated the location of the htail aerodynamic center is know wrt either the aircraft c.g. OR the wing AC, provided the distance of the wing AC is know wrt the nose
y_engine = 0.35*b/2     #y distance of the engine from the c.g.
rho_to = input.rho0
Vbar_v = 0.086              #Vertical tail volume
Sr_over_Sv = 0.303          #Ratio of rudder area over vtail area
Clv_max_r = 1.2             # Maximum lift coefficient vtail at max rudder deflection
Vmin = np.sqrt(MTOW * 2 /(S * rho_to * CLmax))      #or change to Vs, depends on whether the wing is sized for Vs or for CLmax
Vlof = 1.05 * Vmin
V_c = 0.514444*input.V_C    #
vtail_sweep = input.half_chord_sweep_vert # Radians
AR_vtail = input.AR_v
taper_v = input.taper_v

print('make sure to check that the critical mach number of the vtail and htail are higher than that of the wing! use Roskam book II, page 150')
print ('Furthermore, check whether the take-off thrust is properly linked')
#Calculatetes preliminary area of the vertical tail surface and rudder
def S_v(Vbar_v,S,b,x_v):
    S_v = Vbar_v*S*b/x_v
    S_r = Sr_over_Sv*S_v
    return S_v,S_r
S_v,S_r = S_v(Vbar_v,S,b,x_v)

def rudder_sizing(T_OEI=T_OEI,Clv_max_r=Clv_max_r,y_engine=y_engine,rho_to=rho_to,rho_c=rho_c,Vlof=Vlof,V_c=V_c,vtail_sweep=vtail_sweep,x_v=x_v):
    T_moment_OEI = T_OEI*y_engine
    Req_Sv = []
    Rho = [rho_to,rho_c]
    V = [Vlof/np.cos(vtail_sweep),V_c/np.cos(vtail_sweep)]
    for i in range(len(Rho)):
        Req_Sv1 = 1.25*T_moment_OEI/(0.5*Rho[i]*V[i]**2*Clv_max_r*x_v) # maximum rudder deflection 25 degrees
        Req_Sv.append(Req_Sv1)
    return Req_Sv
Req_Sv = rudder_sizing()
print ()
print ('The vertical tail area (Class I) estimate is', np.round(S_v,4), 'm^2. The rudder area equals' , np.round(S_r,3),'m^2.')
print ()

def Vtail_dim_class1(taper_v,AR_vtail,vtail_sweep,S_v):
    b_v_c1 = np.sqrt(AR_vtail*S_v)
    chord_r_c1 = 2*S_v/((1+taper_v)*b_v_c1)
    chord_t_c1 = taper_v*chord_r_c1
    MAC_v_cl1 = chord_r_c1*2/3*((1+taper_v+taper_v**2)/(1+taper_v))
    return b_v_c1,chord_r_c1,chord_t_c1,MAC_v_cl1
b_v_c1,chord_r_c1,chord_t_c1,MAC_v_cl1 = Vtail_dim_class1(taper_v,AR_vtail,vtail_sweep,S_v)

print ()
print ()
print ('Root chord Class I:',np.round(chord_r_c1,4))
print ('Tip chord Class I:',np.round(chord_t_c1,4))
print ('Span Class I:',np.round(b_v_c1,4))
print ('MAC Class I:',np.round(MAC_v_cl1,4))
print ()
print ()
print ('The required vtail area to counteract the OEI moment at take-off is',Req_Sv[0],'[m^2]')
print ('The required vtail area to counteract the OEI moment at cruise is',Req_Sv[1],'[m]')
print ('The largest required area is',np.round(np.max(Req_Sv),4),'at take-off. The required rudder area is',np.round(np.max(Req_Sv)*Sr_over_Sv,4))

"""
def Rudder_def(T_OEI,y_engine,S,b,rho_to):
    qbar_mc = 1/2*rho_to*S*(Vmin*1.2)**2
    N_t_crit = T_OEI*y_engine
    N_D =0.25        # jet engine with high bypass ratio, 0.15 for low
    Cn_delta_r = -Cy_delta_r ()    # Roskam book VI
    r_def_max = ((N_D+1)*N_t_crit)/(qbar_mc*S*b*Cn_delta_r)
    if  np.degrees(r_def_max) > 25:
        print ('error, too large rudder deflection required. Adjust vtail size')
    else:
        print ('the rudder deflection is not too large, the vtail is capable of counteracting the yaw moment')
    return
"""