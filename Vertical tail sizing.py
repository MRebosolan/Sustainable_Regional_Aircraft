import matplotlib.pyplot as plt
import numpy as np
import Class_2_estimation as Cl2
import input
import scissor_plot_wing_shift as sc_shift
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
x_v = 13.253442870709486    # input.x_v  #Distance Aerodynamic centre Vtail to c.g.
vtail_sweep = input.half_chord_sweep_vert # Radians
AR_vtail = 0.9 * input.AR_v
taper_v = input.taper_v
x_cg  = sc_shift.cg_loaded_nose  
x_cg_aft = sc_shift.cg_aft*sc_shift.MAC + sc_shift.xlemac       # x-location of most aft cg


#print('make sure to check that the critical mach number of the vtail and htail are higher than that of the wing! use Roskam book II, page 150')
#print ('Furthermore, check whether the take-off thrust is properly linked')
print ()
print ()
print ()
#Calculatetes preliminary area of the vertical tail surface and rudder
def S_v(Vbar_v,S,b,x_v):
    S_v = Vbar_v*S*b/x_v             # Vtail area
    S_r = Sr_over_Sv*S_v             # Rudder area
    return S_v,S_r
S_v,S_r = S_v(Vbar_v,S,b,x_v)
#print ('initial vertical tail area:',S_v,S_r)
def rudder_design(y_engine,T_OEI,S,b,vtail_sweep,taper_v,AR_vtail,x_cg,rho_to,Vmin):
    stepsize = 0.0025
    running= True
    S_vlist=[S_v]
    S_vlist1=[S_v]
    S_vfinal = [S_v]
    delta_r_calc_list = []
    Cn_delta_r_list = []
    
    while running:  
        b_v_c1 = np.sqrt(AR_vtail*S_vlist[-1])
        chord_r_c1 = 2*S_vlist[-1]/((1+taper_v)*b_v_c1)
        chord_t_c1 = taper_v*chord_r_c1
        MAC_v_cl1 = chord_r_c1*2/3*((1+taper_v+taper_v**2)/(1+taper_v))
        z_v = MAC_v_cl1/chord_r_c1*taper_v*b_v_c1 #vertical distance from AC_vtail to AC_aircraft
        lv = 23.50+.5*chord_r_c1+np.tan(vtail_sweep)*z_v- 1/4*MAC_v_cl1-x_cg_aft

        Nt_crit = y_engine*T_OEI      # Critical engine-out yawing moment
        N_D = 0.25                    # For jet driven aircraft with high by-pass ratio 0.15 for low b.p.r.
        Vmc = 1.2*Vmin                # minimum control speed
        q_mc = 0.5*rho_to*(Vmc)**2    # Dynamic pressure at Vmc
       
        #Keep in mind that we need the max. value, so alpha (which would be the beta, can be roughly between 0 and 30 degrees)
        
       # cf_over_c = 0.3              # Rudder chord/vtail chord
        k_acc = 0.67                  # S Obtained from F. 8.14 Roskam VI
        K_b = 0.97                    # [TODO] Obtained from F. 8.51 and 8.52 Roskam VI
        #a_delta_cl = 0.66            # Obtained from F. 8.53 Roskam VI
        a_delta_CL_cl=    1.10 #1.1   # S (0.55+0.78)/2 # between 0.5-0.78 Alpha_delta_CL/Alpha_delta_cl F. 8.53 Roskam VI
        cl_delta_theory = 4.42        # [1/rad] 4.85 for Joukovsky t/c .09, 5.05 for NACA 0012 Obtained from F. 8.14 Roskam VI
         
        Mach_to = Vmc/np.sqrt(288.15*1.4*287) # S 0.1987 Mach number at sealevel
        beta  = (1-(Mach_to**2))**0.5 # S Compressibility effect Prandtl
        cl_alpha_v = (0.10625*180/np.pi)  # S [1/rad] Obtained from literature NACA0012 (Martijn Source)
        Cl_cl_delta_theory = .945     # S Obtained from F 8.15 Roskam VI, with cl_alpha_v of the NACA0012 of 0.10625 [1/deg]      
    
        A_v = AR_vtail
        A_vf_over_A_v =  1.15         # [TODO] constant value Vtail AR in prsence of fuselage to isolated vtail F 10.14 Roskam VI
        #[TODO] ac_htail/b_v_c1
        A_vhf_over_Avf = 0.84 #1.10   # [TODO] ratio Vtail AR in prsence of fuselage AND htail to fuselage alone F 10.15 Roskam VI
        K_vh = 1.05 #1.04 at Sh =20,1.06 at Sh=22   # [TODO] dependent on Sh/Sv F 10.16 Roskam VI p422
        
        Av_eff = (A_vf_over_A_v)*A_v*(1+K_vh*(A_vhf_over_Avf-1))  #Effective aspect ratio
        delta_c2 = np.radians(25)     # [CHECK for correct sweep angle] semi-chord sweep angle Roskam VI p284                          
        k = cl_alpha_v/(2*np.pi/beta)
        Cl_alpha_v = 2*np.pi*Av_eff/(2+((Av_eff**2*beta**2/(k**2))*(1+(np.tan(delta_c2))**2/beta**2)+4)**0.5)

        Cy_delta_r = (Cl_alpha_v/cl_alpha_v)*(k_acc*K_b)*(a_delta_CL_cl)*(Cl_cl_delta_theory)*cl_delta_theory*(S_vlist[-1]/S)

        Cn_delta_r = []
        for alpha in np.arange(0,25+1,1): # CHECK this Roskam II
            Cn_delta_r1 = -Cy_delta_r * (lv*np.cos(np.radians(alpha))+z_v*np.sin(np.radians(alpha)))/b
            Cn_delta_r.append(Cn_delta_r1)
        #print ('-z_H/b_v = ',-z_v/b_v_c1)

        Cn_delta_r = min(Cn_delta_r,key=abs) # On point, -0.0012 [1/rad] is a good value CHECK whether you need min or max; largest req tail area is what you need
        Cn_delta_r_list.append(Cn_delta_r)

        delta_r_calc = np.degrees(((N_D+1)*Nt_crit)/(q_mc*S*b*Cn_delta_r))    #[rad] rudder deflection, 25 degrees max Roskam book VI p 494
        delta_r_calc_list.append(delta_r_calc)

        if abs(delta_r_calc) > 35:
            S_new = S_vlist[-1]+stepsize
            S_vlist.append(S_new)
            S_vlist1.append(S_new)
            S_vfinal[0] = S_new+2*stepsize #Added 1 additional setpsize for safety
            running = True
       
        elif abs(delta_r_calc) <= 25 and abs(delta_r_calc) > 15 and len(S_vlist1)<20000:
            S_new = S_vlist1[-1]+stepsize
            S_vlist1.append(S_new)
            S_vlist[-1]=S_new
            running = True
        
        else:
            S_vlist[-1] = S_vlist[-1] + stepsize
            print ('The deflection angle equals:',delta_r_calc,'[deg]')
            print ('The required surface area for the vertical tail equals:',S_vlist[-1])
            print ('The required rudder area for the vertical tail equals:',Sr_over_Sv*S_vlist[-1])
            print ('vtial ac to nose:',lv+x_cg_aft)
            print ('vtail rchord',chord_r_c1)
            print ('MAC',MAC_v_cl1)
            print ('span',b_v_c1)
            print (z_v)
            running = False
    
    return S_vlist,S_vlist1,delta_r_calc,S_vfinal,delta_r_calc_list,Cn_delta_r_list

S_vlist,S_vlist1,delta_r_calc,S_vfinal,delta_r_calc_list,Cn_delta_r_list = rudder_design(y_engine, T_OEI, S, b, vtail_sweep, taper_v, AR_vtail, x_cg, rho_to, Vmin)

fig = plt.figure(10)
plt.subplot(121)
plt.plot(S_vlist1,delta_r_calc_list)#,[S_vlist1[0],S_vlist1[-1]],[-25,-25])
plt.xlabel('$S_{v}$ [$m_{2}$]')
plt.ylabel('Rudder deflection $\delta$ [deg]')
plt.grid(b=True, which='major', color='#666666', linestyle='-')
plt.minorticks_on()
plt.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.2)


plt.subplot(122)
plt.plot(Cn_delta_r_list,delta_r_calc_list)#,[S_vlist1[0],S_vlist1[-1]],[-25,-25])
plt.xlabel('$C_{n_{\delta_{r}}}$ [-]')
plt.ylabel('Rudder deflection $\delta$ [$deg$]')
plt.grid(b=True, which='major', color='#666666', linestyle='-')
plt.minorticks_on()
plt.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.2)

fig = plt.figure(11)
plt.plot(S_vlist1,Cn_delta_r_list)#,[S_vlist1[0],S_vlist1[-1]],[-25,-25])
plt.xlabel('Vertical tail area [m2]')
plt.ylabel('Rudder deflection $\delta$ [deg]')
plt.grid(b=True, which='major', color='#666666', linestyle='-')
plt.minorticks_on()
plt.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.2)


plt.show()
print ()
print (S_vfinal[0],S_vfinal[0]*Sr_over_Sv,np.sqrt(S_vfinal[0]*1.7))
print ()
print (S_vfinal[0]/S)

#Htail area = 23.58m**2
#span htail = 9.7m
#print ('make sure that K_vh and A_vhf_over_Avf, cl_alpha_v and Cl_cl_delta_theory have the correct values; it needs to be adjusted by hand')
#print ('CRJ700 has 11.07 [m2] vtail area, if we were to take the relative same thrust, we would have ~15 [m2] ')
