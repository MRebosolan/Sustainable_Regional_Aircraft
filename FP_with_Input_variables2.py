import input
import Class_2_estimation as Cl2
import numpy as np
import matplotlib.pyplot as plt
"""
Responsible person(s): Rick and Tobias

This requires as input:
    c_t                  | Specific fule consumption 
    H                    | Heating value of hydrogen
    rho_c                | Density at cruise altitude
    g                    | Gravtitational acceleration
    S                    | Wing area
    Tto                  | Thrust at take-off conditions
    V_c                  | Cruise speed
    AR                   | Aspect ratio
    CD0                  | Zero drag coefficient
    e                    | Efficiency factor
    MTOW                 | Maximum take-off weight
    Clmax_land           | Clmax in landing configuration
    Clmax_to             | Clmax at take-off 
    CD0_todg             | CD0 at take-off with gear down
    mu                   | Friction coefficient
    mu_br                | Friction coefficent braking
    h_sc                 | Screen height
    gamma_cl             | Climb angle
    gamma_a              | Approach angle
    MLW                  | Maximum landing weight
    h_sc                 | Screen height
    Trev                 | Total thrust reverse
    CD0_landGD           | CD0 in landing configuration with gear down
    
Input variables that are subject to change:
    W                    | Weight depending on where we are in flight
    rho                  | Variable rho, depending the altitude the calculations are performed
    nmax                 | Max load factor as obtained from graphs
    V_nmax               | Velocity at max load factor as obtained from graphs
    rho_to               | In case take-off is not at sealevel

This code gives as outputs:
    
    performance_diagrams:
        Plots thrust and drag against speed
        Plots Power available and power required against speed
    
    max_ROC_s:
        Calculated the maximum steady rate of climb
    
    max_load_factor_steepest_turn:
        Calculates max load factor experienced during steepest turn
    
    steepest_turn_and:
        Radius of steepest turn: (See function for more information)
    
    turn_radius:
        Calculates turn radius
    
    take_off_distances:
        s_to, x_airborne, x_tot | take-of fdistance, airborne distance and total distance respectively
   
    landing_distances:
         x_airborne, x_trans, x_brake, x_gr, x_tot, required_field_length
         airborne distance, transition distance, brake distance,
         groundrune distance, total x distance and required field length
         respectively
         
     optimal_flight_condition:
         VoverF,F,Vopt,Clopt | V/F, fuel consumption, Vopt (optimal speed)
         and Clopt (optimal cl) at a specific moment in time
         
     max_range:
         Calculates max range, gives 2 outputs based on the unified range 
         equation and the jet equation
         Currently not working that well; it is very sensitive to 
         the value of c_t (specific fuel consumption).
         
"""

# Unchanged input parameters
c_t = input.c_t                       #Specific fuel consumption #0.371*0.22481/3600 # 0.371 (TF34-GE-100A) lb/(h*lbf)) #c_t = mf *g/T
H = input.H                           # Heating value of the fuel (Hydrogen) #[J/kg] #heating value of hydrogen (lower and higher)
rho_c = input.rho_c                       #Density at 10km cruise altitude
rho_0 = input.rho0
g = input.g
S = Cl2.S
T_to = 121000 
T_1500m = T_to * (0.974 / rho_0)**(3/4)
Vcr = input.V_C_TAS
A = input.AR
CD0 = input.CD0
e = input.e
MTOW = input.MTOW
CLmax = input.CLmax_land              #Clmax at landing to determine stall speed
CLmax_to = input.CLmax_to
CL_to = CLmax_to/1.1**2 
CD0_togd =  input.CD0_togd            #same as for loading diagrams (taken from adsee)
mu = input.mu                         #from http://www.aerodynamics4students.com/aircraft-performance/take-off-and-landing.php
mu_brake = input.mu_br                #or0.16 #from https://www.icao.int/Meetings/grf2019/Documents/Presentations/GRF2019%20S8%20John%20Gadzinski%20-%20Four%20Winds%20Aerospace.pdf / https://www.stac.aviation-civile.gouv.fr/sites/default/files/stac/manifestation/tra2014/presentations/8-aircraft_breaking.pdf 
h_s = input.h_sc                      #screen height [m]
gamma_climb = input.gamma_cl          #taken from boeing 737  
gamma_app = input.gamma_ap            #gamma approach
MLW = input.MLW                       #Maximum landing weight [N]
T_rev = input.Trev                    #total thrust reverse during braking, TBD
CD0_landGD = input.CD0_landGD

#Parameters to be changed
W = 30000 * g 
rho = 0.974
rho_to = 0.974  #rho_0                        #Can be changed, depending on the runway height Both from 
nmax = 2.6                            # obtain nmax from max_load_factor function
V_nmax = 96.5

rho_land = rho_to                #might be changed
CL_land = CLmax            
CD_land = CD0_landGD + CL_land**2  / (np.pi * A * e)
T_to_eq = T_to / (rho_to / rho_0)**0.75
print(T_to_eq)
print(T_to_eq / MTOW)

#
#
#def performance_diagrams(T_to, CD0, rho_0, rho, S, W, A, e, CLmax, MTOW):
#    """ 
#    This function constructs the performance diagrams for a particular flight condition you specify by changing the parameter rho the density at 
#    your altitude of interest. Assumptions: steady symmetric flight and small angle approximmation
#    
#    T = available thrust at a certain level flight condition (L=W) [N], compensate for altitude if needed for the flight condition of interest
#    CD0 = zero lift drag coefficient at a certain configuration [-]
#    rho = density at altitude of certain flight condition [kg/m^3]
#    S = wing surface area [m^2]
#    W = weight at certain flight condition [N]
#    A = wing aspect ratio [-]
#    e = Oswald efficiency factor at certain configuration
#    
#    Pa = power available [W]
#    Preq = power required [W]
#    D = drag [N]
#    V = velocity range for which the performance diagrams are constructed [m/s]
#    """
#    V = np.arange(1,301,1)
#    T = T_to * (rho / rho_0)**(3/4)
#    Pa = np.array([T*i for i in V])
#    D = np.array([CD0 * 0.5 * rho * i**2 * S + 2 * W**2 / (np.pi * A * e * rho * i**2 * S) for i in V])
#    Preq = np.array([CD0 * 0.5 * rho * i**3 * S + 2 * W**2 / (np.pi * A * e * rho * i * S) for i in V])
#  
#    Vmin = np.sqrt(MTOW * 2 /(S * rho_to * CLmax))
#    Vmin_sl = np.sqrt(MTOW * 2 /(S * 1.225 * CLmax))
#    
#    
#    #D = np.ma.masked_where(D < Vmin, D)
#    Preq = np.ma.masked_where(Vmin < 63, Preq)
#    
#    fig, axs = plt.subplots(1, 2)
#    #fig.suptitle('Jet Performance Diagrams')
#
#    D_SL = np.genfromtxt('SL_drag.csv')
#    #Pa_SL = [T_to * v for v in V]
#    #Preq_SL = [D_SL * v for v in V]
#    
#    axs[0].axhline(y=T_to/1000, label='Thrust at SL', color='limegreen')
#    axs[0].plot(V, D_SL/1000, label='Drag at SL', color='red')
#    axs[0].axvline(x=Vmin_sl, label='Stall speed at SL',linestyle='--', color='blue')
#    
#    axs[0].axhline(y=T/1000, label='Thrust at 1500m ', color='darkgreen')
#    axs[0].plot(V, D/1000, label='Drag at 1500 m ', color='darkred')
#    axs[0].axvline(x=Vmin, label='Stall speed at 1500m', linestyle='--', color='darkblue')
#    axs[0].set(xlabel='Airspeed [m/s]', ylabel='Force [kN]')
#    axs[0].legend()
#    
#    axs[1].plot(V, V*T_to/1000000, label='Power availabe at SL', color='limegreen')
#    axs[1].plot(V, V*D_SL/1000000, label='Power required at SL', color='red')
#    axs[1].axvline(x=Vmin_sl, label='Stall speed at SL',linestyle='--', color='blue')
#    
#    axs[1].plot(V, Pa/1000000, label='Power available at 1500m', color='darkgreen')
#    axs[1].plot(V, Preq/1000000, label='Power required at 1500m', color='darkred')
#    axs[1].axvline(x=Vmin, label='Stall speed at 1500m',linestyle='--', color='darkblue')
#    axs[1].set(xlabel='Airspeed [m/s]', ylabel='Power [MW]')
#    axs[1].legend()
#    
#    np.savetxt('SL_drag.csv',D)
#    
#    axs[0].set_xlim([0,250])
#    axs[1].set_xlim([0,250])
#    axs[0].set_ylim([0,T*1.5/1000])
#    axs[1].set_ylim([0, Pa[-1]*1.2/1000000])
#
#    plt.rc('axes', labelsize=16)    # fontsize of the x and y labels
#    plt.rc('xtick', labelsize=16)    # fontsize of the tick labels
#    plt.rc('ytick', labelsize=16)    # fontsize of the tick labels
#    plt.rc('legend', fontsize=16)    # legend fontsize
#    plt.show()
#    return 
#    
#performance_diagrams(T_1500m, CD0_togd, rho_0, 0.974, S, W, A, e, CLmax_to, MTOW)   
    
 
def max_ROC_s(T_to, CD0, S, A, e, rho_0):
    """ 
    This function computes the maximum steady rate of climb achievable at different flight conditions for which you must manually select 
    the densities (change rho_lst) and aircraft weight (W_lst) for the altitudes of interest.
    
    T_to = Take off thrust at sea-level [N]
    CD0 = zero lift drag coefficient at a certain configuration [-]
    rho_lst = range of densities at corresponding flight altitudes for which you want to calculate the maximum steady rate of climb
    S = wing surface area [m^2]
    W_lst = range of weights at corresponding flight alitudes at which you want to calculate the maximum steady rate of climb [N]
    A = wing aspect ratio [-]
    e = Oswald efficiency factor at certain configuration
    h = range of altitudes at corresponding flight alitudes at which you want to calculate the maximum steady rate of climb 
    
    Pa = power available [W]
    Preq = power required [W]
    D = drag [N]
    V = velocity range for which the performance diagrams are constructed [m/s]
    """
    rho_lst = [1.225, 0.7365, 0.4135] 		#densities at sea level, 5 km, 10 km alitude 
    W_lst = [400000, 355000, 350000]					#weights at the corresponding above alitudes [N]!
    V = np.arange(1,301,1)
    output_lst = []
    for i in range(len(rho_lst)):
        T = (rho_lst[i] / rho_0)**(3/4) * T_to
        Pa = np.array([T*i for i in V])
        Preq = np.array([CD0 * 0.5 * rho_lst[i] * j**3 * S + 2 * W_lst[i]**2 / (np.pi * A * e * rho_lst[i] * j * S) for j in V])          
        Pa_min_Pr = list(np.subtract(Pa,Preq))
        max_ROC = max(Pa_min_Pr)/W_lst[i]
        idx_max_ROC = Pa_min_Pr.index(max(Pa_min_Pr))
        text = 'max steady ROC at altitude with density of', rho_lst[i], '=', max_ROC, 'm/s with corresponding airspeed of', V[idx_max_ROC], 'm/s'
        output_lst.append(text)
    return output_lst

#print(max_ROC_s(T_to, CD0, S, A, e, rho_0))

n = np.linspace(1,2.6,6)

def max_load_factor_steepest_turn(T_to, CD0, rho_0, rho, S, W, A, e, CLmax):
    """ 
    This function plots the maximum achievable load factor in a turn against airspeed. Adjust CLmax to the right aircraft configuration and also take into account
    the structural nmax that you don't exceed that value. Read off nmax at various speed and whit that compute the steepest turn. 
    
    """
    V = np.arange(1,301,1)
    T = T_to * (rho / rho_0)**(3/4)
    V_stall_lst = []
    for j in range(len(n)):
        D = j * np.array([CD0 * 0.5 * rho * i**2 * S + 2 * W**2 / (np.pi * A * e * rho * i**2 * S) for i in V])
        Preq = j * np.array([CD0 * 0.5 * rho * i**3 * S + 2 * W**2 / (np.pi * A * e * rho * i * S) for i in V])
        
        V_stall_lst.append(np.sqrt(n[j] * W * 2 / (S * rho * CLmax)))
        D = np.ma.masked_where(V < V_stall_lst[j], D)
        Preq = np.ma.masked_where(V < V_stall_lst[j], Preq)           
        
        plt.scatter(V, D, s=1, label='Drag [N], n ='+str(round(n[j],1)), color='red')
     
    plt.axhline(y=T, label='Thrust [N]', color='green')   
    plt.legend()
    plt.show()
    return 

#max_load_factor_steepest_turn(T_to, CD0, rho_0, rho, S, W, A, e, CLmax)

def steepest_turn_and(V, nmax, g=g):
    """
    This function computes the radius for the steepest turn. Get the maximum load factor at the speed of interest from the 
    graph constucted by the function max_load_factor_steepest_turn. This function can also be used to compute the minimum turn radius.
    For this, readd off the nmax at various speeds from the above mentioned graph and compute the R with this function for each of the combinations.
    See what yields the minimum turn radius R (should be at a slightly lower speed than at which nmax is achieved in a turn). Finally, the minimum time 
    to turn can be calculated with the minimum turn radius by doing 2piR/V (take V where Rmin is achieved).This speed should be slightly higher than the 
    speed at which minimum turn radius is achieved but slightly lower than speed at which nmax is achieved. All of this is calculated for one specific alitude
    which you can change by varying the input parameter rho.
    """
    R = V**2 / (g * np.sqrt(nmax**2 - 1))
    return R



def turn_radius(V, rate, g=g):
     """ 
     This function computes the turn radius, bank angle and load factor for a particular rate # turn (so for rate use a value incidacting 
     a rate ... turn) at a particular airspeed of interest [m/s]
     angular = angular velocity during the turn [rad/s]
     R = turn radius [m]
     phi = bank angle [rad]
     n = load factor during the turn [-]
     
     """
     angular = rate * np.radians(3)      
     R = V / angular
     phi = np.arctan(V**2 / (g * R))
     n = 1 / np.cos(phi)
     return R, phi, n



def dimensional(c,rho,V,S):             #converts aerdynamic coefficient into aerodynamic force
    return 0.5 * rho * V**2 * S * c

#---------------------------------------------------------
#Take-off and landing performance
#---------------------------------------------------------

"""take-off and landing input dummy parameters"""                                                             #[-]
                        #[kg/m^3]
CD_togd = CD0_togd + CL_to**2 / (np.pi * A * e)                         #[-]

def take_off_distances(CLmax, CD_togd, rho_to, S, CL_to, MTOW, T_to, mu, g, h_s, gamma_climb):
    """
    This function calculates the take-off ground run distance (s_to), airborne distance(x_airborne) and horizontal distance covered until 
    screen height is reached (x_tot) at maximum take off weight. 
    Assumptions: no ground effect, no runway slope, no wind
    
    Vmin = stall speed [m/s]
    Vlof = lift off speed [m/s]
    Vbar = average speed during take off [m/s]
    D_to = average drag force during take off [N]
    L_to = average lift force during take off [N]
    a_bar = average acceleration during take off [m/s^2]
    R = transition radius [m]
    h_trans = height gained during transition [m]
    
    The thrust seems to high as the take off distances are very low, also may reconsider friction coefficients
    """
    Vmin = np.sqrt(MTOW * 2 /(S * rho_to * CLmax))                  #or change to Vs, depends on whether the wing is sized for Vs or for CLmax
    Vlof = 1.05 * Vmin 
    print('Vlof = ', Vlof)                                           
    Vbar = Vlof / np.sqrt(2)                                       
    D_to = dimensional(CD_togd, rho_to, Vbar, S)                    
    L_to = dimensional(CL_to, rho_to, Vbar, S)
    T_bar = T_to / np.sqrt(2)                      
    a_bar = g / MTOW * (T_bar - D_to - mu * (MTOW - L_to))           #
    
    s_to = Vlof**2 / (2 * a_bar)                                    
    print('x_gr_to = ', s_to)
    R = Vlof**2 / (0.15 * g)
    x_trans = R * np.sin(gamma_climb)
    print('x_tr = ', x_trans)
    h_trans = R * (1 - np.cos(gamma_climb))
    if h_trans < h_s:
        x_climb = (h_s - h_trans) / np.tan(gamma_climb)
    else:
        x_climb = 0
    print('x_cl = ', x_climb)
    x_airborne = x_trans + x_climb
    x_tot = s_to + x_airborne
    print('x_to = ', x_tot)
    return s_to, x_airborne, x_tot
    
print(take_off_distances(CLmax, CD_togd, rho_to, S, CL_to, MTOW, T_to, mu, g, h_s, gamma_climb))

def landing_distances(MLW, S, rho_land, CLmax, gamma_app, h_s, CD_land, CL_land, T_rev, mu_brake, g=g):
    R = 1.3**2 * (MLW * 2 /(S * rho_land * CLmax)) / (0.1 * g)
    print('R_la = ', R)
    x_airborne = R * np.sin(gamma_app) + (h_s - (1 - np.cos(gamma_app)) * R) / np.tan(gamma_app)
    print('x_ap = ',R * np.sin(gamma_app) )
    print('x_fl = ', x_airborne - R * np.sin(gamma_app))
    Vmin = np.sqrt(MLW * 2 /(S * rho_land * CLmax))
    V_app = 1.3 * Vmin
    print('V_ap = ', V_app)
    x_trans = 2.6 * Vmin
    print('x_rot = ', x_trans)
    V_bar = V_app / np.sqrt(2)
    D_land = dimensional(CD_land, rho_land, V_bar, S)
    L_land = dimensional(CL_land, rho_land, V_bar, S)
    T_rev_bar = T_rev / np.sqrt(2)
    x_brake = MLW**2 / (2 * g * S) * 2 / rho_land * 1.3**2 / CLmax * 1 / (T_rev_bar + D_land + mu_brake * (MLW - L_land))
    print('x_br = ', x_brake)
    x_gr = x_trans + x_brake
    print('x_gr_la =', x_gr)
    required_field_length = 10 / 6 * x_gr
    print('RFL = ', required_field_length)
    x_tot = x_airborne + x_gr
    print('x_la = ', x_tot)
    
    return x_airborne, x_trans, x_brake, x_gr, x_tot, required_field_length, V_bar, V_app

x_airborne, x_trans, x_brake, x_gr, x_tot, required_field_length, V_bar, Vap = landing_distances(MLW, S, rho_land, CLmax, gamma_app, h_s, CD_land, CL_land, T_rev, mu_brake)
#print(landing_distances(MLW, S, rho_land, CLmax, gamma_app, h_s, CD_land, CL_land, T_rev, mu_brake, g=g))
print('Approach speed: ', Vap)




# for a specific point in time where W will be constant but undefined
def optimal_flight_condition(A,e,CD0,S=S,rho_c=rho_c,g=g,c_t=c_t):    #we seek to maximze V/F or minimize F/V
    """
    This would be the optimal velocity, Cl and Cd at one specific
    moment in time.
    """
    T = T_to * (rho_c / rho_0)**(3/4)
    F =c_t*T
    W = 31500*g    # Still variable value                                
    Clopt = np.sqrt(1/3*CD0*np.pi*A*e)
    Cdopt = 4/3*CD0
    Vopt = np.sqrt(W/S*2/rho_c*1/Clopt)
             
    #FoverV = 1/(1/(F)*np.sqrt(W/S*2/rho*1/Cl_c))  #Fuel flow | max thrust = constant i.e. as sealevel
    #beta = np.arctan(F_over_V)
    VoverF = Vopt/F # Maximize
    return VoverF,F,Vopt,Clopt

VoverF,F,Vopt,Clopt = optimal_flight_condition(A,e,CD0)

#Note that W1/W2 are variable and need to be (re) determined
def max_range(H,Vcr,F,S,A,e,CD0,c_t,g=g,rho_c=rho_c): #only holds at constant altitude  
    """
    eta_t is the total efficiency
    Cl due to operation limitations 
    Cd idem.
    W1 is the initial weight (W4)
    W2 is the final weight (W5)
    Currently the values it prints are quite off.... they are very
    sensitive to the value of c_t. Also, the value of F is rather high
    due to the fact that the heating value H is 3 times higher for
    hydrogen than for kerosene
    """
    T = T_to * (rho_c / rho_0)**(3/4)
    F = c_t*T 
    eta_t = T*Vcr/(F*H/g)
    Cl = np.sqrt(CD0*np.pi*A*e)
    Cd = 2*CD0
    W1 = 33000*g
    W2 = 30000*g

    Range = 2/(c_t*Cd)*np.sqrt(1/S*2/rho_c*Cl)*(np.sqrt(W1)-np.sqrt(W2))    
    Range2 = eta_t*H/g*Cl/Cd*np.log(W1/W2)  
    return Range, Range2,eta_t,Cl,Cd,F

Range,Range2,eta_t,Cl,Cd,F = max_range(H,Vcr,F,S,A,e,CD0,c_t)

#print (Range,Range2,eta_t,Cl,Cd)
#
#    #---------------------------------------------------------
##Cruise economics
##---------------------------------------------------------



W_pay = input.W_payload
R_des = input.Design_range
OEW = input.OEW
cj_c = input.cj_c
LD_c = input.LD_c 
percent = 67
W_fuel_max_pay = Cl2.M_fuel_kg
W_fuel_max = W_fuel_max_pay / percent * 100
W_fuel_exchange = W_fuel_max - W_fuel_max_pay
W_pay_max_fuel = W_pay - W_fuel_exchange

end_i = [0.99651,0.99651,0.998255,0.99302,0.99651,0.99651,0.9971422253661717,0.9939283669854275,0.99651,0.997208]

Mff_D = 1 - W_fuel_max_pay / OEW

end5_D = Mff_D

for i in end_i:
    end5_D = end5_D / i

#
#def payload_range_diagram():
#    R_AB = np.arange(0, R_des + 1, 1)
#    W_AB = (R_des+1) * [W_pay]
#    W_C = W_pay_max_fuel
#    W_D = 0
#    R_C = Vcr / g / cj_c * LD_c * np.log(1/end5_C) / 1000
#    R_D = Vcr / g / cj_c * LD_c * np.log(1/end5_D) / 1000
#    
#    print(R_C, R_D)
#    R_BC = np.arange(R_des, R_C + 1, 1)
#    W_BC = [W_pay + (W_C - W_pay) / (R_C - R_des) * (R-R_des) for R in R_BC]
##
#    R_CD = np.arange(R_C, R_D + 1, 1)
#    W_CD = [W_C + (W_D - W_C) / (R_D - R_C) * (R-R_C) for R in R_CD]
#    print(R_CD)
#    plt.figure()
#    plt.plot(R_AB, W_AB, color='green', label='AB')
#    plt.plot(R_BC, W_BC, color='orange', label='BC')
#    plt.plot(R_CD, W_CD, color='blue', label='CD')
#    plt.legend()
#    plt.show()
#    return R_AB, W_AB, R_BC, W_BC, R_CD, W_CD
    
#R_AB, W_AB, R_BC, W_BC, R_CD, W_CD = payload_range_diagram()   
#
#def payload_range_diagram():
#    R_AB = np.arange(0, R_des + 1, 1)
#    W_AB = (R_des+1) * [W_pay]
#    
#    
#    R_D = Vcr / g / cj_c * LD_c * np.log(1/end5_D) / 1000
#    print(R_D)
#    R_BD = np.arange(R_des, R_D + 1, 1)
#    W_D = 0
#    W_BD = [W_pay - (W_pay - W_D) / (R_D - R_des) * (R - R_des) for R in R_BD]
#    
#    plt.figure()
#    plt.plot(R_AB, W_AB, color='green', label='Max payload')
#    plt.plot(R_BD, W_BD, color='blue', label='Exchange payload')
#    plt.xlim([0,5000])
#    plt.ylim([0,9000])
#    plt.xlabel('Range [km]')
#    plt.ylabel('Payload weight [kg]')
#    plt.rc('axes', labelsize=16)    # fontsize of the x and y labels
#    plt.rc('xtick', labelsize=16)    # fontsize of the tick labels
#    plt.rc('ytick', labelsize=16)    # fontsize of the tick labels
#    plt.rc('legend', fontsize=16)    # legend fontsize
#    plt.legend()
#    plt.show()
#    return 
#
#payload_range_diagram()


    
#
#airport_time = 3600   #assumed turn around time https://www.ikusi.aero/en/blog/what-tasks-are-performed-during-turnaround-time-aircraft
#
#climb_time = 3600 / 3 #calculate the climb time from the max steady rate of climb and assume same for the differences in height
#descent_time = 3600 / 2.5   #add a bit of time w.r.t. climb_time
#
#delta_t = airport_time + climb_time + descent_time
#
#V_block = R_des / (R_des / Vcr + delta_t)
#pax_km_h = V_block * W_pay * 3.6
#print(pax_km_h)

#todo: make payload range diagram, and transport productivity graph by calculating range with range equation

    


    


