import input
import Class_2_estimation as Cl2
import numpy as np
import matplotlib.pyplot as plt

# Unchanged input parameters
g = input.g
c_t = 3.637 /1000/g   #          # [kg/N/s]                      #Specific fuel consumption #0.371*0.22481/3600 # 0.371 (TF34-GE-100A) lb/(h*lbf)) #c_t = mf *g/T
H = input.H                           # Heating value of the fuel (Hydrogen) #[J/kg] #heating value of hydrogen (lower and higher)
rho_c = input.rho_c                       #Density at 10km cruise altitude
rho_0 = input.rho0
S = Cl2.S
T_to = input.Tto 
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
rho_to = 1.05807  #rho_0                #Can be changed, depending on the runway height Both from 
nmax = 2.6                            # obtain nmax from max_load_factor function
V_nmax = 96.5

rho_land = rho_to                #might be changed
CL_land = CLmax            
CD_land = CD0_landGD + CL_land**2  / (np.pi * A * e)
T_to_eq = T_to / (rho_to / rho_0)**0.75


W1 = MTOW*(1-0.03346922107635999)
W2 = MTOW*(1-0.2905605894506181)
Cl_cr = 0.4567
l_over_d_cr = 15.6371
Cd_cr = Cl_cr/l_over_d_cr

from input import atmosphere_calculator as ac
def max_altitude(W1,W2,rho_c):
    rho_max_alt = W2/W1*rho_c
    h = 10000
    while ac(h)[2]/rho_max_alt > 1.00001:
        h = h +0.01
    hmax = h   
    return rho_max_alt, hmax
rho_max_alt, hmax = max_altitude(W1, W2, rho_c)


def performance_diagrams(rho_0, rho, e, clmax=CLmax, mtow=MTOW,s=S,a=A,MTOW=MTOW,cd0=CD0,t_to=T_to):
#   
#   this function constructs the performance diagrams for a particular flight condition you specify by changing the parameter rho the density at 
#    your altitude of interest. assumptions: steady symmetric flight and small angle approximmation
#    
#    t = available thrust at a certain level flight condition (l=w) [n], compensate for altitude if needed for the flight condition of interest
#    cd0 = zero lift drag coefficient at a certain configuration [-]
#    rho = density at altitude of certain flight condition [kg/m^3]
#    s = wing surface area [m^2]
#    w = weight at certain flight condition [n]
#    a = wing aspect ratio [-]
#    e = oswald efficiency factor at certain configuration
#    
#    pa = power available [w]
#    preq = power required [w]
#    d = drag [n]
#    v = velocity range for which the performance diagrams are constructed [m/s]
#   
    w = MTOW
    v = np.arange(1,301,1)
    #t = t_to * (rho / rho_0)**(3/4)
    #pa = np.array([t*i for i in v])
    d = np.array([cd0 * 0.5 * rho * i**2 * s + 2 * w**2 / (np.pi * a * e * rho * i**2 * s) for i in v])
    preq = np.array([cd0 * 0.5 * rho * i**3 * s + 2 * w**2 / (np.pi * a * e * rho * i * s) for i in v])
    
    t_cw1 = t_to * (rho_c / rho_0)#**(3/4)
    t_cw2 = t_to * (rho_max_alt / rho_0)#**(3/4)
    d_cw1 = np.array([cd0 * 0.5 * rho_c * i**2 * s + 2 * W1**2 / (np.pi * a * e * rho_c * i**2 * s) for i in v])
    d_cw2 = np.array([cd0 * 0.5 * rho_max_alt * i**2 * s + 2 * W2**2 / (np.pi * a * e * rho_max_alt * i**2 * s) for i in v])
    #d_cw1 = np.array([Cd_cr*0.5*rho_c*i**2*s for i in v])
    #d_cw2 = np.array([Cd_cr*0.5*rho_max_alt*i**2*s for i in v])

    
    vmin = np.sqrt(mtow * 2 /(s * rho_to * clmax))
    #vmin_sl = np.sqrt(mtow * 2 /(s * 1.225 * clmax))
    vmins1 = np.sqrt(W1*2/(s*rho_c*clmax))
    vcruise = np.sqrt(W2*2/(Cl_cr*rho_max_alt*s))
    
    #d = np.ma.masked_where(d < vmin, d)
    preq = np.ma.masked_where(vmin < 63, preq)
    
    fig, axs = plt.subplots(1)
    #fig.suptitle('jet performance diagrams')

    #d_sl= d_cw1
    #d = d_cw2
    
    #pa_sl = [t_to * v for v in v]
    #preq_sl = [d_sl * v for v in v]
   
    axs.axhline(y=t_cw1/1000, label='thrust at 10000m', color='lightblue')
    axs.plot(v, d_cw1/1000, label='drag at 10000m', color='red')
    axs.axvline(x=vcruise, label='Cruise speed', color='goldenrod')
    
    axs.axhline(y=t_cw2/1000, label='thrust at 12379m ', color='darkgreen')
    axs.plot(v*(np.sqrt(W2/W1)), d_cw2/1000, label='drag at 12379m ', color='darkred')
    axs.axvline(x=vmins1, label='Corresponding stall speeds', color='darkblue')
    axs.set(xlabel='airspeed [m/s]', ylabel='force [kn]')
    
    dmin1 = min(d_cw1)/1000
    dmin2 = min(d_cw2)/1000
    #print(np.where(dmin2 == min(dmin2)))
    
    axs.plot([0,300],[0,23.3857], label='line from origin to determine V/F' ,      linestyle='--', color = 'grey')  #20.0152249
    axs.plot([0,300],[0,27.2685],    linestyle='--', color = 'grey')
    axs.plot([149,174],[dmin2,dmin1],label = 'max endurance: min drag condition', color = 'limegreen')
    axs.legend()
    
    #axs[1].plot(v, v*t_to/1000000, label='power availabe at sl', color='limegreen')
    #axs[1].plot(v, v*d_sl/1000000, label='power required at sl', color='red')
    #axs[1].axvline(x=vmin_sl, label='stall speed at sl',linestyle='--', color='blue')
    
    #axs[1].plot(v, pa/1000000, label='power available at 1500m', color='darkgreen')
    #axs[1].plot(v, preq/1000000, label='power required at 1500m', color='darkred')   
    #axs[1].axvline(x=vmin, label='stall speed at 1500m',linestyle='--', color='darkblue')
    #axs[1].set(xlabel='airspeed [m/s]', ylabel='power [mw]')
    #axs[1].legend()
    
    np.savetxt('sl_drag.csv',d)
    
    axs.set_xlim([0,300])
    #axs[1].set_xlim([0,250])
    axs.set_ylim([0,t_cw1*2/1000])
    #axs[1].set_ylim([0, pa[-1]*1.2/1000000])
    plt.rc('axes', labelsize=8)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=8)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=8)    # fontsize of the tick labels
    plt.rc('legend', fontsize=8)    # legend fontsize
    plt.show()
    return 
   
performance_diagrams(rho_0, rho, e) 
"""
# for a specific point in time where W will be constant but undefined
def optimal_flight_condition(A,e,CD0,S=S,rho_c=rho_c,g=g,c_t=c_t):    #we seek to maximze V/F or minimize F/V

    
    T = T_to * (rho_c / rho_0)**(3/4)
    F =c_t*T
       # Still variable value                                
    
             
    #FoverV = 1/(1/(F)*np.sqrt(W/S*2/rho*1/Cl_c))  #Fuel flow | max thrust = constant i.e. as sealevel
    #beta = np.arctan(F_over_V)
    VoverF = Vopt/F # Maximize
    return VoverF,F,Vopt,Clopt

VoverF,F,Vopt,Clopt = optimal_flight_condition(A,e,CD0)
"""
#Note that W1/W2 are variable and need to be (re) determined
def max_range(H,Vcr,S,A,e,CD0,c_t,g=g,rho_c=rho_c): #only holds at constant altitude  
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
    
    
    Cl = np.sqrt(CD0*np.pi*A*e)
    Cd = 2*CD0
    LD_c = Cl/Cd
    #Clopt = np.sqrt(1/3*CD0*np.pi*A*e)
    #Cdopt = 4/3*CD0
    #LD_c_opt = Clopt/Cdopt
    
    
    W1 = MTOW*(0.9843800695598843) #1-0.03346922107635999
    W2 = MTOW*(0.9556007921574171) #1-0.2905605894506181
    W3 = W1
   
    #Vopt = np.sqrt(W/S*2/rho_c*1/Clopt)
    Vin = np.sqrt(W1*2/(S*rho_c*Cl)) # initial airspeed
    Vfin = np.sqrt(W2*2/(S*rho_c*Cl)) # final airspeed
    Vav = 2*Vin*(1-np.sqrt(W2/W1))/np.log(W1/W2)
    
    T = T_to * (rho_c / rho_0) #44626#
    F = c_t*T 
    eta_t = T*Vav/(F*H/g)
    #range_c_altitude = 2/c_t*np.sqrt(2/rho_c/S*Cl/Cd**2)*(np.sqrt(W1)-np.sqrt(W2))/1000 #km  
    #range_cruise = Vav/c_t*LD_c*np.log(W1/W2)/1000  
    range_unified = (eta_t*H/g*Cl/Cd*np.log(W1/W2))/1000 #km      
    #range_opt = Vopt/c_t*LD_c_opt*np.log(W1/W2)/1000
    R = 2*Vin*Cl/c_t/Cd*(1-np.sqrt(W2/W1))/1000 #Range at constant velocity and angle of attack in km
    E = 1/c_t*LD_c*np.log(W1/W2)
    print('maximum range at constant altitude:', R)
    print('maximum range at constant altitude with unified range equiation:', range_unified)
    print ('maximum endurance in hours equals:',E/3600)
    print (c_t*T/g*E,(W1-W2)/g)
    return
max_range(H,Vcr,S,A,e,CD0,c_t)



