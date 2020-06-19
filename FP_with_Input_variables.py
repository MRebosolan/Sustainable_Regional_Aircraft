import input
import Class_2_estimation as Cl2
import numpy as np
import matplotlib.pyplot as plt

# Unchanged input parameters
c_t = input.c_t                       #Specific fuel consumption #0.371*0.22481/3600 # 0.371 (TF34-GE-100A) lb/(h*lbf)) #c_t = mf *g/T
H = input.H                           # Heating value of the fuel (Hydrogen) #[J/kg] #heating value of hydrogen (lower and higher)
rho_c = input.rho_c                       #Density at 10km cruise altitude
rho_0 = input.rho0
g = input.g
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
    while ac(h)[2]/rho_max_alt > 1.001:
        h = h +1
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
    t = t_to * (rho / rho_0)**(3/4)
    pa = np.array([t*i for i in v])
    d = np.array([cd0 * 0.5 * rho * i**2 * s + 2 * w**2 / (np.pi * a * e * rho * i**2 * s) for i in v])
    preq = np.array([cd0 * 0.5 * rho * i**3 * s + 2 * w**2 / (np.pi * a * e * rho * i * s) for i in v])
    
    t_cw1 = t = t_to * (rho_c / rho_0)
    t_cw2 = t = t_to * (rho_max_alt / rho_0)
    d_cw1 = np.array([cd0 * 0.5 * rho_c * i**2 * s + 2 * W1**2 / (np.pi * a * e * rho_c * i**2 * s) for i in v])
    d_cw2 = np.array([cd0 * 0.5 * rho_max_alt * i**2 * s + 2 * W2**2 / (np.pi * a * e * rho_max_alt * i**2 * s) for i in v])
    
    vmin = np.sqrt(mtow * 2 /(s * rho_to * clmax))
    vmin_sl = np.sqrt(mtow * 2 /(s * 1.225 * clmax))
    
   
    #d = np.ma.masked_where(d < vmin, d)
    preq = np.ma.masked_where(vmin < 63, preq)
    
    fig, axs = plt.subplots(1, 2)
    #fig.suptitle('jet performance diagrams')

    d_sl= np.array([cd0 * 0.5 * rho_0 * i**2 * s + 2 * MTOW**2 / (np.pi * a * e * rho_0 * i**2 * s) for i in v])
    
    #pa_sl = [t_to * v for v in v]
    #preq_sl = [d_sl * v for v in v]
   
    axs[0].axhline(y=t_to/1000, label='thrust at sl', color='limegreen')
    axs[0].plot(v, d_sl/1000, label='drag at sl', color='red')
    
    axs[0].axhline(y=t/1000, label='thrust at 1500m ', color='darkgreen')
    axs[0].plot(v, d/1000, label='drag at 1500 m ', color='darkred')
    axs[0].axvline(x=vmin, label='stall speed at 1500m', linestyle='--', color='darkblue')
    axs[0].set(xlabel='airspeed [m/s]', ylabel='force [kn]')
    axs[0].legend()
    
    axs[1].plot(v, v*t_to/1000000, label='power availabe at sl', color='limegreen')
    axs[1].plot(v, v*d_sl/1000000, label='power required at sl', color='red')
    axs[1].axvline(x=vmin_sl, label='stall speed at sl',linestyle='--', color='blue')
    
    #axs[1].plot(v, pa/1000000, label='power available at 1500m', color='darkgreen')
    #axs[1].plot(v, preq/1000000, label='power required at 1500m', color='darkred')   
    #axs[1].axvline(x=vmin, label='stall speed at 1500m',linestyle='--', color='darkblue')
    #axs[1].set(xlabel='airspeed [m/s]', ylabel='power [mw]')
    #axs[1].legend()
    
    np.savetxt('sl_drag.csv',d)
    
    axs[0].set_xlim([0,250])
    axs[1].set_xlim([0,250])
    axs[0].set_ylim([0,t_to*2/1000])
    axs[1].set_ylim([0, pa[-1]*1.2/1000000])
    plt.rc('axes', labelsize=8)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=8)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=8)    # fontsize of the tick labels
    plt.rc('legend', fontsize=8)    # legend fontsize
    plt.show()
    return 
   
performance_diagrams(rho_0, rho, e) 


"""

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
    t = t_to * (rho / rho_0)
    pa = np.array([t*i for i in v])
    d = np.array([cd0 * 0.5 * rho * i**2 * s + 2 * w**2 / (np.pi * a * e * rho * i**2 * s) for i in v])
    preq = np.array([cd0 * 0.5 * rho * i**3 * s + 2 * w**2 / (np.pi * a * e * rho * i * s) for i in v])
    
    vmin = np.sqrt(mtow * 2 /(s * rho_to * clmax))
    vmin_sl = np.sqrt(mtow * 2 /(s * 1.225 * clmax))
    
   
    #d = np.ma.masked_where(d < vmin, d)
    preq = np.ma.masked_where(vmin < 63, preq)
    
    fig, axs = plt.subplots(1, 2)
    #fig.suptitle('jet performance diagrams')

    d_sl = np.genfromtxt('sl_drag.csv')
    d_sl = np.array([cd0 * 0.5 * rho_0 * i**2 * s + 2 * MTOW**2 / (np.pi * a * e * rho_0 * i**2 * s) for i in v])
    #pa_sl = [t_to * v for v in v]
    #preq_sl = [d_sl * v for v in v]
   
    axs[0].axhline(y=t_to/1000, label='thrust at sl', color='limegreen')
    axs[0].plot(v, d_sl/1000, label='drag at sl', color='red')
    
    axs[0].axhline(y=t/1000, label='thrust at 1500m ', color='darkgreen')
    axs[0].plot(v, d/1000, label='drag at 1500 m ', color='darkred')
    axs[0].axvline(x=vmin, label='stall speed at 1500m', linestyle='--', color='darkblue')
    axs[0].set(xlabel='airspeed [m/s]', ylabel='force [kn]')
    axs[0].legend()
    
    axs[1].plot(v, v*t_to/1000000, label='power availabe at sl', color='limegreen')
    axs[1].plot(v, v*d_sl/1000000, label='power required at sl', color='red')
    axs[1].axvline(x=vmin_sl, label='stall speed at sl',linestyle='--', color='blue')
    
    #axs[1].plot(v, pa/1000000, label='power available at 1500m', color='darkgreen')
    #axs[1].plot(v, preq/1000000, label='power required at 1500m', color='darkred')   
    #axs[1].axvline(x=vmin, label='stall speed at 1500m',linestyle='--', color='darkblue')
    #axs[1].set(xlabel='airspeed [m/s]', ylabel='power [mw]')
    #axs[1].legend()
    
    np.savetxt('sl_drag.csv',d)
    
    axs[0].set_xlim([0,250])
    axs[1].set_xlim([0,250])
    axs[0].set_ylim([0,t*1.5/1000])
    axs[1].set_ylim([0, pa[-1]*1.2/1000000])
    plt.rc('axes', labelsize=16)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=16)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=16)    # fontsize of the tick labels
    plt.rc('legend', fontsize=16)    # legend fontsize
    plt.show()
    return 
   
performance_diagrams(rho_0, rho, e) 
"""








    


