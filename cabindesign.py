# -*- coding: utf-8 -*-
"""
Created on Wed Jun  3 13:43:05 2020

@author: malfl
"""

import numpy as np
import math
from hydrogen_tank_sizing import tank_sizing
from hydrogen_tank_sizing import tank_sizing_fuselage
from fuselage_weight_estimation import W_fuselage_torenbeek
import input

def cabin_design(fractioninfus,fractionintail,HYDROGENVOLUME,top_selecter,podlength=5):
    #PASSENGER SECTION
    #LARGELY BASED ON AIRBUS A220
    
    wseat=18.9*0.0254
    warmrest=2*0.0254
    waisle=0.5
    wclearance=0.01
    seat_pitch=32*0.0254
    emergency_clearance=10*0.0254
    
    h_aisle=2
    
    hheadroom=1.5
    totalabovefloor=0.68
    wall_thickness=0.1 #ADSEE
    floor_thickness=0.1 #ADSEE
    
    three=round(wseat*3+wclearance+warmrest*4,3)
    two=round(wseat*2+wclearance+warmrest*3,3)
    
    total=round(two+three+waisle,3)
    wheadroom=round(total-2*(wclearance+warmrest)-wseat,3)
    wfloor=round(total-0.18,3)
    
    
    
    #two types I exits
    #two types III exits
    paxsectionlength=round(16*seat_pitch+emergency_clearance,3)
    aft_galley_length=1.2 #BIT SMALLER THAN A220
    front_galley_length=2.5 #BIT SMALLER THAN A220
    front_aisle=0.8 #A220
    totalcabinlength=paxsectionlength+aft_galley_length+front_galley_length+front_aisle
    cockpit_length=3.6#BIT SMALLER THAN A220
    
    outer_diameter=3.486 #A220
    inner_diameter=outer_diameter-wall_thickness
    V_tank=HYDROGENVOLUME #CLASS I
    R_tank_fus=outer_diameter/2-0.15
    R_tank_tail=1
    
    #print(totalcabinlength)
    
        
    lf=28 # derived from A220, adapt in input file, total length of fuselage
    
    rho_hydrogen=70
        
    V_tank_cyl=(fractioninfus-fractionintail*fractioninfus)*V_tank
    V_tank_tail=fractionintail*fractioninfus*V_tank
    V_tank_top=(V_tank-V_tank_cyl-V_tank_tail)*top_selecter
    V_tank_pod=(V_tank-V_tank_cyl-V_tank_tail)*(1-top_selecter)
        
    #CYLINDER STORAGE
    t_cyl,m_cyl, tm_cyl, d_cyl,l_cyl=tank_sizing_fuselage(V_tank_cyl,R_tank_fus,1)
    if V_tank_cyl==0:
        t_cyl,m_cyl, tm_cyl, d_cyl,l_cyl=0,0,0,0,0
        
    print('CYLINDER TANK: ','| mass: ',tm_cyl,'| diameter: ',d_cyl,'| length: ',l_cyl)
    
    #TAIL STORAGE
    t_tail,m_tail, tm_tail, d_tail,l_tail=tank_sizing_fuselage(V_tank_tail,R_tank_tail,1)
    if V_tank_tail==0:
        t_tail,m_tail, tm_tail, d_tail,l_tail=0,0,0,0,0
    print('TAIL TANK: ','| mass: ',tm_tail,'| diameter: ',d_tail,'| length: ',l_tail)
    
    #TOP STORAGE
    t_top,m_top,tm_top,d_top,l_top=tank_sizing(V_tank_top,totalcabinlength+l_cyl,2)
    if V_tank_top==0:
        t_top,m_top, tm_top, d_top,l_top=0,0,0,0,0
    print('TOP TANK: ','| mass: ',tm_top,'| diameter: ',d_top,'| length: ',l_top)

    #POD STORAGE
    print(V_tank_pod)
    t_pod,m_pod,tm_pod,d_pod,l_pod=tank_sizing(V_tank_pod,podlength*0.9*2,2) #allow for some cones
    if V_tank_pod==0:
        t_pod,m_pod, tm_pod, d_pod,l_pod=0,0,0,0,0
    print('POD TANK (for each): ','| mass: ',tm_pod/2,'| diameter: ',d_pod,'| tank length: ',l_pod/2)

    
    tm_tanksystem=tm_cyl+tm_tail+tm_top+tm_pod
    
    CGtank=((tm_cyl)*(totalcabinlength+l_cyl/2)+(tm_top)*(totalcabinlength/2+l_cyl/2)+(tm_tail)\
            *(totalcabinlength+l_cyl+l_tail/2)+(tm_pod)*(totalcabinlength/2))/tm_tanksystem
            
    CGfuelfull=((V_tank_cyl)*(totalcabinlength+l_cyl/2)+(V_tank_top)*(totalcabinlength/2+l_cyl/2)\
                +(V_tank_tail)*(totalcabinlength+l_cyl+l_tail/2)+(V_tank_pod)*(totalcabinlength/2))/(V_tank_pod+V_tank_cyl+V_tank_tail+V_tank_top)
    
    CGcomb=(CGtank*tm_tanksystem+CGfuelfull*(V_tank_cyl+V_tank_tail+V_tank_top+V_tank_pod)\
            *rho_hydrogen)/(tm_tanksystem+(V_tank_cyl+V_tank_tail+V_tank_top+V_tank_pod)*rho_hydrogen)


    widthf=outer_diameter
    if (d_top+1.55)<=widthf/2:
        hf=widthf
    else:
        hf=widthf/2+1.55+d_top#CATIA  
        
    print(hf,widthf)
    #AERODYNAMICS ###VERIFIED
    lf+=l_cyl
    lambdaf=lf/(hf) #TORENBEEK
    Perimeter=np.pi*outer_diameter+np.pi*d_top*0.85
    fuselage_area=Perimeter*(lf)*(1-2/lambdaf)**(2/3)*(1+1/lambdaf**2) #TORENBEEK, and extra skin surface due to top
    FFbody=1+2.2/lambdaf**1.5+3.8/lambdaf**3 #https://www.fzt.haw-hamburg.de/pers/Scholz/HOOU/AircraftDesign_13_Drag.pdf,https://arc.aiaa.org/doi/pdf/10.2514/1.47557?casa_token=Ba2QtSu7zucAAAAA:aOfBQWl3BJR2ssM7K9PD_LIeIrlhvPDfImqpjJwciE4oqVkbmIZ-AANSbmtXX6CmqAjgX6VQ0O0
    k_v=7.95*10**(-6) #kinematic viscosity at 206 KELVIN
    MACH=input.mach_cruise
    Reynolds=input.V_C*0.51444*lf/k_v
    Cfturb=0.455/((math.log(Reynolds,10))**2.58*(1+0.144*MACH**2)**0.65)
    CDzerofus=Cfturb*FFbody*fuselage_area/70 #ref area CRJ700
    
    if V_tank_pod!=0:
        lambdaf_pods=podlength/d_pod
        pylon=2*1*0.5
        pod_area=np.pi*d_pod*(podlength-d_pod/2)+2*4*np.pi*d_pod**2/4+pylon#put hemispheres on ends
        FFpods=1+0.35/lambdaf_pods
        Reynolds_pods=input.V_C*0.51444*podlength/k_v
        Cfturb_pods=0.455/((math.log(Reynolds_pods,10))**2.58*(1+0.144*MACH**2)**0.65)
        CDzeropods=Cfturb_pods*FFpods*pod_area/70
        poddrag=2*0.5*input.rho_c*(input.V_C*0.514444)**2*(pod_area)*CDzeropods*(1-top_selecter) #two pods
        print('pods: ',poddrag)
    else: poddrag,Cfturb_pods,CDzeropods,Reynolds_pods,pod_area,lambdaf_pods,FFpods=0,0,0,0,0,0,0
    
    
    fusdrag=0.5*input.rho_c*(input.V_C*0.514444)**2*fuselage_area*(CDzerofus)
    
    
    totdrag=fusdrag+poddrag #ONLY OF FUSELAGE AND PODS
    #AERODYNAMICS
    
    
    

        
        
    lh=(totalcabinlength+l_cyl)/2+3 #GUESS
    fuselage_weight=W_fuselage_torenbeek(input.V_dive, lh, widthf/0.3048, hf/0.3048, fuselage_area/0.3048/0.3048)    
    
    
    return(t_cyl,m_cyl, tm_cyl, d_cyl,l_cyl,t_tail,m_tail, tm_tail, d_tail,l_tail\
           ,t_top,m_top,tm_top,d_top,l_top,t_pod,m_pod,tm_pod,d_pod,l_pod,totalcabinlength,V_tank_cyl, V_tank_tail, V_tank_top,V_tank_pod,\
           tm_tanksystem,CGtank,CGfuelfull,CGcomb,totdrag,fuselage_weight,CDzerofus,FFbody,Cfturb,fuselage_area,CDzeropods,fusdrag,poddrag)






    

