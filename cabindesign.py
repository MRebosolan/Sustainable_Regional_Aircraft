# -*- coding: utf-8 -*-
"""
Created on Wed Jun  3 13:43:05 2020

@author: malfl
"""

import numpy as np

from hydrogen_tank_sizing import tank_sizing
from hydrogen_tank_sizing import tank_sizing_fuselage
from fuselage_weight_estimation import W_fuselage_torenbeek
import input

def cabin_design(fractioninfus,fractionintail,HYDROGENVOLUME,top_selecter,podlength=3):
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
    
    print(totalcabinlength)
    
        
    lf=28 # derived from A220, adapt in input file, total length of fuselage
    
    rho_hydrogen=70
        
    V_tank_cyl=(fractioninfus-fractionintail*fractioninfus)*V_tank
    V_tank_tail=fractionintail*fractioninfus*V_tank
    V_tank_top=(V_tank-V_tank_cyl-V_tank_tail)*top_selecter
    V_tank_pod=V_tank-V_tank_cyl-V_tank_tail*(1-top_selecter)
        
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
    t_pod,m_pod,tm_pod,d_pod,l_pod=tank_sizing(V_tank_pod,3,2)
    if V_tank_pod==0:
        t_pod,m_pod, tm_pod, d_pod,l_pod=0,0,0,0,0
    print('POD TANK: ','| mass: ',tm_pod,'| diameter: ',d_pod,'| length: ',l_pod)

    
    tm_tanksystem=tm_cyl+tm_tail+tm_top+tm_pod
    
    CGtank=((tm_cyl)*(totalcabinlength+l_cyl/2)+(tm_top)*(totalcabinlength/2+l_cyl/2)+(tm_tail)\
            *(totalcabinlength+l_cyl+l_tail/2)+(tm_pod)*(totalcabinlength/2))/tm_tanksystem
            
    CGfuelfull=((V_tank_cyl)*(totalcabinlength+l_cyl/2)+(V_tank_top)*(totalcabinlength/2+l_cyl/2)\
                +(V_tank_tail)*(totalcabinlength+l_cyl+l_tail/2)+(V_tank_pod)*(totalcabinlength/2))/(V_tank_cyl+V_tank_tail+V_tank_top)
    
    CGcomb=(CGtank*tm_tanksystem+CGfuelfull*(V_tank_cyl+V_tank_tail+V_tank_top+V_tank_pod)\
            *rho_hydrogen)/(tm_tanksystem+(V_tank_cyl+V_tank_tail+V_tank_top+V_tank_pod)*rho_hydrogen)
    
    
    lambdaf=lf/outer_diameter #TORENBEEK
    fuselage_area=np.pi*outer_diameter*(lf+l_cyl)*(1-2/lambdaf)**(2/3)*(1+1/lambdaf**2)+l_top*d_top*np.pi #TORENBEEK, and extra skin surface due to top
    lh=(totalcabinlength+l_cyl)/2+3 #GUESS
    widthf=input.widthf
    
    if (d_top+1.55)<=widthf/2:
        hf=widthf
    else:
        hf=widthf/2+1.55+d_top#CATIA
    
    fuselage_weight=W_fuselage_torenbeek(input.V_dive, lh, widthf/0.3048, hf/0.3048, fuselage_area/0.3048/0.3048)
    
    CDzerofus=0.01
    CDzerobubble=0.01*d_top/outer_diameter
    CDzeropods=0.005
    
    fusdrag=0.5*input.rho_c*(input.V_C*0.514444)**2*fuselage_area*(CDzerofus+CDzerobubble)
    poddrag=0.5*input.rho_c*(input.V_C*0.514444)**2*(d_top**2/4*np.pi+np.pi*d_top*l_top)*CDzeropods
    totdrag=fusdrag+poddrag #ONLY OF FUSELAGE AND PODS
        
    
    
    
    return(t_cyl,m_cyl, tm_cyl, d_cyl,l_cyl,t_tail,m_tail, tm_tail, d_tail,l_tail\
           ,t_top,m_top,tm_top,d_top,l_top,totalcabinlength,V_tank_cyl, V_tank_tail, V_tank_top,tm_tanksystem,CGtank,CGfuelfull,CGcomb,totdrag,fuselage_weight)








    

