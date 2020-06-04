# -*- coding: utf-8 -*-
"""
Created on Wed Jun  3 13:43:05 2020

@author: malfl
"""

import numpy as np

from hydrogen_tank_sizing import tank_sizing
from hydrogen_tank_sizing import tank_sizing_fuselage

def cabin_design(fractioninfus,fractionintail,HYDROGENVOLUME):
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
    wall_thickness=0.1
    floor_thickness=0.1
    
    three=round(wseat*3+wclearance+warmrest*4,3)
    two=round(wseat*2+wclearance+warmrest*3,3)
    
    total=round(two+three+waisle,3)
    wheadroom=round(total-2*(wclearance+warmrest)-wseat,3)
    wfloor=round(total-0.18,3)
    
    
    
    #two types I exits
    #two types III exits
    paxsectionlength=round(16*seat_pitch+emergency_clearance,3)
    aft_galley_length=0.9
    front_galley_length=1.8
    totalcabinlength=paxsectionlength+aft_galley_length+front_galley_length
    
    outer_diameter=3.486
    inner_diameter=3.286
    V_tank=HYDROGENVOLUME #CLASS I
    R_tank_fus=outer_diameter/2-0.15
    R_tank_tail=1
    
    
        
    
    rho_hydrogen=70
        
    V_tank_cyl=(fractioninfus-fractionintail*fractioninfus)*V_tank
    V_tank_tail=fractionintail*fractioninfus*V_tank
    V_tank_top=V_tank-V_tank_cyl-V_tank_tail
        
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
    tm_tanksystem=tm_cyl+tm_tail+tm_top
    CGtank=((tm_cyl)*(totalcabinlength+l_cyl/2)+(tm_top)*(totalcabinlength/2+l_cyl/2)+(tm_tail)*(totalcabinlength+l_cyl+l_tail/2))/tm_tanksystem
    CGfuelfull=((V_tank_cyl*rho_hydrogen)*(totalcabinlength+l_cyl/2)+(V_tank_top*rho_hydrogen)*(totalcabinlength/2+l_cyl/2)+(V_tank_tail*rho_hydrogen)*(totalcabinlength+l_cyl+l_tail/2))/(V_tank_cyl+V_tank_tail+V_tank_top)/rho_hydrogen
    CGcomb=(CGtank*tm_tanksystem+CGfuelfull*(V_tank_cyl+V_tank_tail+V_tank_top)*rho_hydrogen)/(tm_tanksystem+(V_tank_cyl+V_tank_tail+V_tank_top)*rho_hydrogen)
    return(t_cyl,m_cyl, tm_cyl, d_cyl,l_cyl,t_tail,m_tail, tm_tail, d_tail,l_tail,t_top,m_top,tm_top,d_top,l_top,totalcabinlength,V_tank_cyl, V_tank_tail, V_tank_top,tm_tanksystem,CGtank,CGfuelfull,CGcomb)






t_cyl,m_cyl, tm_cyl, d_cyl,l_cyl,t_tail,m_tail, tm_tail, d_tail,l_tail,t_top,m_top,tm_top,d_top,l_top,totalcabinlength,V_tank_cyl, V_tank_tail, V_tank_top,tm_tanksystem,CGtank,CGfuel,CGcomb=cabin_design(0.5,0,30)


    

