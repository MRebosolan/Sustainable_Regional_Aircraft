# -*- coding: utf-8 -*-
"""
Created on Wed Jun  3 13:43:05 2020

@author: malfl
"""

import hydrogen_tank_sizing
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
paxsectionlength=round(15*seat_pitch+emergency_clearance,3)
aft_galley_length=1.5
front_galley_length=2
totalcabinlength=paxsectionlength+aft_galley_length+front_galley_length

outer_diameter=3.486
inner_diameter=3.286
V_tank=30 #CLASS I
R_tank_fus=1.5

tankchoice=3
if tankchoice==1:

    #FIRST OPTION ---- completely aft
    tank_sizing_fuselage(V_tank,R_tank_fus,2)
        
    
elif tankchoice==2:
    #SECOND OPTION  ---- also in tail
    fractionintail=0.5
    V_tank_tail=V_tank*fractionintail
    V_tank_fus=V_tank-V_tank_tail
    R_tank_tail=1
    tanksectionlength=V_tank_fus/3.14159/R_tank_fus**2
    tanksectionlengthtail=V_tank_tail/3.14159/R_tank_tail**2
    print(tanksectionlengh,tanksectionlengthtail)
    
    
    
elif tankchoice==3:
    #THIRD OPTION --- on top
    R_tank_top=tank_sizing(V_tank,totalcabinlength,2)[-1]
    print(R_tank_top)

else:
    #FOURTH OPTION --- combination
    fractioninfus=0.5
    V_tank_fus=fractioninfus*V_tank
    V_tank_top=V_tank-V_tank_fus
    R_tank_top=tank_sizing(V_tank_top,totalcabinlength,1)[-1]
    R_tank_fus=1.5
    tanksectionlength=V_tank_fus/3.14159/R_tank_fus**2

    print(R_tank_top,tanksectionlength)
    
    

