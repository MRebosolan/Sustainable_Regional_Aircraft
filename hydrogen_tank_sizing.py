# -*- coding: utf-8 -*-
"""
Created on Fri May 29 09:02:41 2020

@author: malfl
"""


import input 

def tank_sizing(HYDROGENVOLUME, CABIN_LENGTH,N):
    #N is number of tanks
    R=0
    result=0
    while result<HYDROGENVOLUME:
        R+=0.001
        result=N*4/3*3.14159*R**3+3.14159*(CABIN_LENGTH-N*2*R)*R**2

    TANK_DIAMETER=R*2
    TANK_SURFACE_AREA = (CABIN_LENGTH-N*2*R)/N*3.14159*TANK_DIAMETER + 3.14159*TANK_DIAMETER**2 #ONE TANK!!!
    TANK_MATERIAL_DENSITY = 2825 #MONOLITHIC METAL Aluminium alloy 2219 KG/M3

    # TANK THICKNESS COMPUTATION

    DENSITY_LH = input.rho_hydrogen  # Desnity liquid hydrogen 71 KG/M3
    R_GCH = 4157  # Gas constant of liquid hydrogen 4157 Nm/kg K
    T_CT = 13.15  # Temperature is around -260C(13.15K) for cryogenic tank

    PRESSURE_GAS = 145000 # Mantained at 1.45*10^5 pa to minimize boil off

    # CYLINDRICAL TANK WITH HEMISPHERICAL ENDS CAPS

    FOS_TANK = 2  # Safety factor for the tank
    MAX_ALLOWABLE_STRESS = 219000000  # Max allowable stress of the Aluminium alloy 2219 tank: 219 MPa

    TANK_THICKNESS = PRESSURE_GAS * R * FOS_TANK / (2 * MAX_ALLOWABLE_STRESS)
    # print(TANK_THICKNESS,PRESSURE_GAS,R)

    STRUCTURAL_TANK_MASS = TANK_SURFACE_AREA*TANK_MATERIAL_DENSITY*TANK_THICKNESS #ONE TANK!!!
    TOTAL_STRUCTURAL_TANK_MASS =N*STRUCTURAL_TANK_MASS


    #INSULATION
    T_SURROUND= 1
    T_INITINSULATION =1
    GAS_DIFFUSIVITY = -3.119*10**(-6)+3.541*10**(-8)*T_SURROUND+1.679*10**(-10)*T_SURROUND**2
    GAS_VISCOTY = -2.079*10**(-6)+2.777*10**(-8)*T_SURROUND+1.077*10**(-10)*T_SURROUND**2
    R_AD= 9.81*(1/T_SURROUND)*(T_SURROUND-T_INITINSULATION)*TANK_DIAMETER**3/(GAS_DIFFUSIVITY*GAS_VISCOTY)
    UK_P=1
    UK_GC=1
    K_G=1
    N_UD = (0.6+0.387*R_AD**(1/6)/(1+(0.559/(UK_P*UK_GC))**(9/16))**(8/27))**2

    CC_TANKAIR = N_UD*K_G/TANK_DIAMETER

    #Q_CONVECTION = CC_TANKAIR*(T_SURROUND-T_INITINSULATION)

    return(TANK_THICKNESS,STRUCTURAL_TANK_MASS, TOTAL_STRUCTURAL_TANK_MASS, TANK_DIAMETER)
