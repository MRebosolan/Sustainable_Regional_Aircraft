# -*- coding: utf-8 -*-
"""
Created on Fri May 29 09:02:41 2020

@author: malfl
"""


import input 
import matplotlib.pyplot as plt
def tank_sizing(HYDROGENVOLUME,LENGTH,N):
    #N is number of tanks
    R=0
    result=0
    hsc=R/2 #spherical cap height
    asc=R #spherical cap radius
    while result<=HYDROGENVOLUME:
        R+=0.001
        hsc=R/2 #spherical cap height
        asc=R #spherical cap radius
        result=3.14159*(LENGTH-2*hsc*N)*R**2+N*2*3.14159*hsc/6*(3*asc**2+hsc**2) #NOW INCLUDES SPHERICAL CAPS INSTEAD OF HEMISPHERICAL CAPS!
    # print(R)
    TANK_DIAMETER=R*2

    TANK_SURFACE_AREA = (LENGTH-N*2*hsc)/N*3.14159*TANK_DIAMETER + 4*3.14159*R*hsc #ONE TANK!!!


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
    T_SURROUND= 295#K
    BOIL_OFF=0.05*HYDROGENVOLUME*DENSITY_LH/(3600*4)
    
    
    GAS_DIFFUSIVITY = -3.119*10**(-6)+3.541*10**(-8)*T_SURROUND+1.679*10**(-10)*T_SURROUND**2
    GAS_VISCOSITY = -2.079*10**(-6)+2.777*10**(-8)*T_SURROUND+1.077*10**(-10)*T_SURROUND**2
    PR = 0.71

    K_G = 0.02426 #.02041 W / m K
    BOLTZMANN_CONSTANT = 5.67*10**(-8) # units W/m^2 K^4
    EMIS_INSULATION = 0.9 # emissivity of foam
    THERMAL_CONDUCTIVITY_INSULATION = 0.0096  #VARRIES WITH TEMPERATURE
    INSULATION_DENSITY=35.3
    T_LH2 = 10 #TEMPERATURE OF THE CRYOGENIC TANK
    # print(R)
    xilist=[]
    yilist=[]
    
    IWISHTOLOOP=False #SET THIS VALUE TO TRUE IF YOU WANT TO DETERMINE THICKNESS ACCURATELY, OTHERWISE THREE CM IS USED WHICH SHOULD BE SUFFICENT
    
    if IWISHTOLOOP:
        for L in range(1,100):
            INSULATION_THICKNESS=L/1000 #unit is meters, increment from 0,001 to 0,1 meters by a millimeter each time
            running=True
            T_SURROUND= 290#K
            T_INITINSULATION =0#K
            
            while running:
                
                R_AD= 9.81*(1/T_SURROUND)*(T_SURROUND-T_INITINSULATION)*TANK_DIAMETER**3/(GAS_DIFFUSIVITY*GAS_VISCOSITY) #1/T_SURROUND IN KELVIN
                N_UD = (0.6+0.387*R_AD**(1/6)/(1+(0.559/PR)**(9/16))**(8/27))**2
                try:
                    CC_TANKAIR = N_UD*K_G/TANK_DIAMETER #convection coefficient
                except:
                    CC_TANKAIR = 0
                Q_CONVECTION = CC_TANKAIR*(T_SURROUND-T_INITINSULATION)
                Q_RADIATION = EMIS_INSULATION*BOLTZMANN_CONSTANT*(T_SURROUND**4-T_INITINSULATION**4)
                Q_IN = Q_CONVECTION + Q_RADIATION
    
                Q_CONDUCTION = THERMAL_CONDUCTIVITY_INSULATION*(T_INITINSULATION-T_LH2)/(INSULATION_THICKNESS)
                T_INITINSULATION +=0.1 #look for correct insulation temperature for each thickness
                if Q_CONDUCTION>Q_IN:
                    running=False
            
            BOIL_OFF=N*THERMAL_CONDUCTIVITY_INSULATION*TANK_SURFACE_AREA/INSULATION_THICKNESS*(T_INITINSULATION-T_SURROUND)/446592 #BOIL OFF for thickness and insulation T.
            # print(BOIL_OFF,BOIL_OFF*4*3600,INSULATION_THICKNESS,T_INITINSULATION,TANK_SURFACE_AREA)
            xilist.append(INSULATION_THICKNESS)
            yilist.append(-BOIL_OFF*4*3600)
        #print(xilist)#Prints out insulation thicknesses and boil offs
        #print(yilist)#Prints out insulation thicknesses and boil offs
        counter=0
        for boiloff in yilist:
            
            if boiloff<0.005*HYDROGENVOLUME*DENSITY_LH:
                #print('boiloff ',boiloff)
                INSULATION_THICKNESS=xilist[counter]
                #print('insulation thickness ',INSULATION_THICKNESS)
                break
            
            counter+=1
    else:
        INSULATION_THICKNESS=0.03
            
    INSULATION_MASS=INSULATION_THICKNESS*TANK_SURFACE_AREA*N*INSULATION_DENSITY
    STRUCTURAL_TANK_MASS+=INSULATION_MASS
    TOTAL_STRUCTURAL_TANK_MASS =N*STRUCTURAL_TANK_MASS    
    TANK_DIAMETER+=2*INSULATION_THICKNESS+2*TANK_THICKNESS
    #print(HYDROGENVOLUME)
    
    return(TANK_THICKNESS,round(STRUCTURAL_TANK_MASS,2),round(TOTAL_STRUCTURAL_TANK_MASS,1), round(TANK_DIAMETER,3),round(LENGTH,3))
    
    
def tank_sizing_fuselage(HYDROGENVOLUME, R,N):
    #N is number of tanks
    LENGTH=0
    result=0
    hsc=R/5 #SPHERICAL CAP HEIGHT
    asc=R #SPHERICAL CAP RADIUS
    while result<HYDROGENVOLUME:
        LENGTH+=0.001
        result=3.14159*(LENGTH-2*hsc*N)*R**2+N*2*3.14159*hsc/6*(3*asc**2+hsc**2) #NOW INCLUDES SPHERICAL CAPS INSTEAD OF HEMISPHERICAL CAPS!
    # print(R)
    TANK_DIAMETER=R*2
    TANK_SURFACE_AREA =(LENGTH-N*2*hsc)/N*3.14159*TANK_DIAMETER + 4*3.14159*R*hsc#ONE TANK!!!
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
    T_SURROUND= 295#K
    BOIL_OFF=0.05*HYDROGENVOLUME*DENSITY_LH/(3600*4)
    
    
    GAS_DIFFUSIVITY = -3.119*10**(-6)+3.541*10**(-8)*T_SURROUND+1.679*10**(-10)*T_SURROUND**2
    GAS_VISCOSITY = -2.079*10**(-6)+2.777*10**(-8)*T_SURROUND+1.077*10**(-10)*T_SURROUND**2
    PR = 0.71

    K_G = 0.02426 #.02041 W / m K
    BOLTZMANN_CONSTANT = 5.67*10**(-8) # units W/m^2 K^4
    EMIS_INSULATION = 0.9 # emissivity of foam
    THERMAL_CONDUCTIVITY_INSULATION = 0.0096  #VARRIES WITH TEMPERATURE
    INSULATION_DENSITY=35.3
    T_LH2 = 10 #TEMPERATURE OF THE CRYOGENIC TANK
    # print(R)
    xilist=[]
    yilist=[]
    for L in range(1,100):
        INSULATION_THICKNESS=L/1000 #unit is meters, increment from 0,001 to 0,1 meters by a millimeter each time
        running=True
        T_SURROUND= 290#K
        T_INITINSULATION =0#K
        
        while running:
            
            R_AD= 9.81*(1/T_SURROUND)*(T_SURROUND-T_INITINSULATION)*TANK_DIAMETER**3/(GAS_DIFFUSIVITY*GAS_VISCOSITY) #1/T_SURROUND IN KELVIN
            N_UD = (0.6+0.387*R_AD**(1/6)/(1+(0.559/PR)**(9/16))**(8/27))**2
            try:
                CC_TANKAIR = N_UD*K_G/TANK_DIAMETER #convection coefficient
            except:
                CC_TANKAIR = 0
            Q_CONVECTION = CC_TANKAIR*(T_SURROUND-T_INITINSULATION)
            Q_RADIATION = EMIS_INSULATION*BOLTZMANN_CONSTANT*(T_SURROUND**4-T_INITINSULATION**4)
            Q_IN = Q_CONVECTION + Q_RADIATION

            Q_CONDUCTION = THERMAL_CONDUCTIVITY_INSULATION*(T_INITINSULATION-T_LH2)/(INSULATION_THICKNESS)
            T_INITINSULATION +=0.1 #look for correct insulation temperature for each thickness
            if Q_CONDUCTION>Q_IN:
                running=False
        
        BOIL_OFF=N*THERMAL_CONDUCTIVITY_INSULATION*TANK_SURFACE_AREA/INSULATION_THICKNESS*(T_INITINSULATION-T_SURROUND)/446592 #BOIL OFF for thickness and insulation T.
        # print(BOIL_OFF,BOIL_OFF*4*3600,INSULATION_THICKNESS,T_INITINSULATION,TANK_SURFACE_AREA)
        xilist.append(INSULATION_THICKNESS)
        yilist.append(-BOIL_OFF*4*3600)
    #print(xilist)#Prints out insulation thicknesses and boil offs
    #print(yilist)#Prints out insulation thicknesses and boil offs
    counter=0
    for boiloff in yilist:
        
        if boiloff<0.005*HYDROGENVOLUME*DENSITY_LH:
            #print('boiloff ',boiloff)
            INSULATION_THICKNESS=xilist[counter]
            #print('insulation thickness ',INSULATION_THICKNESS)
            break
        
        counter+=1

    INSULATION_MASS=INSULATION_THICKNESS*TANK_SURFACE_AREA*N*INSULATION_DENSITY
        
    STRUCTURAL_TANK_MASS+=INSULATION_MASS
    TOTAL_STRUCTURAL_TANK_MASS =N*STRUCTURAL_TANK_MASS 
    TANK_DIAMETER+=2*INSULATION_THICKNESS+2*TANK_THICKNESS
        
    return(TANK_THICKNESS,round(STRUCTURAL_TANK_MASS,2),round(TOTAL_STRUCTURAL_TANK_MASS,1), round(TANK_DIAMETER,3),round(LENGTH,3))
    
    

