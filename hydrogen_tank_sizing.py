# -*- coding: utf-8 -*-
"""
Created on Fri May 29 09:02:41 2020

@author: malfl
"""


import input 
import matplotlib.pyplot as plt
import numpy as np
def tank_sizing(HYDROGENVOLUME,LENGTH,N,TS0=155,TANK_MATERIAL_DENSITY=1780,EMOD=276*10**9,THERMALEXP=-0.64*10**(-6),MAX_ALLOWABLE_STRESS = 5516000000):
    #N is number of tanks
    R=0
    result=0
    hsc=R/2 #spherical cap height
    asc=R #spherical cap radius
    while result<=HYDROGENVOLUME:
        R+=0.001
        hsc=R/2 #spherical cap height
        asc=R #spherical cap radius

        if (LENGTH-2*hsc*N)>0:
            result=np.pi*(LENGTH-2*hsc*N)*R**2+N*2*np.pi*hsc/6*(3*asc**2+hsc**2) #NOW INCLUDES SPHERICAL CAPS INSTEAD OF HEMISPHERICAL CAPS!
        else:
            result=N*2*3.14159*hsc/6*(3*asc**2+hsc**2)
            
        #print(R,result)
    TANK_DIAMETER=R*2

    TANK_SURFACE_AREA = (LENGTH-N*2*hsc)/N*3.14159*TANK_DIAMETER + 2*3.14159*R*hsc #ONE TANK!!!




    # TANK THICKNESS COMPUTATION

    DENSITY_LH = input.rho_hydrogen  # Desnity liquid hydrogen 71 KG/M3
    R_GCH = 4157  # Gas constant of liquid hydrogen 4157 Nm/kg K
    T_CT = 20  #  cryogenic tank

    PRESSURE_GAS = 145000 # Mantained at 1.45*10^5 pa to minimize boil off

    # CYLINDRICAL TANK WITH HEMISPHERICAL ENDS CAPS

    FOS_TANK = 2  # Safety factor for the tank
    
    extra_T=10
    TDIFFERENCE=T_CT-TS0
    thermal_stress_c=EMOD*THERMALEXP*(TDIFFERENCE)
    thermal_stress_h=EMOD*THERMALEXP*(T_CT+extra_T-TS0)
    thermal_stress_fillup=EMOD*THERMALEXP*(310-TS0)
    EXCEED=False
    #print('Thermal stress fuselage tank: ',thermal_stress)
    
    ###DESIGN CASES
    TANK_THICKNESS_C = PRESSURE_GAS *R/(abs(MAX_ALLOWABLE_STRESS/FOS_TANK-thermal_stress_c)) #PRESSURIZED AND UNDER THERMAL STRESS
    TANK_THICKNESS_H = PRESSURE_GAS *R/(abs(MAX_ALLOWABLE_STRESS/FOS_TANK-thermal_stress_h))
    TANK_THICKNESS_FILLUP = PRESSURE_GAS *R/(abs(MAX_ALLOWABLE_STRESS/FOS_TANK-thermal_stress_fillup))   
    

    TANK_THICKNESS=max(TANK_THICKNESS_FILLUP,TANK_THICKNESS_H,TANK_THICKNESS_C)
    
    if EMOD*THERMALEXP*(max((TS0-20),(310-TS0)))>MAX_ALLOWABLE_STRESS/FOS_TANK:
        print('exceed',MAX_ALLOWABLE_STRESS,EMOD*THERMALEXP*(max((TS0-20),(310-TS0))),EMOD,THERMALEXP)
        EXCEED=True

        
    # print(TANK_THICKNESS,PRESSURE_GAS,R)

    STRUCTURAL_TANK_MASS = TANK_SURFACE_AREA*TANK_MATERIAL_DENSITY*TANK_THICKNESS #ONE TANK!!!
    TOTAL_STRUCTURAL_TANK_MASS =N*STRUCTURAL_TANK_MASS


    #INSULATION
    T_SURROUND= 310#K
    BOIL_OFF=0.05*HYDROGENVOLUME*DENSITY_LH/(3600*4)
    
    
    GAS_DIFFUSIVITY = -3.119*10**(-6)+3.541*10**(-8)*T_SURROUND+1.679*10**(-10)*T_SURROUND**2
    GAS_VISCOSITY = -2.079*10**(-6)+2.777*10**(-8)*T_SURROUND+1.077*10**(-10)*T_SURROUND**2
    PR = 0.71

    K_G = 0.02426 #.02041 W / m K
    BOLTZMANN_CONSTANT = 5.67*10**(-8) # units W/m^2 K^4
    EMIS_INSULATION = 0.9 # emissivity of foam
    THERMAL_CONDUCTIVITY_INSULATION = 0.0096  #VARRIES WITH TEMPERATURE, CHECK NASA
    INSULATION_DENSITY=35.3
    T_LH2 = 10 #TEMPERATURE OF THE CRYOGENIC TANK
    # print(R)
    xilist=[]
    yilist=[]
    
    IWISHTOLOOP=False #SET THIS VALUE TO TRUE IF YOU WANT TO DETERMINE THICKNESS ACCURATELY, OTHERWISE THREE CM IS USED WHICH SHOULD BE SUFFICIENT
    
    if IWISHTOLOOP:
        for L in range(1,100):
            INSULATION_THICKNESS=L/1000 #unit is meters, increment from 0,001 to 0,1 meters by a millimeter each time
            running=True
            T_SURROUND= 310#K
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
            running=True
            T_SURROUND= 205#K
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
            
            BOIL_OFF_COLD=N*THERMAL_CONDUCTIVITY_INSULATION*TANK_SURFACE_AREA/INSULATION_THICKNESS*(T_INITINSULATION-205)/446592 #BOIL OFF for thickness and insulation T.
            # print(BOIL_OFF,BOIL_OFF*4*3600,INSULATION_THICKNESS,T_INITINSULATION,TANK_SURFACE_AREA)
            xilist.append(INSULATION_THICKNESS)
            yilist.append(-BOIL_OFF*6*3600-BOIL_OFF_COLD*4*3600)
        #print(xilist)#Prints out insulation thicknesses and boil offs
        #print(yilist)#Prints out insulation thicknesses and boil offs
        counter=0
        for boiloff in yilist:
            
            if boiloff<0.01*HYDROGENVOLUME/1.072*DENSITY_LH:
                #print('boiloff ',boiloff)
                INSULATION_THICKNESS=xilist[counter]
                print('insulation thickness ',INSULATION_THICKNESS)

                break
            
            counter+=1
    else:
        INSULATION_THICKNESS=0.034
            
    INSULATION_MASS=INSULATION_THICKNESS*TANK_SURFACE_AREA*N*INSULATION_DENSITY
    #print('INSULATION MASS POD: ', INSULATION_MASS)
    STRUCTURAL_TANK_MASS+=INSULATION_MASS
    TOTAL_STRUCTURAL_TANK_MASS =N*STRUCTURAL_TANK_MASS    
    TANK_DIAMETER+=2*INSULATION_THICKNESS+2*TANK_THICKNESS
    #print(HYDROGENVOLUME)
    
    return(TANK_THICKNESS,round(STRUCTURAL_TANK_MASS,2),round(TOTAL_STRUCTURAL_TANK_MASS,1), round(TANK_DIAMETER,3),round(LENGTH,3),EXCEED)
    
    
    
    
    
##############################################################################################################################    
    
#####   ########################################################################################################################
    
    
    
    
    
def tank_sizing_fuselage(HYDROGENVOLUME, R,N,TS0=155,TANK_MATERIAL_DENSITY=1780,EMOD=276*10**9,THERMALEXP=-0.64*10**(-6),MAX_ALLOWABLE_STRESS = 5516000000):
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
    TANK_SURFACE_AREA =(LENGTH-N*2*hsc)/N*3.14159*TANK_DIAMETER + 2*3.14159*R*hsc#ONE TANK!!!

    # TANK THICKNESS COMPUTATION

    DENSITY_LH = input.rho_hydrogen  # Desnity liquid hydrogen 71 KG/M3
    R_GCH = 4157  # Gas constant of liquid hydrogen 4157 Nm/kg K
    T_CT = 20  #  cryogenic tank

    PRESSURE_GAS = 145000 # Mantained at 1.45*10^5 pa to minimize boil off

    # CYLINDRICAL TANK WITH HEMISPHERICAL ENDS CAPS

    FOS_TANK = 2  # Safety factor for the tank
    
    extra_T=10
    TDIFFERENCE=T_CT-TS0
    thermal_stress_c=EMOD*THERMALEXP*(TDIFFERENCE)
    thermal_stress_h=EMOD*THERMALEXP*(T_CT+extra_T-TS0)
    thermal_stress_fillup=EMOD*THERMALEXP*(310-TS0)
    EXCEED=False
    #print('Thermal stress fuselage tank: ',thermal_stress)
    
    ###DESIGN CASES
    TANK_THICKNESS_C = PRESSURE_GAS *R/(abs(MAX_ALLOWABLE_STRESS/FOS_TANK-thermal_stress_c)) #PRESSURIZED AND UNDER THERMAL STRESS
    TANK_THICKNESS_H = PRESSURE_GAS *R/(abs(MAX_ALLOWABLE_STRESS/FOS_TANK-thermal_stress_h))
    TANK_THICKNESS_FILLUP = PRESSURE_GAS *R/(abs(MAX_ALLOWABLE_STRESS/FOS_TANK-thermal_stress_fillup))   
    

    TANK_THICKNESS=max(TANK_THICKNESS_FILLUP,TANK_THICKNESS_H,TANK_THICKNESS_C)
    
    if EMOD*THERMALEXP*(max((TS0-20),(310-TS0)))>MAX_ALLOWABLE_STRESS/FOS_TANK:
        print('exceed',MAX_ALLOWABLE_STRESS,EMOD*THERMALEXP*(max((TS0-20),(310-TS0))),EMOD,THERMALEXP)
        EXCEED=True


    # print(TANK_THICKNESS,PRESSURE_GAS,R)

    STRUCTURAL_TANK_MASS = TANK_SURFACE_AREA*TANK_MATERIAL_DENSITY*TANK_THICKNESS #ONE TANK!!!
    TOTAL_STRUCTURAL_TANK_MASS =N*STRUCTURAL_TANK_MASS


    #INSULATION
    T_SURROUND= 310#K
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
    
    IWISHTOLOOP=True #SET THIS VALUE TO TRUE IF YOU WANT TO DETERMINE THICKNESS ACCURATELY, OTHERWISE THREE CM IS USED WHICH SHOULD BE SUFFICIENT
    
    if IWISHTOLOOP:
        for L in range(1,100):
            INSULATION_THICKNESS=L/1000 #unit is meters, increment from 0,001 to 0,1 meters by a millimeter each time
            running=True
            T_SURROUND= 310#K
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
            yilist.append(-BOIL_OFF*10*3600)
        #print(xilist)#Prints out insulation thicknesses and boil offs
        #print(yilist)#Prints out insulation thicknesses and boil offs
        counter=0
        for boiloff in yilist:
            
            if boiloff<0.01*HYDROGENVOLUME/1.072*DENSITY_LH:
                INSULATION_THICKNESS=xilist[counter]
                print(INSULATION_THICKNESS)
                break
            
            counter+=1
            
            
    else:
        INSULATION_THICKNESS=0.03
        
    INSULATION_MASS=INSULATION_THICKNESS*TANK_SURFACE_AREA*N*INSULATION_DENSITY
    #print('INSULATION MASS FUS: ', INSULATION_MASS)
        
    STRUCTURAL_TANK_MASS+=INSULATION_MASS
    TOTAL_STRUCTURAL_TANK_MASS =N*STRUCTURAL_TANK_MASS 
    TANK_DIAMETER+=2*INSULATION_THICKNESS+2*TANK_THICKNESS
        
    return(TANK_THICKNESS,round(STRUCTURAL_TANK_MASS,2),round(TOTAL_STRUCTURAL_TANK_MASS,1), round(TANK_DIAMETER,3),round(LENGTH,3),EXCEED)
    
    

