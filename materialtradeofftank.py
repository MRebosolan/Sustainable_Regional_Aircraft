# -*- coding: utf-8 -*-
"""
Created on Thu Jun 11 15:19:33 2020

@author: malfl
"""

from hydrogen_tank_sizing import tank_sizing



    


##########ALUMINIUM############
print('#########ALUMINIUM##########')

massbest=99999999

for TS0 in range(155,165,2):
    
    MAX_ALLOWABLE_STRESS = 414000000  # Max allowable stress of the Aluminium alloy 2014 http://asm.matweb.com/search/SpecificMaterial.asp?bassnum=MA2014T6    
    EMOD=72.4*10**9
    THERMALEXP=2.25*10**(-5)
    TANK_MATERIAL_DENSITY = 2800 #MONOLITHIC METAL Aluminium alloy 2219 KG/M3  
    
    TANK_THICKNESS,(STRUCTURAL_TANK_MASS),(TOTAL_STRUCTURAL_TANK_MASS), (TANK_DIAMETER),(LENGTH),EXCEED=tank_sizing(12.5,8,2,TS0,TANK_MATERIAL_DENSITY,EMOD,THERMALEXP,MAX_ALLOWABLE_STRESS)
    TANK_THICKNESSf,(STRUCTURAL_TANK_MASSf),(TOTAL_STRUCTURAL_TANK_MASSf), (TANK_DIAMETERf),(LENGTHf),EXCEEDf=tank_sizing_fuselage(12.5,1.6,1,TS0,TANK_MATERIAL_DENSITY,EMOD,THERMALEXP,MAX_ALLOWABLE_STRESS)
    print('mass and TS0 ',TOTAL_STRUCTURAL_TANK_MASSf+TOTAL_STRUCTURAL_TANK_MASS,TS0 )
    if TOTAL_STRUCTURAL_TANK_MASSf+TOTAL_STRUCTURAL_TANK_MASS<massbest and not EXCEED and not EXCEEDf:
        TS0best=TS0
        massbest=TOTAL_STRUCTURAL_TANK_MASSf+TOTAL_STRUCTURAL_TANK_MASS
        print(massbest,TS0,TANK_THICKNESS)
    
print()
print('t pod: ',TANK_THICKNESS,', mass of both pods: ',TOTAL_STRUCTURAL_TANK_MASS)
print('t fus: ',TANK_THICKNESSf,', mass of fus tank: ',TOTAL_STRUCTURAL_TANK_MASSf)
print('TOTAL: ',TOTAL_STRUCTURAL_TANK_MASSf+TOTAL_STRUCTURAL_TANK_MASS)
print('####################################')

################################################

massbest=99999999
#TITANIUM http://asm.matweb.com/search/SpecificMaterial.asp?bassnum=MTA520
print('#########TITANIUM##########')
for TS0 in range(0,310,5):
    
    EMOD=120*10**9
    THERMALEXP=9.4*10**(-6)
    MAX_ALLOWABLE_STRESS=827000000
    TANK_MATERIAL_DENSITY = 4480
    
    TANK_THICKNESS,(STRUCTURAL_TANK_MASS),(TOTAL_STRUCTURAL_TANK_MASS), (TANK_DIAMETER),(LENGTH),EXCEED=tank_sizing(12.5,8,2,TS0,TANK_MATERIAL_DENSITY,EMOD,THERMALEXP,MAX_ALLOWABLE_STRESS)
    TANK_THICKNESSf,(STRUCTURAL_TANK_MASSf),(TOTAL_STRUCTURAL_TANK_MASSf), (TANK_DIAMETERf),(LENGTHf),EXCEEDf=tank_sizing_fuselage(12.5,1.6,1,TS0,TANK_MATERIAL_DENSITY,EMOD,THERMALEXP,MAX_ALLOWABLE_STRESS)
    print('mass and TS0 ',TOTAL_STRUCTURAL_TANK_MASSf+TOTAL_STRUCTURAL_TANK_MASS,TS0 )
    if TOTAL_STRUCTURAL_TANK_MASSf+TOTAL_STRUCTURAL_TANK_MASS<massbest and not EXCEED and not EXCEEDf:
        TS0best=TS0
        massbest=TOTAL_STRUCTURAL_TANK_MASSf+TOTAL_STRUCTURAL_TANK_MASS
        print(massbest,TS0,TANK_THICKNESS)
print()
print('t pod: ',TANK_THICKNESS,', mass of both pods: ',TOTAL_STRUCTURAL_TANK_MASS)
print('t fus: ',TANK_THICKNESSf,', mass of fus tank: ',TOTAL_STRUCTURAL_TANK_MASSf)
print('TOTAL: ',TOTAL_STRUCTURAL_TANK_MASSf+TOTAL_STRUCTURAL_TANK_MASS)
print('####################################')




print('#########CFRP##########')
for TS0 in range(0,310,5):

    #COMPOSITE https://www.hexcel.com/user_area/content_media/raw/IM7_HexTow_DataSheet.pdf
    MAX_ALLOWABLE_STRESS=5516000000
    THERMALEXP=-0.64*10**(-6)
    EMOD=276*10**9
    TANK_MATERIAL_DENSITY =1780   
    
    TANK_THICKNESS,(STRUCTURAL_TANK_MASS),(TOTAL_STRUCTURAL_TANK_MASS), (TANK_DIAMETER),(LENGTH),EXCEED=tank_sizing(12.5,8,2,TS0,TANK_MATERIAL_DENSITY,EMOD,THERMALEXP,MAX_ALLOWABLE_STRESS)
    TANK_THICKNESSf,(STRUCTURAL_TANK_MASSf),(TOTAL_STRUCTURAL_TANK_MASSf), (TANK_DIAMETERf),(LENGTHf),EXCEEDf=tank_sizing_fuselage(12.5,1.6,1,TS0,TANK_MATERIAL_DENSITY,EMOD,THERMALEXP,MAX_ALLOWABLE_STRESS)
    print('mass and TS0 ',TOTAL_STRUCTURAL_TANK_MASSf+TOTAL_STRUCTURAL_TANK_MASS,TS0 )
    if TOTAL_STRUCTURAL_TANK_MASSf+TOTAL_STRUCTURAL_TANK_MASS<massbest and not EXCEED and not EXCEEDf:
        TS0best=TS0
        massbest=TOTAL_STRUCTURAL_TANK_MASSf+TOTAL_STRUCTURAL_TANK_MASS
        print(massbest,TS0,TANK_THICKNESS)
        
print()
print('t pod: ',TANK_THICKNESS,', mass of both pods: ',TOTAL_STRUCTURAL_TANK_MASS)
print('t fus: ',TANK_THICKNESSf,', mass of fus tank: ',TOTAL_STRUCTURAL_TANK_MASSf)
print('TOTAL: ',TOTAL_STRUCTURAL_TANK_MASSf+TOTAL_STRUCTURAL_TANK_MASS)
print('####################################')

#############################################################################