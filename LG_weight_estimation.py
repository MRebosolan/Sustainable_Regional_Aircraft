


#Parameters for LG weight
#Gear compartment.........Ag.....Bg.....Cg.....Dg
#Main.....................40.0...0.16...0.019..0.000015
#Nose.....................20.0...0.10...0.0....0.000002
#Tail.....................5.0....0.0....0.0031.0.0

def LG_weight(Kgr, MTOW):
    #Already verified via hand calculation
    #Kgr= 1.0 for low wing AC, 1.08 for high wing AC, Wto=AC takeoff weight
    main_LG_weight= Kgr*(40 + 0.16*MTOW**0.75 + 0.019*MTOW + 0.000015*MTOW**1.5)
    nose_lg_weight= Kgr*(20 + 0.10*MTOW**0.75 + 0*MTOW + 0.000002*MTOW**1.5)
    return main_LG_weight, nose_lg_weight