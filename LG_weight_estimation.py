

#Parameters for LG weight
#Gear compartment.........Ag.....Bg.....Cg.....Dg
#Main.....................40.0...0.16...0.019..0.000015
#Nose.....................20.0...0.10...0.0....0.000002
#Tail.....................5.0....0.0....0.0031.0.0

def LG_weight(Kgr, Wto, Ag, Bg, Cg, Dg):
    #Already verified via hand calculation
    #Kgr= 1.0 for low wing AC, 1.08 for high wing AC, Wto=AC takeoff weight
    LG_weight= Kgr*(Ag + Bg*Wto**0.75 + Cg*Wto + Dg*Wto**1.5)
    return LG_weight
