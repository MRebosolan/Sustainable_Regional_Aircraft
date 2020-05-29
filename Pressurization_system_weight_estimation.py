"""responsible: jorn"""
def pressure_system_weight(lpax):
    #lpax=length of passenger cabin in ft
    Wapi= 6.75 * lpax**1.28
    return Wapi