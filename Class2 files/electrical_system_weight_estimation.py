#electrical system

def electrical_gd(w_fs, w_iae): # weights are fuel system and instrumentation system
    return 1163*((w_fs+w_iae)/1000)**0.506

def electrical_torenbeek(V_pax): #V_pax = cabin volume in ft^3
    return 10.8 *(V_pax**0.7)*(1- 0.018* (V_pax**0.35))