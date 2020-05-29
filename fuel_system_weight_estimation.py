#fuel system
'''
class 2
input: number of tanks, specific fuel weight, fuel weight
output: fuel system weight
responsible: jorn
'''

def W_fuelsystem (N_t, K_fsp, W_f):
    N_e = 2 #number of engines
    # N_t = 4 # number of fuel tanks
    #K_fsp = density in lbs/gal
    Wfs = 80 * (N_e + N_t -1) + 15*(N_t**0.5)*((W_f/K_fsp)**0.333)
    return Wfs

