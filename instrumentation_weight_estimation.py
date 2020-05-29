'''
input: MTOW, number of pilots, range, empty weight (i guess, not sure)
output: three ways of instrumentation weight calcs
responsible: jorn
'''
#instrumentation

def instrumentation_torenbeek(MTOW):
    N_e = 2
    return 120+ 20*N_e + 0.006*MTOW
    

def instrumentation_2 (W_e, R):
    return 0.575*(W_e**0.556)*(R**0.25)

def instrumentation_gd(MTOW, N_pil):
    N_e = 2
    flight = N_pil*(15+0.032*MTOW/1000)
    engine = N_e *(5+ 0.006*MTOW/1000)
    other = 0.15*(MTOW/1000)+0.012*MTOW
    return flight+engine+other