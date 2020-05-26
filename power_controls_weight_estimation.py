#power controls

def engine_controls(lf, b):
    N_e = 2 #number of engines
    
    return 88.46*((lf + b) * N_e/100)**0.294

def starting(W_e, pneumatic = True):
    if pneumatic:
        return 9.33*(W_e/1000)**1.078
    else:
        return 38.93*(W_e/1000)**0.918
    
def total(lf, b, W_e, pneumatic = True):
    return engine_controls(lf,b), starting(W_e, pneumatic)