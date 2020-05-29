#nacelle weight
from math import *

def W_nacelle_gd (A_inlet, ln, p2): #already accounted for two engines, A_inlet is area in ft2, ln = nacelle length from inlet to compressor, p2 = max static pressure (between 15 and 50 psi)
    return 7.435*2*(sqrt(A_inlet)*ln*p2)**0.731

def W_nacelle_torenbeek(T_TO):
    return 0.065*T_TO