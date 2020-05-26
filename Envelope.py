import numpy as np


###########################################################
# Parameters

MTOW = 36000#kg
V_S, V_C, V_D, V_B, V_A = 0, 0, 0, 0, 0
V = 0 #Quadratic variable
G_W, rho, CNmax, S, CLmax, CLa = 0, 0.00080329261, 2.1, 88.04879, 0, 0
kc = 33
nlimpos, nlimneg = 0, 0
nlim = 4.4 #usually 4.4 but check 
Kg = 0 #Gust alleviation factor
Ude = 0 #depends on altitude in ft (above or below 20k)
mug = 0 
cbar = 0

###########################################################
#Conversions
MTOW = 2.20462*MTOW
G_W = 2.20462*G_W

##########################################################
#Calculations
# Gross Weight (G_W)
# Same as ramp weight but also considering taxiing fuel

class Envelope:
    
    
    def G_W(MTOW):
        return MTOW/0.99
    
    # Determination of Design Limit Load Factor (nlim pos min)
    
    def nlimpos(G_W):
        if 2.1+24000/(G_W+10000)<(3.8*MTOW):
            if 2.1+24000/(G_W+10000)>2.5:
                return 2.1+24000/(G_W+10000)
    #nlimpos<3.8 @W_TO
    #nlimpos>2.5 always
    #nlim=4.4
    
    def nlimneg(nlim):
        return 0.4*nlim
    
    # Preiliminary Design Assumption
    
    def CNmax(CLmax):
        return 1.1*CLmax
    
    # Determination of Stall Speed
    
    def V_S(G_W,S,rho,CNmax):
        return np.sqrt((2(G_W/S)/rho*CNmax))
    
    # Determination of Design Cruise Speed (V_C minimum)
    
    def V_C(kc,G_W,S):
        return np.sqrt(kc*(G_W/S))
        
    def V_B(V_C):
        return V_C-43#kts
    
    #Construction of Gust Load Factor Lines
    
    def mug(G_W,S,rho,cbar,CLa):
        return 2*(G_W/S)/(rho*cbar*CLa)
    
    def Kg(mug):
        return 0.88*mug/(5.3+mug)
    
    def nlim(Kg,Ude,V,Cla,G_W,S):
        return 1 + (Kg*Ude*V*CLa)/(498*(G_W/S))
    
    # Determination of Design Manoeuvering Speed (V_A minimum)
    
    def V_A(V_S,nlim):
        return np.sqrt(V_S*nlim)
    
    #Determination of Design Diving Speed (V_D minimum)
    
    def V_D(V_C):
        return 1.25*V_C




