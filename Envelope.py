import numpy as np
import matplotlib.pyplot as plt

###########################################################
# Parameters

MTOW = 36000#kg
V_S, V_C, V_D, V_B, V_A = 0, 0, 0, 0, 0
V = 0 #Quadratic variable
G_W, rho, S, CLmax, CLa, CNmax = 0, 0.00080329261, 88.04879, 1.4, 0, 0 
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

G_W=MTOW/0.99

# Determination of Design Limit Load Factor (nlim pos min)

nlimpos=2.1+24000/(G_W+10000)
#        if 2.1+24000/(G_W+10000)<(3.8*MTOW):
#            if 2.1+24000/(G_W+10000)>2.5:
#                return 
#nlimpos<3.8 @W_TO
#nlimpos>2.5 always
#nlim=4.4

nlimneg=-0.4*nlim

# Preiliminary Design Assumption

CNmax=1.1*CLmax

# Determination of Stall Speed

#V_S=np.sqrt((2(G_W/S)/rho*CNmax))
V_S=130#kts

# Determination of Design Cruise Speed (V_C minimum)

V_C=np.sqrt(kc*(G_W/S))
V_C=200#kts
    
V_B=V_C-43#kts
V_B=190#kts

# Determination of Design Manoeuvering Speed (V_A minimum)

V_A=np.sqrt(V_S*nlim)
V_A=160#kts

#Determination of Design Diving Speed (V_D minimum)

V_D=1.25*V_C
V_D=240#kts    

#Construction of Gust Load Factor Lines

def mug(G_W,S,rho,cbar,CLa):
    return 2*(G_W/S)/(rho*cbar*CLa)

def Kg(mug):
    return 0.88*mug/(5.3+mug)

def nlim(Kg,Ude,V,Cla,G_W,S):
    return 1 + (Kg*Ude*V*CLa)/(498*(G_W/S))

################################################################
# Plotting

quadlist1=[0,1,nlimpos]
quadlist2=[0,V_S,V_A]
quad=np.polyfit(quadlist1,quadlist2,2)    
constlist= np.arange(0,300)
constlist4= np.arange(V_A,V_D)


plt.plot(constlist,nlimpos*np.ones(300),'--')
plt.plot(constlist,nlimneg*np.ones(300),'--')
plt.plot(constlist4,nlimpos*np.ones(80))
plt.plot(V_D*np.ones(4),np.arange(nlimpos,nlimneg+1,-1.0))
plt.plot(np.arange(V_C,V_D),np.arange(nlimneg,nlimneg+1,1/40))

plt.axhline(y=0, color='k')
plt.axvline(x=0, color='k')



plt.show()




