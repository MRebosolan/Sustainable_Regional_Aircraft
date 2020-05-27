import numpy as np
import matplotlib.pyplot as plt

###########################################################
# Parameters

MTOW = 36000#kg
V_S, V_C, V_D, V_B, V_A, V_H = 0, 0, 0, 0, 0, 0
V = 0 #Quadratic variable
G_W, rho, S, CLmax, CLa, CNmax = 91.391484, 0.00237, 871.88, 1.8, 0, 0 
kc = 33
nlimpos, nlimneg = 0, 0
nlim = 4.4 #usually 4.4 but check 
Kg = 0 #Gust alleviation factor
Ude = 0 #depends on altitude in ft (above or below 20k)
mug = 0 
cbar = 0

###########################################################
#Conversions to Imperial 
MTOW = 2.20462*MTOW 
G_W = 2.20462*G_W 
S = 10.7639*S

##########################################################
#Calculations
# Gross Weight (G_W)
# Same as ramp weight but also considering taxiing fuel

G_W=MTOW/0.99
#S=MTOW*10.7639/(4375.84)

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

V_S=(np.sqrt(2*(91.391484)/(rho*CNmax)))*0.592484
#V_S=137#kts

# Determination of Design Cruise Speed (V_C minimum)

V_C=np.sqrt(91.391484)*kc
#V_C=186.91#kts
    
V_B=V_C-43#kts

V_H=175#kts

# Determination of Design Manoeuvering Speed (V_A minimum)

V_A=np.sqrt(nlimpos)*V_S
#V_A=425#kts

#Determination of Design Diving Speed (V_D minimum)

V_D=1.25*V_C
#V_D=369#kts    

#Construction of Gust Load Factor Lines

def mug(G_W,S,rho,cbar,CLa):
    return 2*(G_W/S)/(rho*cbar*CLa)

def Kg(mug):
    return 0.88*mug/(5.3+mug)

def nlim(Kg,Ude,V,Cla,G_W,S):
    return 1 + (Kg*Ude*V*CLa)/(498*(G_W/S))

################################################################
# Plotting

# Quadratics
quadlist1=[0,V_S,V_A]
quadlist2=[0,1,nlimpos]
quadlist3=[0,V_S,V_H]
quadlist4=[0,-0.6,-1]

quad = np.polyfit(quadlist1, quadlist2, 2)
quad2 = np.polyfit(quadlist3, quadlist4, 2)
f = np.poly1d(quad)
f2 = np.poly1d(quad2)

x_new = np.linspace(quadlist1[0], quadlist1[-1], 50)
x_new2 = np.linspace(quadlist3[0], quadlist3[-1], 50)
y_new = f(x_new)
y_new2 = f2(x_new2)


# Linears
constlist= np.arange(0,300)
constlist2= np.arange(V_A,V_D)
lits=np.arange(0,V_C)

plt.plot(constlist,nlimpos*np.ones(300),'--',color = 'r')
plt.plot(np.arange(0,V_C),np.linspace(1,nlimpos,187),'--',color = 'g')
plt.plot(np.arange(0,V_C),np.linspace(1,-1,187),'--',color = 'g')
plt.plot(np.arange(0,V_D),np.linspace(1,nlimpos,369),'--',color = 'g')
plt.plot(np.arange(0,V_D),np.linspace(1,-1,369),'--',color = 'g')
plt.plot(V_C*np.ones(100),np.linspace(nlimpos,0,100),'--',color = 'g')
plt.plot(V_A*np.ones(100),np.linspace(nlimpos,0,100),'--',color = 'g')
plt.plot(V_H*np.ones(100),np.linspace(0,-1,100),'--',color = 'g')
plt.plot(constlist,-1*np.ones(300),'--', color = 'r')
plt.plot(constlist2,nlimpos*np.ones(V_D-V_A),'b')
plt.plot(V_D*np.ones(100),np.linspace(nlimpos,0,100),'b')
plt.plot(np.arange(V_C,V_D),np.arange(-1,0,1/74),'b')
plt.plot(np.linspace(V_H,V_C,100),-1*np.ones(100),'b')
plt.plot(x_new,y_new,'b')
plt.plot(x_new2,y_new2,'b')
plt.text(V_A,nlimpos,'A')
plt.text(V_H,-1,'H')
plt.text(V_C,-1,'F')
plt.text(V_D,0,'E')
plt.text(V_D,nlimpos,'D')
plt.text(V_C,nlimpos,'C')

plt.ylim(-3,5)


plt.axhline(y=0, color='k')
plt.axvline(x=0, color='k')



plt.show()




