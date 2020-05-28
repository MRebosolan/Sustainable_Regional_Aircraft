import numpy as np
import matplotlib.pyplot as plt

###########################################################
# Parameters

MTOW     = 36000     #kg
G_WvoerS    = 91.391484 #psf
rho      = 0.00237   #slugs/ft^3
S        = 871.88    #ft^2
CLmax    = 1.8
CLmaxneg = -1.0
kc       = 33        #Varies 36-29 as W/S varies from 20-100 psf
nlim     = 4.4       #usually 4.4 but check 
Kg       = 0.781     #Gust alleviation factor
UdeB     = 85        #Depends on altitude in ft (above 20k ft) [ROSKAM V pg 38]
UdeC     = 65        #Depends on altitude in ft (above 20k ft) [ROSKAM V pg 38]
UdeD     = 35        #Depends on altitude in ft (above 20k ft) [ROSKAM V pg 38]
CLa      = 5.44      #rad^-1
mug      = 0 
cbar     = 0

###########################################################
#Conversions to Imperial 
MTOW = 2.20462*MTOW 
S = 10.7639*S

##########################################################
#Calculations
# Gross Weight (G_W)
# Same as ramp weight but also considering taxiing fuel

G_W=MTOW/0.99

# Determination of Design Limit Load Factor (nlim pos min)

if 2.1+24000/(G_W+10000)>2.5:
    nlimpos=2.1+24000/(G_W+10000)
else:
    nlimpos = 2.5
    
#nlimpos<3.8 @W_TO
#nlimpos>2.5 always
#nlim=4.4

nlimneg=-0.4*nlim

# Preiliminary Design Assumption

CNmax=1.1*CLmax
CNmaxneg=1.1*CLmaxneg

# Determination of Stall Speed

V_S =(np.sqrt(2*(G_WvoerS)/(rho*CNmax)))*0.592484
V_S2=(V_S)*0.9 


# Determination of Design Cruise Speed (V_C minimum)

V_B=161.6#kts

#V_B=196.33#kts

V_C=np.sqrt(G_WvoerS)*kc

#V_C=V_B+43#kts
#V_C=186.91#kts

V_H=156#kts

# Determination of Design Manoeuvering Speed (V_A minimum)

V_A=np.sqrt(nlimpos)*V_S
#V_A=425#kts

#Determination of Design Diving Speed (V_D minimum)

V_D=1.25*V_C
#V_D=369#kts    

#Construction of Gust Load Factor Lines

def mug(G_W,S,rho,cbar,CLa):
    return 2*(G_WvoerS)/(rho*cbar*CLa)


nlimgustBslope =  (Kg*UdeB*CLa)/(498*(G_WvoerS))
nlimgustCslope =  (Kg*UdeC*CLa)/(498*(G_WvoerS))
nlimgustDslope =  (Kg*UdeD*CLa)/(498*(G_WvoerS))

################################################################################################
# Plotting

#  ----------------Quadratics-------------------------

quadlist1=[0,V_S,V_A]
quadlist2=[0,1,nlimpos]
quadlist3=[0,V_S,V_H]
quadlist4=[0,-0.6,-1]
quadlist5=[0,V_S2,V_B]
quadlist6=[0,1,1+nlimgustBslope*V_B]

quad  = np.polyfit(quadlist1, quadlist2, 2)
quad2 = np.polyfit(quadlist3, quadlist4, 2)
quad3 = np.polyfit(quadlist5, quadlist6, 2)

f  = np.poly1d(quad)
f2 = np.poly1d(quad2)
f3 = np.poly1d(quad3)

x_new  = np.linspace(quadlist1[0], quadlist1[-1], 50)
x_new2 = np.linspace(quadlist3[0], quadlist3[-1], 50)
x_new3 = np.linspace(quadlist5[0], quadlist5[-1], 50)

y_new  = f(x_new)
y_new2 = f2(x_new2)
y_new3 = f3(x_new3)

# ---------------Linears-----------------------------

constlist= np.arange(0,300)
constlist2= np.arange(V_A,V_D)

# -------- Line Intersection Point--------------------

def line_intersection(line1, line2):
    xdiff = (line1[0][0] - line1[1][0], line2[0][0] - line2[1][0])
    ydiff = (line1[0][1] - line1[1][1], line2[0][1] - line2[1][1])

    def det(a, b):
        return a[0] * b[1] - a[1] * b[0]

    div = det(xdiff, ydiff)
    if div == 0:
       raise Exception('lines do not intersect')

    d = (det(*line1), det(*line2))
    x = det(d, xdiff) / div
    y = det(d, ydiff) / div
    return x, y


#----------------------------------------------------------------------------------------------

plt.figure(0)

plt.plot(constlist,nlimpos*np.ones(300),'--',color = 'r')
plt.plot(constlist,-1*np.ones(300),'--', color = 'r')

#plt.plot(np.linspace(0,V_C,300),np.linspace(1,nlimpos,300),'--',color = 'g')
#plt.plot(np.linspace(0,V_C,300),np.linspace(1,-1,300),'--',color = 'g')
#plt.plot(np.linspace(0,V_D,400),np.linspace(1,nlimpos,400),'--',color = 'g')
#plt.plot(np.linspace(0,V_D,400),np.linspace(1,0,400),'--',color = 'g')

plt.plot(np.linspace(0,V_B,300),np.linspace(1,1+nlimgustBslope*V_B,300),'--',color = 'y')
plt.plot(np.linspace(0,V_B,300),np.linspace(1,1-nlimgustBslope*V_B,300),'--',color = 'y')
plt.plot(np.linspace(0,V_C,400),np.linspace(1,1+nlimgustCslope*V_C,400),'--',color = 'y')
plt.plot(np.linspace(0,V_C,300),np.linspace(1,1-nlimgustCslope*V_C,300),'--',color = 'y')
plt.plot(np.linspace(0,V_D,400),np.linspace(1,1+nlimgustDslope*V_D,400),'--',color = 'y')
plt.plot(np.linspace(0,V_D,300),np.linspace(1,1-nlimgustDslope*V_D,300),'--',color = 'y')

plt.plot(V_C*np.ones(100),np.linspace(nlimpos,0,100),'--',color = 'g')
plt.plot(V_A*np.ones(100),np.linspace(nlimpos,0,100),'--',color = 'g')
plt.plot(V_H*np.ones(100),np.linspace(0,-1,100),'--',color = 'g')
plt.plot(np.linspace(V_A,V_D,200),nlimpos*np.ones(200),'b')
plt.plot(V_D*np.ones(100),np.linspace(nlimpos,0,100),'b')
plt.plot(np.linspace(V_C,V_D,100),np.arange(-1,0,1/100),'b')
plt.plot(np.linspace(V_H,V_C,100),-1*np.ones(100),'b')
plt.plot(x_new,y_new,'b')
plt.plot(x_new2,y_new2,'b')
plt.text(V_A,nlimpos,'A')
plt.text(V_H,-1,'H')
plt.text(V_C,-1,'F')
plt.text(V_D,0,'E')
plt.text(V_D,nlimpos,'D')
plt.text(V_C,nlimpos,'C')

plt.ylim(-2,4)


plt.axhline(y=0, color='k')
plt.axvline(x=0, color='k')

# ---------------------------------------------------------------------------------------------

plt.figure(1)

plt.plot(constlist,nlimpos*np.ones(300),'--',color = 'r')
plt.plot(constlist,-1*np.ones(300),'--', color = 'r')

plt.plot(np.linspace(0,V_B,300),np.linspace(1,1+nlimgustBslope*V_B,300),'--',color = 'r')
plt.plot(np.linspace(0,V_B,300),np.linspace(1,1-nlimgustBslope*V_B,300),'--',color = 'r')
plt.plot(np.linspace(0,V_C,400),np.linspace(1,1+nlimgustCslope*V_C,400),'--',color = 'r')
plt.plot(np.linspace(0,V_C,300),np.linspace(1,1-nlimgustCslope*V_C,300),'--',color = 'r')
plt.plot(np.linspace(0,V_D,400),np.linspace(1,1+nlimgustDslope*V_D,400),'--',color = 'r')
plt.plot(np.linspace(0,V_D,300),np.linspace(1,1-nlimgustDslope*V_D,300),'--',color = 'r')


plt.plot(np.linspace(V_B,V_C,300),np.linspace(1+nlimgustBslope*V_B,1+nlimgustCslope*V_C,300),color = 'b')
plt.plot(np.linspace(V_C,V_D,300),np.linspace(1+nlimgustCslope*V_C,1+nlimgustDslope*V_D,300),color = 'b')
plt.plot(np.linspace(V_B,V_C,300),np.linspace(1-nlimgustBslope*V_B,1-nlimgustCslope*V_C,300),color = 'b')
plt.plot(np.linspace(V_C,V_D,300),np.linspace(1-nlimgustCslope*V_C,1-nlimgustDslope*V_D,300),color = 'b')
plt.plot(np.linspace(0,V_B,300),np.linspace(1,1-nlimgustBslope*V_B,300),color = 'b')
plt.plot(V_D*np.ones(100),np.linspace(1+nlimgustDslope*V_D,1-nlimgustDslope*V_D,100),color = 'b')
plt.plot(x_new3,y_new3,'b')


#plt.plot(V_C*np.ones(100),np.linspace(nlimpos,0,100),'--',color = 'g')
#plt.plot(V_A*np.ones(100),np.linspace(nlimpos,0,100),'--',color = 'g')
#plt.plot(V_H*np.ones(100),np.linspace(0,-1,100),'--',color = 'g')
#plt.plot(np.linspace(V_A,V_D,200),nlimpos*np.ones(200),'b')
#plt.plot(V_D*np.ones(100),np.linspace(nlimpos,0,100),'b')
#plt.plot(np.linspace(V_C,V_D,100),np.arange(-1,0,1/100),'b')
#plt.plot(np.linspace(V_H,V_C,100),-1*np.ones(100),'b')
#plt.plot(x_new,y_new,'b')
#plt.plot(x_new2,y_new2,'b')
#plt.text(V_A,nlimpos,'A')
#plt.text(V_H,-1,'H')
#plt.text(V_C,-1,'F')
#plt.text(V_D,0,'E')
#plt.text(V_D,nlimpos,'D')
#plt.text(V_C,nlimpos,'C')

plt.ylim(-2,4)


plt.axhline(y=0, color='k')
plt.axvline(x=0, color='k')

#------------------------------------------------------------------------------------

plt.figure(2)

plt.plot(constlist,nlimpos*np.ones(300),'--',color = 'r')
plt.plot(constlist,-1*np.ones(300),'--', color = 'r')

plt.plot(np.linspace(0,V_B,300),np.linspace(1,1+nlimgustBslope*V_B,300),'--',color = 'r')
plt.plot(np.linspace(0,V_B,300),np.linspace(1,1-nlimgustBslope*V_B,300),'--',color = 'r')
plt.plot(np.linspace(0,V_C,400),np.linspace(1,1+nlimgustCslope*V_C,400),'--',color = 'r')
plt.plot(np.linspace(0,V_C,300),np.linspace(1,1-nlimgustCslope*V_C,300),'--',color = 'r')
plt.plot(np.linspace(0,V_D,400),np.linspace(1,1+nlimgustDslope*V_D,400),'--',color = 'r')
plt.plot(np.linspace(0,V_D,300),np.linspace(1,1-nlimgustDslope*V_D,300),'--',color = 'r')


plt.plot(np.linspace(V_B,V_C,300),np.linspace(1+nlimgustBslope*V_B,1+nlimgustCslope*V_C,300),color = 'b')
plt.plot(np.linspace(V_C,V_D,300),np.linspace(1+nlimgustCslope*V_C,1+nlimgustDslope*V_D,300),color = 'b')
plt.plot(np.linspace(V_B,V_C,300),np.linspace(1-nlimgustBslope*V_B,1-nlimgustCslope*V_C,300),color = 'b')
plt.plot(np.linspace(V_C,V_D,300),np.linspace(1-nlimgustCslope*V_C,1-nlimgustDslope*V_D,300),color = 'b')
plt.plot(np.linspace(0,V_B,300),np.linspace(1,1-nlimgustBslope*V_B,300),color = 'b')
plt.plot(V_D*np.ones(100),np.linspace(1+nlimgustDslope*V_D,1-nlimgustDslope*V_D,100),color = 'b')
plt.plot(x_new3,y_new3,'b')


plt.plot(V_C*np.ones(100),np.linspace(nlimpos,0,100),'--',color = 'g')
plt.plot(V_A*np.ones(100),np.linspace(nlimpos,0,100),'--',color = 'g')
plt.plot(V_H*np.ones(100),np.linspace(0,-1,100),'--',color = 'g')
plt.plot(np.linspace(V_A,V_D,200),nlimpos*np.ones(200),'b')
plt.plot(V_D*np.ones(100),np.linspace(nlimpos,0,100),'b')
plt.plot(np.linspace(V_C,V_D,100),np.arange(-1,0,1/100),'b')
plt.plot(np.linspace(V_H,V_C,100),-1*np.ones(100),'b')
plt.plot(x_new,y_new,'b')
plt.plot(x_new2,y_new2,'b')
plt.text(V_A,nlimpos,'A')
plt.text(V_H,-1,'H')
plt.text(V_C,-1,'F')
plt.text(V_D,0,'E')
plt.text(V_D,nlimpos,'D')
plt.text(V_C,nlimpos,'C')

plt.ylim(-2,4)


plt.axhline(y=0, color='k')
plt.axvline(x=0, color='k')

plt.show()







