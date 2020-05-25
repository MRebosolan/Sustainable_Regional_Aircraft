# Parameter

V_S = 0
G_W = 0 
rho = 0
CNmax = 0
S = 0
CLmax = 0 
V_C = 0 
kc = 33
nlimpos = 0
nlim = 0 #usually 4.4 but check 
nlimneg = 0
V_A = 0
V_B =0 
V = 0 #Quadratic variable
Kg = 0 #Gust alleviation factor
Ude = 0 #depends on altitude in ft (above or below 20k)
mug = 0 
CLa = 0
cbar = 0

# Determination of Stall Speed

V_S = (2(G_W/S)/rho*CNmax)**(1/2)

# Preiliminary Design Assumption

CLmax = 1.1*CNmax

# Determination of Design Cruise Speed (V_C minimum)

V_C = kc*(G_W/S)**(1/2)
V_C = V_B+43#kts

# Determination of Design Limit Load Factor (nlim pos min)

nlimpos = 2.1+24000/(G_W+10000)
#nlimpos<3.8
#nlim=4.4

nlimneg = 0.4*nlim

#Construction of Gust Load Factor Lines

mug = 2*(G_W/S)/(rho*cbar*CLa)
Kg = 0.88*mug/(5.3+mug)
nlim = 1 + (Kg*Ude*V*CLa)/(498*(G_W/S))

# Determination of Design Manoeuvering Speed (V_A minimum)

V_A = V_S*nlim**(1/2)



