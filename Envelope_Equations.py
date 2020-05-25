# Variables

Vs = 0
G_W = 0 
rho = 0
CNmax = 0
S = 0
CLmax = 0 
Vc = 0 
kc = 33
nlimpos = 0
nlim = 4.4
nlimneg = 0

# Determination of Stall Speed

Vs = (2(G_W/S)/rho*CNmax)**(1/2)

# Preiliminary Design Assumption

CLmax = 1.1*CNmax

# Determination of Design Cruise Speed (Vc minimum)

Vc = kc*(G_W/S)**(1/2)

# Determination of Design Limit Load Factor (nlim pos min)

nlimpos = 2.1+24000/(G_W+10000)
#nlimpos<3.8
#nlim=4.4

nlimneg = 0.4*nlim




