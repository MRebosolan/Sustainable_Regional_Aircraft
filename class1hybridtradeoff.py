# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 11:52:53 2020

@author: malfl
"""
from Class_1_estimation import CLASS1WEIGHTHYBRID

mtowlist=[]
xlist=[]
oewlist=[]
kerosenelist=[]
hydrogenlist=[]
tfuellist=[]
tmasslist=[]
hvollist=[]
tdiameterlist=[]
energylist=[]
cjclist=[]
emissionslist=[]
emissionsratiolist=[]


import matplotlib.pyplot as plt

for i in range(99,101):
    outputc1h=CLASS1WEIGHTHYBRID(i/100,1)
    mtowlist.append(outputc1h[0])
    oewlist.append(outputc1h[1])
    kerosenelist.append(outputc1h[5])
    hydrogenlist.append(outputc1h[6])
    cjclist.append(outputc1h[10])
    
    tfuellist.append(outputc1h[2])
    hvollist.append(outputc1h[7])
    tdiameterlist.append(outputc1h[8])
    tmasslist.append(outputc1h[9])
    energylist.append(kerosenelist[-1]*42.8+hydrogenlist[-1]*122.8)
    emissionslist.append(cf(tfuellist[-1], i/100, 1-i/100, input.NOx_H2, input.GWP)[0])
    emissionsratiolist.append(cf(tfuellist[-1], i/100, 1-i/100, input.NOx_H2, input.GWP)[1])
    xlist.append(i)


plt.subplot(2,2,1)
plt.plot(xlist,mtowlist,label='MTOW')
plt.plot(xlist,oewlist,label='OEW')
#plt.plot([0,100],[28992,28992],label='MTOW kerosene reserve frac')
#plt.plot([0,100],[18273,18273],label='OEW kerosene reserve frac')

plt.ylabel('WEIGHT [kg]')
plt.xlabel('%MASS OF HYDROGEN IN MIXTURE')
plt.legend()

plt.subplot(2,2,2)
plt.plot(xlist,kerosenelist,label='kerosene mass')
plt.plot(xlist,hydrogenlist,label='hydrogen mass')
plt.plot(xlist,tfuellist,label='total fuel mass')
#plt.plot([0,100],[660,660],label='kerosene in kerosene reserve frac')
#plt.plot([0,100],[2084,2084],label='hydrogen in kerosene reserve frac')
#plt.plot([0,100],[2744,2744],label='total fuel in kerosene reserve frac')

plt.ylabel('WEIGHT [kg]')
plt.xlabel('%MASS OF HYDROGEN IN MIXTURE')
plt.legend()
plt.subplot(2,2,3)
plt.plot(xlist,hvollist,label='Hydrogen volume')
plt.ylabel('Volume [CUBIC METERS]')
plt.xlabel('%MASS OF HYDROGEN IN MIXTURE')
plt.legend()

plt.subplot(2,2,4)
plt.plot(xlist,cjclist,label='SFC')
plt.ylabel('kg/N/s')
plt.xlabel('%MASS OF HYDROGEN IN MIXTURE')
plt.legend()

#plt.subplot(2,2,1)
#plt.plot(xlist,tmasslist,label='Tank mass')
#plt.ylabel('kgs')
#plt.xlabel('%MASS OF HYDROGEN IN MIXTURE')
#plt.legend()
#
##plt.subplot(3,3,6)
##plt.plot(xlist,energylist,label='total E carried')
##plt.plot(xlist, [i*122.8 for i in hydrogenlist],label='E hydrogen')
##plt.plot(xlist, [i*42.8 for i in kerosenelist],label='E kerosene')
##plt.ylabel('MJ')
##plt.xlabel('%MASS OF HYDROGEN IN MIXTURE')
##plt.legend()
#
#plt.subplot(2,2,2)
#plt.plot(xlist,[i*input.hydrogen_cost for i in hydrogenlist],label='Hydrogen cost')
#plt.plot(xlist,[i*0.6/0.81 for i in kerosenelist],label='kerosene cost')
#plt.plot(xlist,[sum(x) for x in zip([i*input.hydrogen_cost for i in hydrogenlist], [i*0.6/0.81 for i in kerosenelist])],label='total cost')
#plt.ylabel('US DOLLARS')
#plt.xlabel('%MASS OF HYDROGEN IN MIXTURE')
#plt.legend()
#
#plt.subplot(2,2,3)
#plt.plot(xlist,emissionslist,label='emissions')
#plt.xlabel('%MASS OF HYDROGEN IN MIXTURE')
#plt.legend()
#
#plt.subplot(2,2,4)
#plt.plot(xlist,emissionsratiolist,label='CF RATIO wrt CRJ700')
#plt.xlabel('%MASS OF HYDROGEN IN MIXTURE')
#plt.legend()
#plt.show()

