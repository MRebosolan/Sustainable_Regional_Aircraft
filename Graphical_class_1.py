# -*- coding: utf-8 -*-
"""
Created on Wed May 27 17:46:12 2020

@author: malfl
"""
import class_1_estimation


mtowlist=[]
xlist=[]
oewlist=[]
kerosenelist=[]
hydrogenlist=[]
tfuellist=[]
tmasslist=[]
hvollist=[]
tdiameterlist=[]

import matplotlib.pyplot as plt

for i in range(0,101):
    outputc1h=CLASS1WEIGHTHYBRID(i/100,1)[0]
    print(CLASS1WEIGHTHYBRID(i/100,1)[1])
    mtowlist.append(outputc1h[0])
    oewlist.append(outputc1h[1])
    kerosenelist.append(outputc1h[5])
    hydrogenlist.append(outputc1h[6])
    tfuellist.append(outputc1h[2])
    hvollist.append(outputc1h[7])
    tdiameterlist.append(outputc1h[8])
    tmasslist.append(outputc1h[9])
    xlist.append(i)

plt.subplot(2,3,1)
plt.plot(xlist,mtowlist,label='MTOW')  
plt.plot(xlist,oewlist,label='OEW')
plt.plot([0,100],[28992,28992],label='MTOW kerosene reserve frac')
plt.plot([0,100],[18273,18273],label='OEW kerosene reserve frac')

plt.ylabel('WEIGHT [kg]')
plt.xlabel('%MASS OF HYDROGEN IN MIXTURE')
plt.legend()

plt.subplot(2,3,2)
plt.plot(xlist,kerosenelist,label='kerosene mass')
plt.plot(xlist,hydrogenlist,label='hydrogen mass')
plt.plot(xlist,tfuellist,label='total fuel mass')
plt.plot([0,100],[660,660],label='kerosene in kerosene reserve frac')
plt.plot([0,100],[2084,2084],label='hydrogen in kerosene reserve frac')
plt.plot([0,100],[2744,2744],label='total fuel in kerosene reserve frac')

plt.ylabel('WEIGHT [kg]')
plt.xlabel('%MASS OF HYDROGEN IN MIXTURE')
plt.legend()
plt.subplot(2,3,3)
plt.plot(xlist,hvollist,label='Hydrogen volume')
plt.ylabel('Volume [CUBIC METERS]')
plt.xlabel('%MASS OF HYDROGEN IN MIXTURE')
plt.legend()

plt.subplot(2,3,4)
plt.plot(xlist,tdiameterlist,label='Tank diameter')
plt.ylabel('Meters')
plt.xlabel('%MASS OF HYDROGEN IN MIXTURE')
plt.legend()
plt.subplot(2,3,5)
plt.plot(xlist,tmasslist,label='Tank mass')
plt.ylabel('kgs')
plt.xlabel('%MASS OF HYDROGEN IN MIXTURE')
plt.legend()
plt.show() 