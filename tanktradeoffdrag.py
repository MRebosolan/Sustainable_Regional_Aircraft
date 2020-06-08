# -*- coding: utf-8 -*-
"""
Created on Fri Jun  5 11:31:42 2020

@author: malfl
"""

from cabindesign import cabin_design

finfuslist=[]
cdlist=[]


cd2list=[]
totdraglist=[]
lcyllist=[]
d_podlist=[]
tmf_podlist=[]
cylmass_list=[]


for finfus in range(1,21):
    
    if finfus==1:
        finfus=0.01
    finfus/=20
    

   
    t_cyl,m_cyl, tm_cyl, d_cyl,l_cyl,t_tail,m_tail, tm_tail, d_tail,l_tail\
   ,t_top,m_top,tm_top,d_top,l_top,t_pod,m_pod,tm_pod,d_pod,l_pod,totalcabinlength,V_tank_cyl, V_tank_tail, V_tank_top,V_tank_pod,\
   tm_tanksystem,CGtank,CGfuelfull,CGcomb,totdrag,fuselage_weight,CDzerofus,FFbody,Cfturb,fuselage_area,CDzeropods,fusdrag,poddrag,empennage_length=cabin_design(finfus,0,30,1)  
    cdlist.append(totdrag)
    finfuslist.append(1-finfus)
    cylmass_list.append(tm_cyl+V_tank_cyl*70)
    
    t_cyl,m_cyl, tm_cyl, d_cyl,l_cyl,t_tail,m_tail, tm_tail, d_tail,l_tail\
   ,t_top,m_top,tm_top,d_top,l_top,t_pod,m_pod,tm_pod,d_pod,l_pod,totalcabinlength,V_tank_cyl, V_tank_tail, V_tank_top,V_tank_pod,\
   tm_tanksystem,CGtank,CGfuelfull,CGcomb,totdrag,fuselage_weight,CDzerofus,FFbody,Cfturb,fuselage_area,CDzeropods,fusdrag,poddrag,empennage_length=cabin_design(finfus,0,30,0)  
    cd2list.append(fusdrag)
    lcyllist.append(l_cyl)
    totdraglist.append(totdrag)
    d_podlist.append(d_pod)
    tmf_podlist.append(tm_pod/2+V_tank_pod*70/2)
    if finfus==0:
        podlength=l_pod/2
    
    print(totalcabinlength)
    
    
import matplotlib.pyplot as plt

plt.subplot(231)
plt.plot(finfuslist,cdlist,label='total drag toptank')
plt.plot(finfuslist,cd2list,label='fus drag of pod config')
plt.plot(finfuslist,totdraglist,label='total drag pods')
plt.legend()

plt.subplot(232)
plt.plot(finfuslist,lcyllist,label='length cylindrical tank')
plt.legend()

plt.subplot(233)
plt.plot(finfuslist,d_podlist,label='diameter pod, length is: '+str(podlength))
plt.legend()

plt.subplot(234)
plt.plot(finfuslist,tmf_podlist,label='weight 1 pod')
plt.legend()

plt.subplot(235)
plt.plot(finfuslist,cylmass_list,label='weight cyl')
plt.legend()

plt.show()
    
#PODS WILL ALSO INFLUENCE FLOW OVER WINGS