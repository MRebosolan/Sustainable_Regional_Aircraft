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
for finfus in range(0,11):
    finfus/=10
    t_cyl,m_cyl, tm_cyl, d_cyl,l_cyl,t_tail,m_tail, tm_tail, d_tail,l_tail\
           ,t_top,m_top,tm_top,d_top,l_top,totalcabinlength,V_tank_cyl, V_tank_tail, V_tank_top,\
           tm_tanksystem,CGtank,CGfuelfull,CGcomb,totdrag,fuselage_weight,CDzerofus,FFbody,Cfturb,fuselage_area,CDzeropods,fusdrag,poddrag=cabin_design(finfus,0,25,1)  
    cdlist.append(totdrag)
    finfuslist.append(1-finfus)
    t_cyl,m_cyl, tm_cyl, d_cyl,l_cyl,t_tail,m_tail, tm_tail, d_tail,l_tail\
           ,t_top,m_top,tm_top,d_top,l_top,totalcabinlength,V_tank_cyl, V_tank_tail, V_tank_top,\
           tm_tanksystem,CGtank,CGfuelfull,CGcomb,totdrag,fuselage_weight,CDzerofus,FFbody,Cfturb,fuselage_area,CDzeropods,fusdrag,poddrag=cabin_design(finfus,0,25,0)  
    cd2list.append(fusdrag)
    lcyllist.append(l_cyl)
    totdraglist.append(totdrag)
    
import matplotlib.pyplot as plt

plt.subplot(121)
plt.plot(finfuslist,cdlist,label='total drag toptank')
plt.plot(finfuslist,cd2list,label='drag fuselage of pods')
plt.plot(finfuslist,totdraglist,label='total drag pods')
plt.legend()
plt.subplot(122)
plt.plot(finfuslist,lcyllist)
plt.show()
    
