# -*- coding: utf-8 -*-
"""
Created on Fri Jun  5 11:31:42 2020

@author: malfl
"""

finfuslist=[]
cdlist=[]

for finfus in range(0,11):
    finfus/=100
    t_cyl,m_cyl, tm_cyl, d_cyl,l_cyl,t_tail,m_tail, tm_tail, d_tail,l_tail\
           ,t_top,m_top,tm_top,d_top,l_top,totalcabinlength,V_tank_cyl, V_tank_tail, V_tank_top,\
           tm_tanksystem,CGtank,CGfuelfull,CGcomb,totdrag,fuselage_weight,CDzerofus,FFbody,Cfturb,fuselage_area=cabin_design(finfus,0,25,1)  
    cdlist.append(CDzerofus)
    finfuslist.append(1-finfus)

import matplotlib.pyplot as plt

plt.plot(finfuslist,cdlist)
plt.show()
    
