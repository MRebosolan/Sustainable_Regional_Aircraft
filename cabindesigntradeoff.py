# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 11:55:00 2020

@author: malfl
"""
import input
import numpy as np
from cabindesign import cabin_design
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from fuselage_weight_estimation import W_fuselage_torenbeek

fig = plt.figure()
ax = fig.add_subplot(221, projection='3d')
#something for plotting in 3D
def set_axes_radius(ax, origin, radius):
    ax.set_xlim3d([origin[0] - radius, origin[0] + radius])
    ax.set_ylim3d([origin[1] - radius, origin[1] + radius])
    ax.set_zlim3d([origin[2] - radius, origin[2] + radius])

def set_axes_equal(ax):
    '''Make axes of 3D plot have equal scale so that spheres appear as spheres,
    cubes as cubes, etc..  This is one possible solution to Matplotlib's
    ax.set_aspect('equal') and ax.axis('equal') not working for 3D.

    Input
      ax: a matplotlib axis, e.g., as output from plt.gca().
    '''

    limits = np.array([
        ax.get_xlim3d(),
        ax.get_ylim3d(),
        ax.get_zlim3d(),
    ])

    origin = np.mean(limits, axis=1)
    radius = 0.5 * np.max(np.abs(limits[:, 1] - limits[:, 0]))
    set_axes_radius(ax, origin, radius)

tm_list=[]
CG_list=[]
fuslist=[]
taillist=[]
outer_diameter=3.486
fuselageweightlist=[]
draglist=[]
#t_cyl,m_cyl, tm_cyl, d_cyl,l_cyl,t_tail,m_tail, tm_tail, d_tail,l_tail,t_top,m_top,tm_top,\
#        d_top,l_top,totalcabinlength,V_tank_cyl, V_tank_tail, V_tank_top,tm_tanksystem,CGtank,CGfuel,CGcomb=cabin_design(0.5,0.5,25)
        
for fractioninfus in range(0,11,3):
    fractioninfus/=10
    for fractionintail in range(0,11,3):
        fractionintail/=10
        
        t_cyl,m_cyl, tm_cyl, d_cyl,l_cyl,t_tail,m_tail, tm_tail, d_tail,l_tail\
           ,t_top,m_top,tm_top,d_top,l_top,t_pod,m_pod,tm_pod,d_pod,l_pod,totalcabinlength,V_tank_cyl, V_tank_tail, V_tank_top,V_tank_pod,\
           tm_tanksystem,CGtank,CGfuelfull,CGcomb,totdrag,fuselage_weight,CDzerofus,FFbody,Cfturb,fuselage_area,CDzeropods,fusdrag,poddrag=cabin_design(fractioninfus,fractionintail,25)
        outer_diameter=3.486
        
        print(fractioninfus,fractionintail)
        tm_list.append(tm_tanksystem)
        CG_list.append(CGcomb)
        fuslist.append(fractioninfus)
        taillist.append(fractionintail)
        draglist.append(totdrag)
        fuselageweightlist.append(fuselageweight)

        



        
x=np.array(fuslist)
y=np.array(taillist)
z=np.array(tm_list)
zskin=np.array(fuselageweightlist)
zCG=np.array(CG_list)
zdrag=np.array(draglist)










#alpha=np.linspace(0,2*np.pi,20)
#x=np.cos(alpha)*outer_diameter/2
#y=np.sin(alpha)*outer_diameter/2
#z=np.outer(np.array([0,totalcabinlength]),np.ones(20))
#
#ztop=np.outer(np.array([0,totalcabinlength+l_cyl]),np.ones(20))
#zfustank=np.outer(np.array([totalcabinlength,totalcabinlength+l_cyl]),np.ones(20))
#ztail=np.outer(np.array([totalcabinlength+l_cyl,totalcabinlength+l_cyl+l_tail]),np.ones(20))
#
#xtail=np.cos(alpha)*(d_tail/2)
#ytail=np.sin(alpha)*(d_tail/2)+outer_diameter/2-d_tail/2
#
#xtop=np.cos(alpha)*(d_top/2+0.01) #Adding a little bit of clearance
#ytop=np.sin(alpha)*(d_top/2+0.01)+1.55+d_top/2
#
#fig = plt.figure()
#ax = fig.gca(projection='3d')
#ax.set_aspect('equal')         # important!

# ...draw here...
ax.scatter(x, y, z, color='b')
ax.set_xlabel('fuselage fraction')
ax.set_ylabel('tail fraction of fuselage')
ax.set_zlabel('total tank mass')

ax = fig.add_subplot(222, projection='3d')
ax.scatter(x, y, zskin, color='b')
ax.set_xlabel('fuselage fraction')
ax.set_ylabel('tail fraction of fuselage')
ax.set_zlabel('fuselage weight')

ax = fig.add_subplot(223, projection='3d')
ax.scatter(x, y, zCG, color='b')
ax.set_xlabel('fuselage fraction')
ax.set_ylabel('tail fraction of fuselage')
ax.set_zlabel('CG')

ax = fig.add_subplot(224, projection='3d')
ax.scatter(x, y, zdrag, color='b')
ax.set_xlabel('fuselage fraction')
ax.set_ylabel('tail fraction of fuselage')
ax.set_zlabel('total drag')

#ax.plot_surface(xtop, ytop, ztop, color='r')
#ax.plot_surface(x, y, zfustank, color='g')
#ax.plot_surface(xtail, ytail, ztail, color='b')
#
#ax.set_aspect('equal')
#
#set_axes_equal(ax)             # important!