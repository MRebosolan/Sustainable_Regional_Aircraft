# -*- coding: utf-8 -*-
"""
Created on Thu Jun  4 11:55:00 2020

@author: malfl
"""

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
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

outer_diameter=3.486
alpha=np.linspace(0,2*np.pi,40)
x=np.cos(alpha)*outer_diameter/2
y=np.sin(alpha)*outer_diameter/2
z=np.outer(np.array([0,totalcabinlength]),np.ones(40))

ztop=np.outer(np.array([0,totalcabinlength+l_cyl]),np.ones(40))
zfustank=np.outer(np.array([totalcabinlength,totalcabinlength+l_cyl]),np.ones(40))
ztail=np.outer(np.array([totalcabinlength+l_cyl,totalcabinlength+l_cyl+l_tail]),np.ones(40))

xtail=np.cos(alpha)*(d_tail/2)
ytail=np.sin(alpha)*(d_tail/2)+outer_diameter/2-d_tail/2

xtop=np.cos(alpha)*(d_top/2+0.01) #Adding a little bit of clearance
ytop=np.sin(alpha)*(d_top/2+0.01)+2

fig = plt.figure()
ax = fig.gca(projection='3d')
ax.set_aspect('equal')         # important!

# ...draw here...
ax.plot_surface(x, y, z, color='b')
ax.plot_surface(xtop, ytop, ztop, color='r')
ax.plot_surface(x, y, zfustank, color='g')
ax.plot_surface(xtail, ytail, ztail, color='b')

ax.set_aspect('equal')

set_axes_equal(ax)             # important!