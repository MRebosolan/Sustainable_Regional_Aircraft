# -*- coding: utf-8 -*-
"""
Created on Wed Jun  3 13:43:05 2020

@author: malfl
"""
<<<<<<< HEAD

=======
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
>>>>>>> 90088981bcf15a41058c792902db0a3ecc85d8e9
import hydrogen_tank_sizing
#PASSENGER SECTION
#LARGELY BASED ON AIRBUS A220

wseat=18.9*0.0254
warmrest=2*0.0254
waisle=0.5
wclearance=0.01
seat_pitch=32*0.0254
emergency_clearance=10*0.0254

h_aisle=2

hheadroom=1.5
totalabovefloor=0.68
wall_thickness=0.1
floor_thickness=0.1

three=round(wseat*3+wclearance+warmrest*4,3)
two=round(wseat*2+wclearance+warmrest*3,3)

total=round(two+three+waisle,3)
wheadroom=round(total-2*(wclearance+warmrest)-wseat,3)
wfloor=round(total-0.18,3)



#two types I exits
#two types III exits
paxsectionlength=round(15*seat_pitch+emergency_clearance,3)
aft_galley_length=1.5
front_galley_length=2
totalcabinlength=paxsectionlength+aft_galley_length+front_galley_length

outer_diameter=3.486
inner_diameter=3.286
<<<<<<< HEAD
V_tank=30 #CLASS I
R_tank_fus=1.5

<<<<<<< HEAD
tankchoice=3
if tankchoice==1:

    #FIRST OPTION ---- completely aft
    tank_sizing_fuselage(V_tank,R_tank_fus,2)
=======
tankchoice=4
if tankchoice==1:

    #FIRST OPTION ---- completely aft
    thickness,mass, totalmass, diameter,length=tank_sizing_fuselage(V_tank,R_tank_fus,2)
    print('mass: ',totalmass,'diameter: ',diameter,'length: ',length)
>>>>>>> b29fe89e17984d855544738417507b62b9e5c078
        
    
elif tankchoice==2:
    #SECOND OPTION  ---- also in tail
    fractionintail=0.5
    V_tank_tail=V_tank*fractionintail
    V_tank_fus=V_tank-V_tank_tail
    R_tank_tail=1
<<<<<<< HEAD
    tanksectionlength=V_tank_fus/3.14159/R_tank_fus**2
    tanksectionlengthtail=V_tank_tail/3.14159/R_tank_tail**2
    print(tanksectionlengh,tanksectionlengthtail)
=======
    
    tankthicknessfus,structmasstankfus,totstructmasstankfus,tankdiameterfus,lengthfus=tank_sizing_fuselage(V_tank_fus,R_tank_fus,1)
    tankthicknesstail,structmasstanktail,totstructmasstanktail,tankdiametertail,lengthtail=tank_sizing_fuselage(V_tank_tail,R_tank_tail,1)
    print('FUSELAGE ','mass: ',totstructmasstankfus,'diameter: ',tankdiameterfus,'length: ',lengthfus)
    print('TAIL ','mass: ',totstructmasstanktail,'diameter: ',tankdiametertail,'length: ',lengthtail)
    
    
    
    
>>>>>>> b29fe89e17984d855544738417507b62b9e5c078
    
    
    
elif tankchoice==3:
    #THIRD OPTION --- on top
<<<<<<< HEAD
    R_tank_top=tank_sizing(V_tank,totalcabinlength,2)[-1]
    print(R_tank_top)
=======
    tankthickness,masstank,totalmass,DIAMETER_tank_top,length=tank_sizing(V_tank,totalcabinlength,2)
    print('mass: ',totalmass,'diameter: ',DIAMETER_tank_top,'length: ',length)
>>>>>>> b29fe89e17984d855544738417507b62b9e5c078

else:
    #FOURTH OPTION --- combination
    fractioninfus=0.5
    V_tank_fus=fractioninfus*V_tank
    V_tank_top=V_tank-V_tank_fus
<<<<<<< HEAD
    R_tank_top=tank_sizing(V_tank_top,totalcabinlength,1)[-1]
    R_tank_fus=1.5
    tanksectionlength=V_tank_fus/3.14159/R_tank_fus**2

    print(R_tank_top,tanksectionlength)
=======
    
    tankthickness,masstank,totalmass,DIAMETER_tank_top,length=tank_sizing(V_tank_top,totalcabinlength,2)
    print('mass: ',totalmass,'diameter: ',DIAMETER_tank_top,'length: ',length)
    
    R_tank_fus=1.5
    TANK_THICKNESS,STRUCTURAL_TANK_MASS, TOTAL_STRUCTURAL_TANK_MASS, TANK_DIAMETER,LENGTH=tank_sizing_fuselage(V_tank_fus,R_tank_fus,2)
    print(TOTAL_STRUCTURAL_TANK_MASS,TANK_DIAMETER,LENGTH)    

>>>>>>> b29fe89e17984d855544738417507b62b9e5c078
    
=======

V_tank=30 #CLASS I
R_tank_fus=1.5
R_tank_tail=1

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
    
    


fractioninfus=0 #FRACTION OF FUEL IN FUSELAGE (CYLINDER AND TAIL)
fractionintail=0 #FRACTION OF FUSELAGE FUEL IN TAIL
    
V_tank_cyl=(fractioninfus-fractionintail*fractioninfus)*V_tank
V_tank_tail=fractionintail*fractioninfus*V_tank
V_tank_top=V_tank-V_tank_cyl-V_tank_tail
    
#CYLINDER STORAGE
t_cyl,m_cyl, tm_cyl, d_cyl,l_cyl=tank_sizing_fuselage(V_tank_cyl,R_tank_fus,1)
print('CYLINDER TANK: ','| mass: ',tm_cyl,'| diameter: ',d_cyl,'| length: ',l_cyl)

#TAIL STORAGE
t_tail,m_tail, tm_tail, d_tail,l_tail=tank_sizing_fuselage(V_tank_tail,R_tank_tail,1)
print('TAIL TANK: ','| mass: ',tm_tail,'| diameter: ',d_tail,'| length: ',l_tail)

#TOP STORAGE
t_top,m_top,tm_top,d_top,l_top=tank_sizing(V_tank_top,totalcabinlength+l_cyl,2)
print('TOP TANK: ','| mass: ',tm_top,'| diameter: ',d_top,'| length: ',l_top)


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

>>>>>>> 90088981bcf15a41058c792902db0a3ecc85d8e9
    

