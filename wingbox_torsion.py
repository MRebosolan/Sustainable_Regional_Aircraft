# Winbox for bending stress
import numpy as np
import input
import Wingbox_design
import matplotlib.pyplot as plt

x_array = Wingbox_design.x_array[2:]
moment_array = Wingbox_design.moment_array
moment_array2 = Wingbox_design.moment_array2

# coordinates for the airfoil
yield_stress_material = 200  # start of with steel for now

# Read from file

f1 = open('airfoil2.txt', 'r')
lines1 = f1.readlines()

xcoord1 = []
ycoord1 = []
ycoord2 = []



for i in range(0, 103):
    xcoord1.append(float(lines1[i].split()[0]))
    ycoord1.append(float(lines1[i].split()[1]))
    ycoord2.append(-1 * float(lines1[i].split()[1]))

xcoord1 = xcoord1[::-1]
ycoord1 = ycoord1[::-1]
ycoord2 = ycoord2[::-1]

# plt.figure(1)
# plt.grid(True, which="major", color="#999999")
# plt.grid(True, which="minor", color="#DDDDDD", ls="--")
# plt.minorticks_on()
# plt.plot(xcoord1, ycoord1, color='r')
# plt.plot(xcoord1, ycoord2, color='r')
# plt.xlim(0, 1)
# plt.ylim(-0.3, 0.3)
# plt.text(0.0, 0.0, 'LE')
# plt.text(1.0, 0.0, 'TE')
# plt.ylabel('y/c [-]')
# plt.xlabel('x/c [-]')
#
# # plt.show()
#------------------------------------------------------------------------------------------------------------------
#INPUTS

chord_length = 2  # chord length in meters
t_d = 0.01 #THICkNESS OF AIRFOIL 10cm
number_points = 20 #NUMBER OF BOOMS  (ON TOP SIDE FOR NOW)
moment_cs = 450000 #MOMENT OF CROSS SECTION

#------------------------------------------------------------------------------------------------------------------

n = 1000 / number_points
boom_locationx = []
boom_locationy = []
boom_area = []
boom_moi = []
stress_boom_upper = []
stress_boom_lower = []

#make airfoil symetrical and remove negative values at end

for i in range(len(xcoord1[:-3])):
    if int(xcoord1[i] * 1000) % n == 0:
        boom_locationx.append(xcoord1[i])
        boom_locationy.append(ycoord1[i])

#compute the boom area and moment of inertia

for i in range(len(boom_locationx)-1):
    if i > 0:
        boom_12dx = (boom_locationx[i - 1] - boom_locationx[i]) * chord_length
        boom_12dy = (boom_locationy[i - 1] - boom_locationy[i]) * chord_length
        boom_23dx = (boom_locationx[i + 1] - boom_locationx[i]) * chord_length
        boom_23dy = (boom_locationy[i + 1] - boom_locationy[i]) * chord_length

        b_1 = np.sqrt((boom_12dx)**2+(boom_12dy)**2)
        b_2 = np.sqrt((boom_23dx)**2+(boom_23dy)**2)

        area_boom = t_d*b_1*(2+boom_locationy[i-1]/boom_locationy[i])/6+t_d*b_2*(2+boom_locationy[i+1]/boom_locationy[i])/6
        moi_boom = area_boom*(boom_locationy[i]*chord_length)**2

        boom_area.append(area_boom)
        boom_moi.append(moi_boom)

# total moment of inertia of the structure (2x top and bottom)
moi_boom_total = sum(2*boom_moi)

for i in range(len(boom_locationx)-1):
    if i>0:
        boom_stress = moment_cs*boom_locationy[i]*chord_length / moi_boom_total
        stress_boom_upper.append(boom_stress)
        stress_boom_lower.append(-boom_stress)

print(stress_boom_upper)
print(stress_boom_lower)

# Find the shear flows

# Finally look at torsion

