# Winbox for bending stress
import input
import Wingbox_design

x_array = Wingbox_design.x_array[2:]
moment_array = Wingbox_design.moment_array
moment_array2 = Wingbox_design.moment_array2

 # coordinates for the airfoil
yield_stress_material = 200 # start of with steel for now

# Read from file

f1 = open('airfoil2.txt', 'r')
lines1 = f1.readlines()

xcoord1 = []
ycoord1 = []
xcoord2 = []
ycoord2 = []

for i in range(0, 103):

    xcoord1.append(float(lines1[i].split()[0]))
    ycoord1.append(float(lines1[i].split()[1]))



for i in range(103, 205):
    xcoord2.append(float(lines1[i].split()[0]))
    ycoord2.append(float(lines1[i].split()[1]))



xcoord2.insert(-1, 1.0)
xcoord2 = xcoord2[::-1]
ycoord2.insert(-1, 0.0)
xcoord2 = xcoord2[::-1]

print(xcoord1)
print(xcoord2)



