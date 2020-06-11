import numpy as np

#------------------------------------------------------------------------------------------------------------------
#INPUTS

R_f = 2 # Fuselage radius
t_f = 0.002 #THICkNESS OF FUSELAGE
number_boom = 10 #NUMBER OF BOOMS
moment_cs = 1100000 #MOMENT OF FUSELAGE CROSS SECTION
shear_cs = 150000 #SHEAR FORCE OF THE CROSS SECTION

#------------------------------------------------------------------------------------------------------------------

#make airfoil symmetrical and remove negative values at end
location_angle1 = 360/(2*number_boom)
delta_angle = 360/number_boom

boom_locationx = []
boom_locationy = []
boom_moi = []
boom_area = []
boom_stress = []

for i in range(number_boom):
    boom_locationx.append(R_f * np.cos(np.deg2rad(location_angle1+delta_angle*i)))
    boom_locationy.append(R_f * np.sin(np.deg2rad(location_angle1+delta_angle*i)))

#with these two angle we can find the boom_locationx and the boom_locationy

b = 2*np.pi*R_f/number_boom


for i in range(len(boom_locationx)):
    if i == 0:
        area_boom = t_f * b / 6 + t_f * b * (2 + boom_locationy[i + 1] / boom_locationy[i]) / 6
        moi_boom = area_boom * boom_locationy[i] ** 2
        boom_area.append(area_boom)
        boom_moi.append(moi_boom)
    elif i == len(boom_locationx)-1:
        area_boom = t_f * b * (2 + boom_locationy[i - 1] / boom_locationy[i]) / 6 + t_f * b / 6
        moi_boom = area_boom * boom_locationy[i]  ** 2
        boom_area.append(area_boom)
        boom_moi.append(moi_boom)
    if i > 0 and i < len(boom_locationx)-1:
        area_boom = t_f*b*(2+boom_locationy[i-1]/boom_locationy[i])/6+t_f*b*(2+boom_locationy[i+1]/boom_locationy[i])/6
        moi_boom = area_boom*boom_locationy[i]**2
        boom_area.append(area_boom)
        boom_moi.append(moi_boom)

moi_boom_total = sum(boom_moi)

for i in range(len(boom_locationx)):
        stress_boom = moment_cs*boom_locationy[i] / moi_boom_total
        boom_stress.append(stress_boom)

# print(boom_stress)

# SHEAR FLOW OF THE CROSS SECTION  V_y*Boom_area*boom_locationy/total moi
boom_deltashear = []
for i in range(len(boom_locationx)):
    deltashear_boom = shear_cs*boom_area[i]*boom_locationy[i]/moi_boom_total
    boom_deltashear.append(deltashear_boom)

# print(boom_deltashear)


shear_flow1 = []
shear_flow2 = []
shear_stress = []

shearflow = 0
shear_flow1.append(shearflow)

for i in range(len(boom_locationx[int(number_boom/4):])):
    shearflow = shearflow + boom_deltashear[int(number_boom/4)+i]
    shear_flow1.append(shearflow)

for i in range(len(boom_locationx[:int(number_boom/4)])-1):
    shearflow = shearflow + boom_deltashear[i]
    shear_flow2.append(shearflow)

shear_flow = shear_flow2+shear_flow1
shear_stress = np.array(shear_flow)/t_f


