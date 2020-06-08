import numpy as np

#------------------------------------------------------------------------------------------------------------------
#INPUTS

R_f = 2 # Fuselage radius
t_f = 0.002 #THICkNESS OF FUSELAGE
number_boom = 8 #NUMBER OF BOOMS
moment_cs = 800000 #MOMENT OF FUSELAGE CROSS SECTION

#------------------------------------------------------------------------------------------------------------------

#make airfoil symetrical and remove negative values at end
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

print(boom_area)
print(moi_boom_total)
print(boom_stress)


