from numpy import *
import matplotlib.pyplot as plt
from input import *


#Fuel_use = 3515.4  # Total fuel mass at design range
#H2_ff = arange(0., 1.1, 0.1)  # Hydrogen fuel fraction


# GWP
# Global Warming Potential (equivalent emission ratio to CO2 on 100 year scale (CO2-eq))
# For CO2, H20, Nox
# Altitudes 0 to 15 km
GWP_alt = array([[1., 0., -7.1],
                 [1., 0., -7.1],
                 [1., 0., -7.1],
                 [1., 0., -4.3],
                 [1., 0., -1.5],
                 [1., 0., 6.5],
                 [1., 0., 14.5],
                 [1., 0., 37.5],
                 [1., 0., 60.5],
                 [1., 0, 64.7],
                 [1., 0.34, 57.7],
                 [1., 0.43, 46.5],
                 [1., 0.53, 25.6],
                 [1., 0.62, 4.6],
                 [1., 0.72, 0.6]])



# plt.plot(eq, EI_NOxx)
# plt.show()

# ppmNOx = zeros((4, 5))
# for i in range(len(P3)):
#     for j in range(len(fa)):
#         ppmNOx[i, j] = A * (P3[i])**0 * exp(T3/547) * fa[j]**1.6876 * (100*dPP)**(-0.56)
#
# for i in range(0, 4):
#     plt.plot(eq, ppmNOx[i, :])
# plt.show()

def cf(H2_ff, NOx_H2, GWP):
    Fuel_use_H2 = H2_ff * Fuel_use  # Fuel mass H2 at design range
    Fuel_use_ker = (1 - H2_ff) * Fuel_use  # Fuel mass Kerosene at design range

    # kg particles/ kg JET-A1
    perkgJetA1 = [3.16, 1.24, 0.02311]  # CO2, H20, N0x

    perkgH2 = [0, 1, NOx_H2]  # CO2, H20, N0x

    CF_CRJ = 0  # Carbon Footprint CRJ700
    CF_concept = 0  # Carbon Footprint Concept

    # CRJ700 Total emissions (kg CO2-eq)
    for i in range(len(GWP)):
        CO2eq_paxkm = (perkgJetA1[i] * Fuel_use_CRJ * GWP[i]) / (Pax_CRJ * Range_CRJ)
        CF_CRJ += CO2eq_paxkm

    print("CRJ=",CF_CRJ)
    # Concept Total emissions (kg CO2-eq)
    for i in range(len(GWP)):
        CO2eq_paxkm_ker = (perkgJetA1[i] * Fuel_use_ker * GWP[i]) / (Npax * Design_range)
        CO2eq_paxkm_H2 = (perkgH2[i] * Fuel_use_H2 * GWP[i]) / (Npax * Design_range)
        CF_concept += CO2eq_paxkm_ker + CO2eq_paxkm_H2
    print("Concept=", CF_concept)

    Ratio_CF = CF_concept / CF_CRJ  # Ratio of the Carbon Footprints, target = max 0.75

    return Ratio_CF

CF_list = []

for i in range(len(H2_ff)):
   CF = cf(H2_ff[i], NOx_H2, GWP)
   CF_list.append(CF)

# for i in range(1,4):
#    plt.plot(GWP[:, i], GWP[:, 0])

plt.plot(H2_ff, CF_list)
plt.ylabel("CF Concept / CF CRJ")
plt.xlabel("Ratio H2/ Total fuel")
plt.ylim([0, 1])
plt.xlim([0, 1])
plt.grid()
plt.show()
















