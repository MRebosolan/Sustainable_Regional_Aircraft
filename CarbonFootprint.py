from numpy import *
import matplotlib.pyplot as plt

# CRJ700 parameters
Range_CRJ = 2593  # design range
Pax_CRJ = 78  # Number of passengers
Fuel_use_CRJ = 4740  # Fuel mass at design range
Cruise_alt_max_CRJ = 12497  # Max operating altitude

# Our parameters
Range = 2000  # design range
Pax = 75  # Number of passengers
Fuel_use = 3515.4  # Total fuel mass at design range
H2_ff = arange(0., 1.1, 0.1)  # Hydrogen fuel fraction
Cruise_alt = 10000  # Max operating altitude

# kg particles/ kg H2


# fuel to air vs NOx ppm
phi_NOx_H2 = [[0.1, 1.5],
              [0.2, 4.],
              [0.3, 8.],
              [0.4, 10],
              [0.5, 10.6]]


def cf(H2ff, phi_NOx_H2):
    Fuel_use_H2 = H2_ff * Fuel_use  # Fuel mass H2 at design range
    Fuel_use_ker = (1 - H2_ff) * Fuel_use  # Fuel mass Kerosene at design range

    # kg particles/ kg JET-A1
    perkgJetA1 = [3.16, 1.24, 0.02311]  # CO2, H20, N0x

    perkgH2 = [0, 1, phi_NOx_H2]  # CO2, H20, N0x

    # Global Warming Potential (equivalent emission ratio to CO2 on 100 year scale (CO2-eq))
    GWP = [1.0, 0.48, 36.05]  # CO2, H20, N0x

    CF_CRJ = 0  # Carbon Footprint CRJ700
    CF_concept = 0  # Carbon Footprint Concept

    # CRJ700 Total emissions (kg CO2-eq)
    for i in range(len(GWP)):
        CO2eq_paxkm = (perkgJetA1[i] * Fuel_use_CRJ * GWP[i]) / (Pax_CRJ * Range_CRJ)
        CF_CRJ += CO2eq_paxkm

    # Concept Total emissions (kg CO2-eq)
    for i in range(len(GWP)):
        CO2eq_paxkm_ker = (perkgJetA1[i] * Fuel_use_ker * GWP[i]) / (Pax * Range)
        CO2eq_paxkm_H2 = (perkgH2[i] * Fuel_use_H2 * GWP[i]) / (Pax * Range)
        CF_concept += CO2eq_paxkm_ker + CO2eq_paxkm_H2

    Ratio_CF = CF_concept / CF_CRJ  # Ratio of the Carbon Footprints, target = max 0.75

    return Ratio_CF


CF_list = []

for j in range(len(phi_NOx_H2)):

    for i in range(len(H2_ff)):
        row = phi_NOx_H2[j]
        NOx = row[1]
        CF = cf(H2_ff[i], NOx)
        CF_list.append(CF)
        plt.plot(H2_ff, CF)

plt.ylabel("CF Concept / CF CRJ")
plt.xlabel("Ratio H2/ Total fuel")
plt.legend()
plt.grid()
plt.show()















