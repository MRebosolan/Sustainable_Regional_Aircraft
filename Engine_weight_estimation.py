import numpy as np
from math import log

#Regression on data obtained from http://www.jet-engine.net/civtfspec.html
#R Squared = 0.954354


thrust = np.array([6500, 7580, 6442, 7500, 7000, 8650, 13790, 14500, 18000, 3700, 4304, 14330])
dry_weights= np.array([1364, 1586, 1586, 1311, 1385, 1625, 2350, 2470, 3800, 754, 852, 2445])

A, B = np.polyfit(np.log(thrust), dry_weights,1)


def engine_weight(dry_thrust_SL, num_engines):
    engine_weight= A * log(dry_thrust_SL) + B
    return engine_weight * num_engines

calculated_thrust = []
for t in thrust:
    calculated_thrust.append(engine_weight(t, 1))

correlation_matrix = np.corrcoef(thrust, calculated_thrust)



