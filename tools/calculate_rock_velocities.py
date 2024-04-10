#!/usr/bin/env python3

import pide
import numpy as np
import matplotlib.pyplot as plt

temp = np.array([273,600,1000]) #setting up temperature array
a = pide.pide() #creating the initial object
a.set_temperature(temp)
a.set_pressure([1,2,3])
a.set_solid_phase_method('rock')
a.set_composition_solid_rock(granite = [1.0,0.2,1.0],granulite = [0,0.8,0])
velocities = a.calculate_seismic_velocities()
print(velocities)



