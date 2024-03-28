#!/usr/bin/env python3

import os,sys

core_path_ext = os.path.join(os.path.dirname(os.path.abspath(__file__)) , '../pide')

sys.path.append(core_path_ext)

import pide
import numpy as np
import matplotlib.pyplot as plt

temp = np.array([273,600,1000]) #setting up temperature array
a = pide.pide() #creating the initial object
a.set_temperature(temp)
a.set_pressure([1,2,3])

a.set_composition_solid_mineral(ol = [0.6,0.6,0.6], opx = [0.4,0.4,0.3], quartz = [0.0,0.0,0.1])
velocities = a.calculate_seismic_velocities()



