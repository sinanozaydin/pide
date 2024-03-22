#!/usr/bin/env python3

import os,sys

core_path_ext = os.path.join(os.path.dirname(os.path.abspath(__file__)) , '../pide')

sys.path.append(core_path_ext)

import pide
import numpy as np
import matplotlib.pyplot as plt

temp = np.arange(1073,1473,5)

a = pide.pide() #creating the initial object
a.set_temperature(temp) #settin temperature array in K
a.set_pressure(1)
a.list_mantle_water_solubilities('opx')
a.set_mantle_water_solubility(ol = 1,opx = 3)
a.set_parameter('ti_ol', 0.1)
a.set_parameter('al_opx', 5)
a.calculate_mineral_water_solubility(mineral_name = 'opx', method = 'array')

ax = plt.subplot(111)
ax.plot(temp,a.max_opx_water)

ax.plot(np.array(temp_liu)+273.15,water_liu, 'o')
plt.show()
