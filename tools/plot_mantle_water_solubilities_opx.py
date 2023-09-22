#!/usr/bin/env python3

import os,sys

core_path_ext = os.path.join(os.path.dirname(os.path.abspath(__file__)) , '../SEL')

sys.path.append(core_path_ext)

import SEL
import numpy as np
import matplotlib.pyplot as plt

temp = np.arange(1073,1473,5)

a = SEL.SEL() #creating the initial object
a.set_temperature(temp) #settin temperature array in K
a.set_pressure(1)
a.list_mantle_water_solubilities('opx')
a.set_mantle_water_solubility(ol = 1,opx = 3)
a.set_parameter('ti_ol', 0.1)
a.set_parameter('al_opx', 5)
a.calculate_mineral_water_solubility(mineral_name = 'opx', method = 'array')
print(a.max_opx_water)
ax = plt.subplot(111)
ax.plot(temp,a.max_opx_water[0])
temp_liu = [810.9859154929578,1214.9295774647887]
water_liu = [48.4652358596021,124.38594530143827]
ax.plot(np.array(temp_liu)+273.15,water_liu, 'o')
plt.show()
