#!/usr/bin/env python3

import os,sys

import pide
import numpy as np
import matplotlib.pyplot as plt

temp = np.arange(800,1600,5)

a = pide.pide() #creating the initial object
a.set_temperature(temp) #settin temperature array in K
a.set_pressure(1)
a.list_mantle_water_solubilities('ol')
a.set_mantle_water_solubility(ol = 2,opx = 1)
a.set_parameter('ti_ol', 0.1)
a.calculate_mineral_water_solubility(mineral_name = 'ol', method = 'array')


ax = plt.subplot(111)

ax.plot(temp,a.max_ol_water[0])
plt.show()
