#!/usr/bin/env python3

import os,sys

core_path_ext = os.path.join(os.path.dirname(os.path.abspath(__file__)) , '../SEEL')

sys.path.append(core_path_ext)

import SEEL
import numpy as np
import matplotlib.pyplot as plt

temp = np.arange(1200,1800,5)

a = SEEL.SEEL() #creating the initial object
a.set_temperature(temp) #settin temperature array in K
a.set_pressure(1)
a.list_mantle_water_solubilities('ol')
a.set_mantle_water_solubility(ol = 4)
a.set_parameter('ti_ol',0.1)
a.calculate_bulk_mantle_water_solubility(method = 'array')
