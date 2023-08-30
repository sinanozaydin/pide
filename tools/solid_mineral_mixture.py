#!/usr/bin/env python3

import os,sys

core_path_ext = os.path.join(os.path.dirname(os.path.abspath(__file__)) , '../SEEL')

sys.path.append(core_path_ext)

#ol:olivine
#opx:orthopyroxene
#cpx:clinopyroxene
#garnet:garnet
#mica:mica
#amp:amphibole
#quartz:quartz
#plag:plagioclase
#kfelds:k-feldspar
#sulphide:sulphide
#graphite:graphite
#mixture:mixture
#other:other


import SEEL
import numpy as np
import matplotlib.pyplot as plt

temp = np.arange(600,1300,5) #setting up temperature array
a = SEEL.SEEL() #creating the initial object
a.set_composition_solid_mineral(ol_frac = [0.6],opx_frac = [0.25],cpx_frac = [0.1], garnet_frac = [0.05]) #setting composition
a.set_temperature(temp) #settin temperature array in K
a.set_pressure(1.0)
a.calculate_conductivity(method = 'array')